"""End-to-end orchestration: for each sequence file, run codon/GC optimization,
detect hypermutable sites, re-optimize while avoiding them, and write out
per-file CSVs plus a Word comparison report.
"""

from pathlib import Path
from os import path

import numpy as np
import pandas as pd

from eso.detection.common_motifs import load_common_motifs
from eso.detection.dispatch import (
    collapse_recombination_sites,
    collapse_slippage_sites,
    find_recombination_candidates,
    find_slippage_candidates,
    recombination_sites_for_constraints,
    slippage_sites_for_constraints,
)
from eso.detection.methylation import load_motifs, find_motif_sites
from eso.io_utils import file_opener, file_stem, relevant_file_paths, test_input
from eso.optimize import optimization_engine
from eso.report import create_word_document_with_highlighted_differences
from eso.sequence_utils import parse_region


def suspect_site_extractor(target_seq, compute_motifs, num_sites, motifs_path=None,
                            common_motifs=None, recombination_mode='thorough',
                            slippage_mode='default'):
    """Detect recombination and slippage sites (and, if `compute_motifs`, methylation
    motif sites) in `target_seq`. Returns a dict of dataframes keyed by
    'df_recombination', 'df_slippage', and optionally 'df_motifs' - each
    collapsed to one representative per distinct site and limited to
    `num_sites`, for reporting - plus 'df_recombination_raw'/'df_slippage_raw',
    a reduced-but-not-collapsed view (see
    eso.detection.slippage.slippage_sites_for_constraints) meant for building
    correction constraints from instead.

    Detection itself only runs once per category; the report and
    constraint-building views are both derived from that single run. Use the
    _raw dataframes (not the collapsed ones) to build correction constraints:
    collapsing (as the report view does) can silently drop a genuinely
    distinct, only-partially-overlapping hotspot entirely, leaving the region
    it uniquely covered with no fix. `num_sites` is report-only for the same
    reason - limiting it would reintroduce that exact silent-coverage-loss
    risk for whichever candidates it cuts.

    recombination_mode: see eso.detection.dispatch.find_recombination_sites -
        "thorough" (default, Levenshtein-tolerant) or "fast" (exact-match only).
    slippage_mode: see eso.detection.dispatch.find_slippage_sites -
        "default" or "fast" (equivalent sensitivity; "default" is also faster).
    common_motifs: list of str or None
        Names from eso.detection.common_motifs.COMMON_MOTIFS (currently
        "dam", "dcm") to include alongside any `motifs_path` file. At least
        one of `motifs_path`/`common_motifs` is required if `compute_motifs`.
    """
    df_recombination_candidates = find_recombination_candidates(target_seq, mode=recombination_mode)
    df_slippage_candidates = find_slippage_candidates(target_seq, mode=slippage_mode)

    sites_collector = {
        'df_recombination': collapse_recombination_sites(
            df_recombination_candidates, num_sites, mode=recombination_mode),
        'df_slippage': collapse_slippage_sites(df_slippage_candidates, num_sites, mode=slippage_mode),
        'df_recombination_raw': recombination_sites_for_constraints(df_recombination_candidates, mode=recombination_mode),
        'df_slippage_raw': slippage_sites_for_constraints(df_slippage_candidates, mode=slippage_mode),
    }

    if compute_motifs:
        relevant_motifs = list(load_motifs(motifs_path)) if motifs_path else []
        if common_motifs:
            relevant_motifs = relevant_motifs + load_common_motifs(common_motifs)
        sites_collector['df_motifs'] = find_motif_sites(target_seq, num_sites, relevant_motifs)

    return sites_collector


def reoptimize_until_stable(
        curr_seq, compute_motifs, num_sites, motifs_path, common_motifs,
        recombination_mode, slippage_mode, mini_gc, maxi_gc, method, organism_name,
        custom_score_fn, custom_score_minimize, orf_regions=(), exclusion_regions=(),
        max_rounds=5, on_round_start=None):
    """Repeatedly detects hotspots in `curr_seq` and re-optimizes to avoid
    them, continuing for as long as each fresh detection pass still finds
    something.

    A single fixed optimize -> detect -> re-optimize recipe (as this used to
    be, run exactly once) has a real gap: the final re-optimization's own
    synonymous-codon substitutions - chosen to satisfy GC/CAI/avoidance
    objectives - are never re-screened, so a chosen codon can introduce a
    brand-new tandem repeat or near-duplicate that the tool would report as
    clean. This closes that gap by looping detect+re-optimize until a
    detection pass finds nothing at all (the common case: one extra round
    confirms nothing new was introduced) or `max_rounds` is reached (a
    backstop against a pathological sequence where fixing one hotspot always
    creates another, so this can't loop forever).

    `on_round_start(round_index, sites)`, if given, is called right after
    each detection pass, before that round's re-optimization - lets a caller
    report live progress (e.g. eso_desktop's status screen) without this
    function knowing anything about how progress is displayed.

    Returns (final_sequence, obj_description, total_num_edits,
    cumulative_sites, rounds_used):
      - obj_description is the last round's DNAChisel objectives summary, or
        None if no round ever found anything (curr_seq is unchanged from
        what was passed in - callers computing a CAI-after value should fall
        back to their own CAI-before description in that case, since nothing
        was edited).
      - cumulative_sites is shaped like suspect_site_extractor's return
        value, but concatenated across every round that found anything - the
        full history of what was found and corrected, not just the last
        round (which, at convergence, finds nothing by definition).
    """
    empty = {'df_recombination': pd.DataFrame(), 'df_slippage': pd.DataFrame(),
              'df_recombination_raw': pd.DataFrame(), 'df_slippage_raw': pd.DataFrame()}
    if compute_motifs:
        empty['df_motifs'] = pd.DataFrame()
    cumulative = {key: [] for key in empty}

    obj_description = None
    total_num_edits = 0
    rounds_used = 0

    for round_index in range(1, max_rounds + 1):
        sites = suspect_site_extractor(
            curr_seq, compute_motifs, num_sites, motifs_path, common_motifs=common_motifs,
            recombination_mode=recombination_mode, slippage_mode=slippage_mode)
        rounds_used = round_index

        if on_round_start is not None:
            on_round_start(round_index, sites)

        df_motifs = sites.get('df_motifs', pd.DataFrame())
        found_anything = (
            len(sites['df_recombination_raw']) > 0 or len(sites['df_slippage_raw']) > 0
            or len(df_motifs) > 0
        )
        if not found_anything:
            break

        for key in cumulative:
            df = sites.get(key, pd.DataFrame())
            if len(df) > 0:
                cumulative[key].append(df)

        curr_seq, obj_description, num_edits = optimization_engine(
            curr_seq, df_recombination=sites['df_recombination_raw'], df_slippage=sites['df_slippage_raw'],
            df_motifs=df_motifs, mini_gc=mini_gc, maxi_gc=maxi_gc, method=method, organism_name=organism_name,
            custom_score_fn=custom_score_fn, custom_score_minimize=custom_score_minimize,
            orf_regions=orf_regions, exclusion_regions=exclusion_regions)
        total_num_edits += num_edits

    cumulative_sites = {
        key: (pd.concat(dfs, ignore_index=True) if dfs else empty[key].copy())
        for key, dfs in cumulative.items()
    }
    return curr_seq, obj_description, total_num_edits, cumulative_sites, rounds_used


def _extract_cai(objectives_text_summary, num_codons):
    """Parse the CAI objective's score out of DNAChisel's summary text, or
    return None if there wasn't one (e.g. `organism_name` wasn't recognized,
    so no codon-optimization objective was ever added - DNAChisel then
    summarizes as "===> No specifications", which has no ':' to parse and
    used to crash this with an unhandled IndexError).
    """
    first_line = objectives_text_summary.split('\n')[0]
    if ':' not in first_line:
        return None
    cai_score = float(first_line.split(':')[1].strip())
    return np.exp(cai_score / num_codons).round(4)


def backend(data, file, output_path, compute_motifs, num_sites, motifs_path,
            optimize, mini_gc, maxi_gc, method, organism_name, indexes,
            recombination_mode='thorough', slippage_mode='default', common_motifs=None,
            custom_score_fn=None, custom_score_minimize=False):
    """Run the two-pass optimization (CAI/GC only, then + hotspot avoidance) over
    every sequence record in `data`, and write out CSVs + a Word report to
    `output_path/<file_stem>/`.
    """
    recombination_collector = []
    recombination_raw_collector = []
    slippage_collector = []
    slippage_raw_collector = []
    motifs_collector = []
    sequences_for_doc = []

    filename_indexes = file_stem(file[0])
    curr_output_path = path.join(output_path, filename_indexes)
    Path(curr_output_path).mkdir(parents=True, exist_ok=True)

    final_results = []

    for ii, record in enumerate(data):
        curr_seq = str(record.seq).upper()
        original_seq = curr_seq
        seq_indexes = str(ii)

        if (filename_indexes, seq_indexes) not in indexes:
            orf_regions = ()
            exclusion_regions = ()
        else:
            relevant_index_data = indexes[(filename_indexes, seq_indexes)]
            orf_regions = parse_region(relevant_index_data[0])
            exclusion_regions = parse_region(relevant_index_data[1])

        num_codons = sum((orf[1] - orf[0]) / 3 for orf in orf_regions) or len(curr_seq) // 3

        maximal_cai = None
        if optimize:
            curr_seq, obj_description, _ = optimization_engine(
                curr_seq, mini_gc=mini_gc, maxi_gc=maxi_gc, method=method, organism_name=organism_name,
                custom_score_fn=custom_score_fn,
                custom_score_minimize=custom_score_minimize,
                orf_regions=orf_regions, exclusion_regions=exclusion_regions)
            if custom_score_fn is None:
                maximal_cai = _extract_cai(obj_description, num_codons)

        if optimize:
            # Loops detect+re-optimize until a fresh detection pass finds
            # nothing at all (or a round cap is hit) - a single fixed pass
            # can't tell whether its own edits introduced a new hotspot; see
            # reoptimize_until_stable's docstring.
            curr_seq, obj_description_reopt, num_edits, cumulative_sites, _ = reoptimize_until_stable(
                curr_seq, compute_motifs, num_sites, motifs_path, common_motifs,
                recombination_mode, slippage_mode, mini_gc, maxi_gc, method, organism_name,
                custom_score_fn, custom_score_minimize, orf_regions, exclusion_regions)
        else:
            cumulative_sites = suspect_site_extractor(
                curr_seq, compute_motifs, num_sites, motifs_path, common_motifs=common_motifs,
                recombination_mode=recombination_mode, slippage_mode=slippage_mode)
            obj_description_reopt, num_edits = None, 0

        # reporting (CSV) uses the collapsed, num_sites-limited view; the
        # optimizer is given the raw, uncollapsed, unlimited view instead, so a
        # genuinely distinct hotspot that only partially overlaps a
        # higher-scoring one still gets its own correction constraint - see
        # suspect_site_extractor's docstring and docs/detector-comparisons.md.
        # Both views are now the UNION across every re-optimization round, not
        # just the first - see reoptimize_until_stable's docstring.
        df_recombination = cumulative_sites['df_recombination']
        if len(df_recombination) > 0:
            df_recombination.loc[:, 'sequence_number'] = str(ii)
            recombination_collector.append(df_recombination)
        df_recombination_raw = cumulative_sites['df_recombination_raw']
        if len(df_recombination_raw) > 0:
            df_recombination_raw.loc[:, 'sequence_number'] = str(ii)
            recombination_raw_collector.append(df_recombination_raw)

        df_slippage = cumulative_sites['df_slippage']
        if len(df_slippage) > 0:
            df_slippage.loc[:, 'sequence_number'] = str(ii)
            slippage_collector.append(df_slippage)
        df_slippage_raw = cumulative_sites['df_slippage_raw']
        if len(df_slippage_raw) > 0:
            df_slippage_raw.loc[:, 'sequence_number'] = str(ii)
            slippage_raw_collector.append(df_slippage_raw)

        if compute_motifs:
            df_motifs = cumulative_sites['df_motifs']
            if len(df_motifs) > 0:
                df_motifs.loc[:, 'sequence_number'] = str(ii)
                motifs_collector.append(df_motifs)

        if optimize:
            # nothing was found/edited across every round -> curr_seq is
            # exactly pass 1's output, and its CAI is exactly maximal_cai;
            # reuse that instead of inventing an "after" description for a
            # sequence that was never touched a second time.
            obj_description = obj_description_reopt if obj_description_reopt is not None else obj_description

            with open(path.join(curr_output_path, 'final_sequence.txt'), "w", encoding="utf-8") as text_file:
                if custom_score_fn is None and maximal_cai is not None:
                    cai_constrained = _extract_cai(obj_description, num_codons)
                    text_file.write('The maximal CAI of gene (with no constraints) objective:\n')
                    text_file.write(f'{maximal_cai}\n')
                    text_file.write('The CAI of gene (after constraints) objective:\n')
                    text_file.write(f'{cai_constrained}\n')
                elif custom_score_fn is None:
                    text_file.write(
                        "No codon-usage objective was applied (organism_name wasn't recognized).\n")
                else:
                    text_file.write('Optimized using a custom score function instead of CAI/tAI.\n')
                text_file.write('The number of codons edited due to hypermutable site constraints:\n')
                text_file.write(f'{num_edits}\n')
                text_file.write('The final sequence is:\n')
                for line_start in range(0, len(curr_seq), 70):
                    text_file.write(curr_seq[line_start:line_start + 70] + '\n')

        sequences_for_doc.append((f"{filename_indexes}_{ii}", original_seq, curr_seq))
        final_results.append((ii, curr_seq))

    if recombination_collector:
        pd.concat(recombination_collector, ignore_index=True).to_csv(
            path.join(curr_output_path, 'recombination_sites.csv'), index=False)

    if optimize and recombination_raw_collector:
        # every candidate actually given a correction constraint during
        # optimization, not just the collapsed "one representative per
        # distinct site" view above - since those two views can now
        # legitimately differ (see suspect_site_extractor's docstring and
        # docs/detector-comparisons.md), this is what to check against
        # final_sequence.txt if a diff shows an edit with no corresponding
        # row in recombination_sites.csv.
        pd.concat(recombination_raw_collector, ignore_index=True).to_csv(
            path.join(curr_output_path, 'recombination_sites_corrected.csv'), index=False)

    if slippage_collector:
        pd.concat(slippage_collector, ignore_index=True).to_csv(
            path.join(curr_output_path, 'slippage_sites.csv'), index=False)

    if optimize and slippage_raw_collector:
        # see the recombination_sites_corrected.csv comment above.
        pd.concat(slippage_raw_collector, ignore_index=True).to_csv(
            path.join(curr_output_path, 'slippage_sites_corrected.csv'), index=False)

    if compute_motifs and motifs_collector:
        pd.concat(motifs_collector, ignore_index=True).to_csv(
            path.join(curr_output_path, 'motif_sites.csv'), index=False)

    if sequences_for_doc:
        create_word_document_with_highlighted_differences(sequences_for_doc, curr_output_path)

    return final_results


def main(input_folder=None, output_path=None, compute_motifs=False, num_sites=np.inf,
         motifs_path=None, common_motifs=None, optimize=True, mini_gc=0.3, maxi_gc=0.7,
         method='use_best_codon', organism_name='not_specified', indexes=None,
         recombination_mode='thorough', slippage_mode='default', custom_score_fn=None,
         custom_score_minimize=False):
    """Optimize every FASTA/GenBank file in `input_folder`, writing per-file CSVs
    of detected hotspots and the optimized sequence into `output_path`.

    Parameters
    ----------
    input_folder: str
        Directory to scan for .fasta/.fna/.ffn/.faa/.frn/.fa/.gb/.gbk/.genbank
        files (optionally gzipped), directly inside or one level under it.
    output_path: str
        Directory to write results into (one subdirectory per input file).
    compute_motifs: bool
        Whether to also detect methylation motif sites (needs `motifs_path` and/or
        `common_motifs`).
    num_sites: int or float('inf')
        Max number of hotspots to report (in the CSVs/collapsed dataframes)
        per category. Default: all. Does NOT limit how many correction
        constraints are built during optimization - every detected candidate
        is always given a constraint, regardless of `num_sites`, since
        limiting that too could silently leave part of a genuine hotspot
        unconstrained (see eso.pipeline.suspect_site_extractor).
    motifs_path: str or None
        Path to a MEME-minimal-format PSSM file.
    common_motifs: list of str or None
        Names from eso.detection.common_motifs.COMMON_MOTIFS (currently "dam",
        "dcm" - E. coli's methylation systems) to include alongside any
        `motifs_path` file, with no file needed at all. At least one of
        `motifs_path`/`common_motifs` is required if compute_motifs=True.
    optimize: bool
        Whether to codon/GC-optimize and avoid hotspots, vs. just detect them.
    mini_gc, maxi_gc: float in [0, 1]
        Allowed GC-content range within any 50nt window.
    method: {"use_best_codon", "match_codon_usage", "harmonize_rca"}
        Codon optimization strategy.
    organism_name: str
        Host organism for codon optimization (see eso.optimize._codon_optimization_objectives).
    indexes: dict
        Maps (file_stem, seq_index_str) -> (orf_region_string, exclusion_region_string),
        1-indexed and inclusive, e.g. {("my_gene", "0"): ("1-6, 51-68", "1-6, 50-68")}.
        Omit or pass {} to treat entire sequences as the ORF with no exclusions.
    recombination_mode: {"thorough", "fast"}
        See eso.detection.dispatch.find_recombination_sites. "thorough" (default)
        catches near-duplicate hotspots and stays roughly linear (confirmed
        practical up to 1,000,000nt - see docs/detector-comparisons.md);
        "fast" is exact-match only, 19-34x faster, for workloads where that
        speed gap actually matters.
    slippage_mode: {"default", "fast"}
        See eso.detection.dispatch.find_slippage_sites. Both detect identical
        hotspots; "default" is also faster at every length tested.
    custom_score_fn: callable(str) -> float, or None
        If given, replaces CAI/tAI (organism_name/method) with this scoring
        function - see eso.custom_score.CustomScore. Most users should use
        the `--custom-score-file` CLI flag / eso.custom_score.load_custom_score_from_file
        instead of passing a function directly here.
    custom_score_minimize: bool
        See eso.optimize.optimization_engine; only used if custom_score_fn is given.

    Returns
    -------
    (message, results) where message is 'Success!' or a validation error, and
    results is a list of (file, seq_index, optimized_sequence) tuples.
    """
    indexes = indexes or {}
    if output_path is None:
        output_path = path.join(input_folder or '.', 'output')

    files = relevant_file_paths(input_folder=input_folder)
    message = test_input(mini_gc, maxi_gc, indexes, files)

    if message != 'Success!':
        return message, []

    final_results = []
    for file in files:
        data = file_opener(file)
        curr_results = backend(
            data, file, output_path, compute_motifs, num_sites, motifs_path,
            optimize=optimize, mini_gc=mini_gc, maxi_gc=maxi_gc, method=method,
            organism_name=organism_name, indexes=indexes, recombination_mode=recombination_mode,
            slippage_mode=slippage_mode, common_motifs=common_motifs, custom_score_fn=custom_score_fn,
            custom_score_minimize=custom_score_minimize)
        final_results.extend((file, seq_index, seq) for seq_index, seq in curr_results)

    return message, final_results
