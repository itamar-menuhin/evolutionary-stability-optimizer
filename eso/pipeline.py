"""End-to-end orchestration: for each sequence file, run codon/GC optimization,
detect hypermutable sites, re-optimize while avoiding them, and write out
per-file CSVs plus a Word comparison report.
"""

from pathlib import Path
from os import path

import numpy as np
import pandas as pd

from eso.detection.common_motifs import load_common_motifs
from eso.detection.dispatch import find_recombination_sites, find_slippage_sites
from eso.detection.methylation import load_motifs, find_motif_sites
from eso.io_utils import file_opener, relevant_file_paths, test_input
from eso.optimize import optimization_engine
from eso.report import create_word_document_with_highlighted_differences
from eso.sequence_utils import parse_region


def suspect_site_extractor(target_seq, compute_motifs, num_sites, motifs_path=None,
                            common_motifs=None, recombination_mode='thorough',
                            slippage_mode='default'):
    """Detect recombination and slippage sites (and, if `compute_motifs`, methylation
    motif sites) in `target_seq`. Returns a dict of dataframes keyed by
    'df_recombination', 'df_slippage', and optionally 'df_motifs'.

    recombination_mode: see eso.detection.dispatch.find_recombination_sites -
        "thorough" (default, Levenshtein-tolerant) or "fast" (exact-match only).
    slippage_mode: see eso.detection.dispatch.find_slippage_sites -
        "default" or "fast" (equivalent sensitivity; "default" is also faster).
    common_motifs: list of str or None
        Names from eso.detection.common_motifs.COMMON_MOTIFS (currently
        "dam", "dcm") to include alongside any `motifs_path` file. At least
        one of `motifs_path`/`common_motifs` is required if `compute_motifs`.
    """
    sites_collector = {
        'df_recombination': find_recombination_sites(target_seq, num_sites, mode=recombination_mode),
        'df_slippage': find_slippage_sites(target_seq, num_sites, mode=slippage_mode),
    }

    if compute_motifs:
        relevant_motifs = list(load_motifs(motifs_path)) if motifs_path else []
        if common_motifs:
            relevant_motifs = relevant_motifs + load_common_motifs(common_motifs)
        sites_collector['df_motifs'] = find_motif_sites(target_seq, num_sites, relevant_motifs)

    return sites_collector


def _extract_cai(objectives_text_summary, num_codons):
    cai_score = float(objectives_text_summary.split('\n')[0].split(':')[1].strip())
    return np.exp(cai_score / num_codons).round(4)


def backend(data, file, output_path, compute_motifs, num_sites, motifs_path,
            optimize, mini_gc, maxi_gc, method, organism_name, indexes,
            recombination_mode='thorough', slippage_mode='default', common_motifs=None,
            custom_score_fn=None, custom_score_window=None, custom_score_minimize=False):
    """Run the two-pass optimization (CAI/GC only, then + hotspot avoidance) over
    every sequence record in `data`, and write out CSVs + a Word report to
    `output_path/<file_stem>/`.
    """
    recombination_collector = []
    slippage_collector = []
    motifs_collector = []
    sequences_for_doc = []

    filename_indexes = path.basename(file[0]).split('.')[0]
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
                custom_score_fn=custom_score_fn, custom_score_window=custom_score_window,
                custom_score_minimize=custom_score_minimize,
                orf_regions=orf_regions, exclusion_regions=exclusion_regions)
            if custom_score_fn is None:
                maximal_cai = _extract_cai(obj_description, num_codons)

        curr_sites_collector = suspect_site_extractor(
            curr_seq, compute_motifs, num_sites, motifs_path, common_motifs=common_motifs,
            recombination_mode=recombination_mode, slippage_mode=slippage_mode)

        df_recombination = curr_sites_collector['df_recombination']
        if len(df_recombination) > 0:
            df_recombination.loc[:, 'sequence_number'] = str(ii)
            recombination_collector.append(df_recombination)

        df_slippage = curr_sites_collector['df_slippage']
        if len(df_slippage) > 0:
            df_slippage.loc[:, 'sequence_number'] = str(ii)
            slippage_collector.append(df_slippage)

        df_motifs = pd.DataFrame()
        if compute_motifs:
            df_motifs = curr_sites_collector['df_motifs']
            if len(df_motifs) > 0:
                df_motifs.loc[:, 'sequence_number'] = str(ii)
                motifs_collector.append(df_motifs)

        if optimize:
            curr_seq, obj_description, num_edits = optimization_engine(
                curr_seq, df_recombination=df_recombination, df_slippage=df_slippage, df_motifs=df_motifs,
                mini_gc=mini_gc, maxi_gc=maxi_gc, method=method, organism_name=organism_name,
                custom_score_fn=custom_score_fn, custom_score_window=custom_score_window,
                custom_score_minimize=custom_score_minimize,
                orf_regions=orf_regions, exclusion_regions=exclusion_regions)

            with open(path.join(curr_output_path, 'final_sequence.txt'), "w") as text_file:
                if custom_score_fn is None:
                    cai_constrained = _extract_cai(obj_description, num_codons)
                    text_file.write('The maximal CAI of gene (with no constraints) objective:\n')
                    text_file.write(f'{maximal_cai}\n')
                    text_file.write('The CAI of gene (after constraints) objective:\n')
                    text_file.write(f'{cai_constrained}\n')
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

    if slippage_collector:
        pd.concat(slippage_collector, ignore_index=True).to_csv(
            path.join(curr_output_path, 'slippage_sites.csv'), index=False)

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
         custom_score_window=None, custom_score_minimize=False):
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
        Max number of hotspots to report/constrain per category. Default: all.
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
    custom_score_window, custom_score_minimize:
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
            custom_score_window=custom_score_window, custom_score_minimize=custom_score_minimize)
        final_results.extend((file, seq_index, seq) for seq_index, seq in curr_results)

    return message, final_results
