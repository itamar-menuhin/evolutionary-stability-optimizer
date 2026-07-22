"""DNAChisel-based sequence optimization: enforce translation, GC-content
windows, codon usage, and avoidance of detected hypermutable sites.
"""

import pandas as pd
import dnachisel
from dnachisel.DnaOptimizationProblem import NoSolutionError

from eso.codon_usage import CODON_USAGE_TABLES
from eso.constraints import (
    convert_df_to_constraints,
    exclusion_site_correcter,
    recombination_to_multiple_avoidance_sites,
)
from eso.custom_score import CustomScore
from eso.detection.slippage import modify_df_slippage


def _codon_optimization_objectives(organism_name, orf_regions, method):
    """Build DNAChisel CodonOptimize objectives.

    `organism_name` may be one of eso.codon_usage.CODON_USAGE_TABLES' custom
    tables (not present in the Kazusa database), or any species/TaxID
    supported by python-codon-tables.
    """
    if organism_name in CODON_USAGE_TABLES:
        codon_usage_table = CODON_USAGE_TABLES[organism_name]()
        return [
            dnachisel.CodonOptimize(location=orf, codon_usage_table=codon_usage_table.copy(), method=method)
            for orf in orf_regions
        ]

    try:
        return [dnachisel.CodonOptimize(species=organism_name, location=orf, method=method) for orf in orf_regions]
    except Exception:
        print('unknown organism, no codon optimization!')
        return []


def optimization_engine(
    seq,
    mini_gc=0.3,
    maxi_gc=0.7,
    window_size_gc=50,
    method='use_best_codon',
    organism_name='not_specified',
    custom_score_fn=None,
    custom_score_window=None,
    custom_score_minimize=False,
    df_recombination=None,
    df_slippage=None,
    df_motifs=None,
    orf_regions=(),
    exclusion_regions=(),
):
    """Optimize `seq` for codon usage, GC content, and (if hotspot dataframes
    are given) avoidance of detected recombination/slippage/methylation sites,
    while preserving the amino-acid translation over `orf_regions`.

    Citation: Jack, Leonard, Mishler, Renda, Leon, Suarez & Barrick (2015).
    "Predicting the genetic stability of engineered DNA sequences with the
    EFM Calculator." ACS Synthetic Biology. DOI: 10.1021/acssynbio.5b00068

    Parameters
    ----------
    seq: str
        DNA sequence (ACGT alphabet).
    mini_gc, maxi_gc: float in [0, 1]
        Allowed GC-content range within any `window_size_gc`-nt window.
    method: {"use_best_codon", "match_codon_usage", "harmonize_rca"}
        Codon optimization strategy (see DNAChisel's CodonOptimize).
    organism_name: str
        Host organism for codon optimization: one of
        eso.codon_usage.CODON_USAGE_TABLES' keys, a python-codon-tables
        species name/TaxID, or "not_specified" to skip codon optimization.
    custom_score_fn: callable(str) -> float, or None
        If given, replaces the CodonOptimize (CAI/tAI-style) objective with
        eso.custom_score.CustomScore wrapping this function (higher is
        better sequence, unless custom_score_minimize=True). `organism_name`
        and `method` are then ignored. See eso.custom_score.CustomScore for
        the windowed-vs-global tradeoff.
    custom_score_window: int or None
        If given, custom_score_fn is called on each successive
        `custom_score_window`-nt chunk and the results summed (fast, but
        only correct if the score truly decomposes that way). If None
        (default), custom_score_fn is called once on the whole sequence on
        every trial mutation during optimization (always correct, but can be
        very slow on long sequences - a warning is raised when this mode is
        used).
    custom_score_minimize: bool
        If True, treats a lower custom_score_fn value as better.
    df_recombination, df_slippage, df_motifs: pandas.DataFrame or None
        Detected hotspots to avoid (from eso.detection.*); pass None or an
        empty dataframe to skip a given constraint type.
    orf_regions: sequence of (start, end) tuples
        Regions to keep in-frame and translation-preserving. Defaults to the
        whole sequence, trimmed to a multiple of 3.
    exclusion_regions: sequence of (start, end) tuples
        Regions that must not be modified.

    Returns
    -------
    (final_sequence, objectives_summary, num_edits)
    """
    if df_recombination is None:
        df_recombination = pd.DataFrame()
    if df_slippage is None:
        df_slippage = pd.DataFrame()
    if df_motifs is None:
        df_motifs = pd.DataFrame()

    if len(orf_regions) == 0:
        new_last_index = ((len(seq) - 1) // 3) * 3
        orf_regions = [(0, new_last_index)]

    if custom_score_fn is not None:
        obj = [CustomScore(custom_score_fn, window=custom_score_window, minimize=custom_score_minimize)]
    else:
        obj = _codon_optimization_objectives(organism_name, orf_regions, method)

    cnst = [dnachisel.EnforceGCContent(mini=mini_gc, maxi=maxi_gc, window=window_size_gc)]
    for orf in orf_regions:
        cnst.append(dnachisel.EnforceTranslation(location=orf))

    cnst.extend(dnachisel.AvoidChanges(location=region) for region in exclusion_regions)

    if not df_recombination.empty:
        df_rec = recombination_to_multiple_avoidance_sites(df_recombination, exclusion_regions)
        df_rec = exclusion_site_correcter(df_rec, exclusion_regions)
        cnst.extend(convert_df_to_constraints(df_rec))

    if not df_slippage.empty:
        df_slip = modify_df_slippage(df_slippage)
        df_slip = exclusion_site_correcter(df_slip.copy(), exclusion_regions)
        cnst.extend(convert_df_to_constraints(df_slip))

    if not df_motifs.empty:
        df_mot = df_motifs.copy()[['start_index', 'end_index', 'actual_site']].rename(
            columns={'start_index': 'start', 'end_index': 'end', 'actual_site': 'sequence'})
        df_mot.loc[:, 'start'] = df_mot['start'].astype(int)
        # eso.detection.methylation's end_index is INCLUSIVE (the index of the
        # motif's last nucleotide), but everything downstream of here -
        # exclusion_site_correcter, convert_df_to_constraints, and DNAChisel's
        # own Location - uses EXCLUSIVE end (matching Python slicing). Without
        # the +1, every motif's AvoidPattern location was exactly one
        # nucleotide too short to ever contain its own pattern, so
        # DNAChisel could never find it there and always reported the
        # constraint as trivially satisfied - methylation-motif avoidance
        # silently did nothing. Confirmed directly: before this fix, a GATC
        # motif passed as df_motifs survived optimization completely
        # untouched (0 edits); see tests/test_optimize.py.
        df_mot.loc[:, 'end'] = df_mot['end'].astype(int) + 1
        df_mot = exclusion_site_correcter(df_mot, exclusion_regions)
        cnst.extend(convert_df_to_constraints(df_mot))

    problem = None
    flag = 0
    while flag < 60:  # retry, dropping unsatisfiable constraints, up to 60 times
        problem = dnachisel.DnaOptimizationProblem(sequence=str(seq), constraints=cnst, objectives=obj)
        try:
            problem.resolve_constraints()
            break
        except NoSolutionError as e:
            cnst.remove(e.constraint)
            flag += 1
    else:
        raise NoSolutionError(f"More than 60 hard constraints were not satisfied ({flag}).", problem=problem)

    problem.optimize()
    obj_description = problem.objectives_text_summary()
    num_edits = problem.number_of_edits()
    final_sequence = str(problem.sequence).upper()

    return final_sequence, obj_description, num_edits
