"""DNAChisel-based sequence optimization: enforce translation, GC-content
windows, codon usage, and avoidance of detected hypermutable sites.
"""

import warnings

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

# Substring DNAChisel raises when an AvoidPattern's location doesn't align to
# the codon grid under EnforceTranslation: DnaOptimizationProblem.resolve_constraint
# can compute a "jointly mutable" mutation-space region (accounting for codon
# interdependencies) that ends up not overlapping the constraint's own declared
# location at all. AvoidPattern.localized() correctly returns None for that
# (they genuinely don't overlap) - but the DNAChisel caller doesn't check for
# None before calling .evaluate() on it. Confirmed via direct reproduction with
# a non-codon-aligned 2nt repeat and, separately, a 1nt homopolymer run; a
# codon-aligned repeat (e.g. a 3nt unit) never triggers it.
#
# An earlier fix widened every avoid-window outward to whole codons before
# building the AvoidPattern, which does dodge this crash - but it also changes
# what's being asked for: "avoid this pattern at this spot" becomes "avoid this
# pattern ANYWHERE in the whole codon", which for a single-nucleotide pattern
# in a homopolymer run can make the constraint unsatisfiable for every codon in
# the run (e.g. all-Phe TTT/TTC codons always contain a T), silently dropping
# every constraint and leaving the site untouched - no crash, but no effect
# either (confirmed directly: a "ATG"+"T"*12+"TAA" homopolymer produced
# num_edits=0 with that approach). So instead of forcing this by widening the
# window, the retry loop below simply catches the crash like any other
# unsatisfiable-constraint case and drops whichever constraints are actually
# still failing, exactly as it already did for DNAChisel's own
# NoSolutionError-with-no-named-constraint case.
_LOCALIZED_NONE_CRASH_MESSAGE = "'NoneType' object has no attribute 'evaluate'"


def _warn_dropped_constraint(constraint):
    warnings.warn(
        f"Could not satisfy constraint {constraint} - dropping it and continuing without it. "
        "This most often happens when disrupting a detected hypermutable site would require a "
        "change that conflicts with another hard constraint, most commonly translation "
        "preservation (e.g. every synonymous codon at that position still contains the pattern "
        "being avoided, such as a homopolymer landing on a Met/Trp codon with no alternative). "
        "The site this constraint was protecting was left unmodified in the final sequence.",
        stacklevel=3,
    )


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
        and `method` are then ignored. Called once per ORF, on the whole ORF,
        for every trial mutation during optimization - can be slow on long
        sequences or an expensive custom_score_fn (a warning is raised).
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
        # Scoped to each ORF individually, matching how _codon_optimization_objectives
        # already scopes CodonOptimize to `location=orf` - without this, custom scoring
        # would apply across the whole sequence including any non-ORF flanks (UTRs,
        # locked regions), inconsistent with codon-usage scoring's own ORF-only behavior.
        obj = [
            CustomScore(custom_score_fn, location=orf, minimize=custom_score_minimize)
            for orf in orf_regions
        ]
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
    while flag < 60:  # retry, dropping unsatisfiable/crashing constraints, up to 60 times
        problem = dnachisel.DnaOptimizationProblem(sequence=str(seq), constraints=cnst, objectives=obj)
        try:
            problem.resolve_constraints()
            break
        except NoSolutionError as e:
            if e.constraint is not None:
                _warn_dropped_constraint(e.constraint)
                cnst.remove(e.constraint)
                flag += 1
                continue
            drop_and_retry = True
        except AttributeError as e:
            # See _LOCALIZED_NONE_CRASH_MESSAGE above: a genuine DNAChisel-internal
            # crash on non-codon-aligned AvoidPattern locations. Only swallow this
            # exact crash - anything else with this type is a real bug, re-raise it.
            if _LOCALIZED_NONE_CRASH_MESSAGE not in str(e):
                raise
            drop_and_retry = True

        if drop_and_retry:
            # Either DNAChisel's own final consistency check
            # (perform_final_constraints_check, run at the end of
            # resolve_constraints) failed without identifying a single culprit
            # constraint - NoSolutionError.constraint defaults to None, seen in
            # practice with several individually-resolvable but densely packed
            # AvoidPattern constraints that regress each other by the time
            # solving reaches the last one - or resolve_constraint itself
            # crashed on a non-codon-aligned AvoidPattern location. Either way,
            # `cnst.remove(None)` isn't possible here, so instead re-evaluate
            # every constraint ourselves and drop whichever ones are still
            # actually failing, so the retry loop can keep making progress the
            # same way it does for the normal case. Constraints like
            # EnforceGCContent have location=None until
            # initialized_on_problem() fills it in (as a copy, not in-place) -
            # evaluate that initialized copy, not the raw constraint, which
            # would crash the same way on its own unset location.
            failing = [
                c for c in cnst
                if not c.initialized_on_problem(problem, role='constraint').evaluate(problem).passes
            ]
            if not failing:
                raise
            for constraint in failing:
                _warn_dropped_constraint(constraint)
                cnst.remove(constraint)
            flag += 1
    else:
        raise NoSolutionError(f"More than 60 hard constraints were not satisfied ({flag}).", problem=problem)

    problem.optimize()
    obj_description = problem.objectives_text_summary()
    num_edits = problem.number_of_edits()
    final_sequence = str(problem.sequence).upper()

    return final_sequence, obj_description, num_edits
