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
from eso.sequence_utils import validate_dna_alphabet

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

#: Constraint types representing a hard, user-specified requirement - not a
#: "nice to avoid if possible" hotspot/enzyme-avoidance pattern - that the
#: retry-and-drop logic below must never itself remove, even when DNAChisel's
#: solver can't satisfy everything at once. EnforceTranslation isn't listed
#: here: confirmed directly against DNAChisel's own source
#: (dnachisel.EnforceTranslation.enforced_by_nucleotide_restrictions is True)
#: that it's enforced by restricting the mutation space itself - only
#: synonymous edits are ever reachable in the first place - so it can never
#: actually be named as the failing constraint to begin with, unlike
#: EnforceGCContent (ordinary post-hoc resolution, confirmed via the same
#: source check) - the constraint this list was written for. Deliberately NOT
#: also protected here via the ambiguous-fallback's blanket evaluate() (which
#: isn't guarded by that same structural restriction, unlike the named-
#: culprit path): confirmed directly (see the regression this caused in
#: tests/test_optimize.py's retry-budget-at-scale test) that EnforceTranslation
#: can transiently evaluate as failing mid-crash-recovery on a densely
#: conflicting sequence even though the reading frame itself is never actually
#: at risk there, and treating it as protected in that path forces this loop's
#: much slower one-at-a-time drop strategy (below) onto an ordinary dense
#: hotspot-vs-hotspot conflict that doesn't need it - actually reachable, not
#: just theoretical, and it burns through the retry budget for no benefit.
_PROTECTED_CONSTRAINT_TYPES = (dnachisel.EnforceGCContent,)

#: Above this many dropped constraints, the shrink/verify pass (see the end
#: of the retry loop below) is skipped entirely, falling back to trusting
#: the greedy drop order as final - the same behavior this verification was
#: added on top of. Restoring each dropped constraint costs one full
#: resolve_constraints() call, and that cost compounds with the retry
#: loop's own already-quadratic worst case (many rounds, each re-evaluating
#: every remaining constraint, when a protected constraint like
#: EnforceGCContent stays implicated for a long stretch): benchmarked
#: directly on a synthetic near-0%-GC poly-T stress sequence, going from
#: ~125 to ~250 dropped constraints raised total runtime from ~10s to ~45s,
#: and ~250 to ~475 reached ~180s - a real, compounding cost, not a flat
#: overhead, that a genuinely dense conflict can hit in practice. Below this
#: cap, minimality is verified directly rather than just aimed for by the
#: greedy order; above it, correctness rigor is deliberately traded for
#: practical runtime - some risk of an unnecessarily-dropped low-severity
#: constraint remains in that regime, same as before this verification
#: existed at all.
_SHRINK_VERIFICATION_CAP = 100


def _constraints_to_drop(cnst, problem):
    """Decide which constraint(s) in `cnst` to drop to make progress on a
    currently-unresolved `problem`. Returns a (possibly empty) list.

    A constraint that isn't currently failing is never a candidate at all, so
    a low-risk site that isn't actually part of the present conflict is never
    touched just because it happens to be the least severe thing in `cnst`.

    If none of the currently-failing constraints is a hard, protected
    requirement (see _PROTECTED_CONSTRAINT_TYPES) - e.g. several mutually-
    conflicting hotspot/enzyme-avoidance constraints crashing DNAChisel's
    solver, none of them a hard requirement - every currently-failing
    constraint is returned at once, exactly like the original fallback: there
    is no hard requirement to protect here, so there's no reason to spend one
    retry-budget round per constraint re-discovering the same failures one at
    a time (this matters in practice - a single dense/repetitive sequence can
    have well over 100 such constraints crash/conflict at once).

    If a protected constraint IS among what's currently failing, only the
    SINGLE least-severe (ascending own risk score - see
    eso.constraints.convert_df_to_constraints) currently-failing, non-protected
    constraint is returned - retrying the whole solve after each single
    removal, so nothing is dropped unless it's still actually part of the
    conflict once re-evaluated, and the least risky candidate is always tried
    before a riskier one. An empty list here (no droppable constraint is
    currently failing alongside the protected one) means a genuine,
    irreducible conflict: even removing every non-protected constraint can't
    rescue the protected one.
    """
    failing = [
        c for c in cnst
        if not c.initialized_on_problem(problem, role='constraint').evaluate(problem).passes
    ]
    if not failing:
        return []
    protected_failing = [c for c in failing if isinstance(c, _PROTECTED_CONSTRAINT_TYPES)]
    droppable_failing = [c for c in failing if not isinstance(c, _PROTECTED_CONSTRAINT_TYPES)]
    if not protected_failing:
        return droppable_failing
    if not droppable_failing:
        return []
    return [min(droppable_failing, key=lambda c: getattr(c, 'eso_severity', float('inf')))]


def _try_restore(cnst_without, constraint, obj, seq):
    """Attempt a single resolve_constraints() with `constraint` added back to
    `cnst_without` (a fresh list, not mutated). Returns the resulting
    resolved problem if it succeeds, or None if it doesn't (either a
    NoSolutionError, or the same localized-None crash the main retry loop
    below already knows how to recognize and treat as "still doesn't work" -
    any OTHER exception is a real bug and propagates, matching the main
    loop's own philosophy). This is deliberately a single attempt, not a
    nested retry loop - a "best-effort, provably safe" check of whether this
    one constraint specifically is still necessary given everything else
    currently in play, not a full re-optimization.
    """
    candidate_cnst = cnst_without + [constraint]
    trial_problem = dnachisel.DnaOptimizationProblem(sequence=str(seq), constraints=candidate_cnst, objectives=obj)
    try:
        trial_problem.resolve_constraints()
    except NoSolutionError:
        return None
    except AttributeError as e:
        if _LOCALIZED_NONE_CRASH_MESSAGE not in str(e):
            raise
        return None
    return trial_problem


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


#: sentinel meaning "skip codon-usage optimization" - the documented
#: default. Deliberately distinct from an actually-unrecognized organism
#: name below: skipping via this sentinel is normal, expected behavior and
#: shouldn't warn; an unrecognized organism name is a likely user mistake
#: (a typo, or the wrong TaxID) that silently produces the exact same
#: "no codon optimization happened" outcome, so it's worth surfacing.
NOT_SPECIFIED = 'not_specified'


def _codon_optimization_objectives(organism_name, orf_regions, method, codon_usage_table=None):
    """Build DNAChisel CodonOptimize objectives.

    `organism_name` may be one of eso.codon_usage.CODON_USAGE_TABLES' custom
    tables (not present in the Kazusa database), NOT_SPECIFIED to skip codon
    optimization entirely, or any species/TaxID supported by
    python-codon-tables. Ignored entirely if `codon_usage_table` is given -
    see optimization_engine's docstring.
    """
    if codon_usage_table is not None:
        return [
            dnachisel.CodonOptimize(location=orf, codon_usage_table=codon_usage_table.copy(), method=method)
            for orf in orf_regions
        ]

    if organism_name == NOT_SPECIFIED:
        return []

    if organism_name in CODON_USAGE_TABLES:
        codon_usage_table = CODON_USAGE_TABLES[organism_name]()
        return [
            dnachisel.CodonOptimize(location=orf, codon_usage_table=codon_usage_table.copy(), method=method)
            for orf in orf_regions
        ]

    try:
        return [dnachisel.CodonOptimize(species=organism_name, location=orf, method=method) for orf in orf_regions]
    except Exception:
        # Genuinely unrecognized (not the NOT_SPECIFIED sentinel above) -
        # likely a typo or wrong TaxID, silently producing the same
        # no-codon-optimization outcome as an intentional skip. A real
        # warning (not a print) so a caller can detect this programmatically
        # (warnings.catch_warnings()) instead of it only ever reaching a
        # terminal a human happens to be watching.
        warnings.warn(
            f"'{organism_name}' isn't a recognized organism/species name or TaxID, and isn't one of "
            "eso.codon_usage.CODON_USAGE_TABLES' bundled tables - codon-usage optimization was skipped "
            "entirely for this run, exactly as if organism_name had been left as its default "
            f"({NOT_SPECIFIED!r}). Check the spelling, or try the TaxID instead of the species name.",
            stacklevel=3,
        )
        return []


def optimization_engine(
    seq,
    mini_gc=0.3,
    maxi_gc=0.7,
    window_size_gc=50,
    method='use_best_codon',
    organism_name='not_specified',
    codon_usage_table=None,
    custom_score_fn=None,
    custom_score_minimize=False,
    df_recombination=None,
    df_slippage=None,
    df_motifs=None,
    orf_regions=(),
    exclusion_regions=(),
    avoid_hairpins=False,
    hairpin_stem_size=20,
    hairpin_window=200,
    avoid_enzymes=(),
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
        Ignored if `codon_usage_table` is given.
    codon_usage_table: dict, or None
        A codon-usage table of the form {'*': {'TAA': 0.33, ...}, 'K': {...},
        ...} - see eso.codon_usage.load_custom_codon_table_from_file for
        loading one from a CSV file (most users should use the
        `--codon-usage-table-file` CLI flag / that function instead of
        building this dict by hand). If given, this is used directly for
        codon-usage optimization instead of `organism_name` (both
        eso.codon_usage.CODON_USAGE_TABLES lookups and python-codon-tables
        species lookups are skipped entirely). Ignored if custom_score_fn is
        also given, same as `organism_name`.
    custom_score_fn: callable(str) -> float, or None
        If given, replaces the CodonOptimize (CAI/tAI-style) objective with
        eso.custom_score.CustomScore wrapping this function (higher is
        better sequence, unless custom_score_minimize=True). `organism_name`,
        `codon_usage_table`, and `method` are then ignored. Called once per ORF, on the whole ORF,
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
    avoid_hairpins: bool
        If True, adds DNAChisel's built-in AvoidHairpins constraint (IDT's
        secondary-structure guideline: avoid a stem whose reverse complement
        recurs nearby) over the whole sequence - a standard gene-design
        concern (strong secondary structure, especially over a start
        codon/RBS, can suppress translation initiation) this library didn't
        previously expose at all despite it being one call away in an
        already-installed dependency.
    hairpin_stem_size, hairpin_window: int
        See dnachisel.AvoidHairpins; only used if avoid_hairpins=True.
    avoid_enzymes: sequence of str
        Restriction enzyme names (e.g. "BsaI", "EcoRI" - see
        dnachisel.list_common_enzymes() for the full recognized list) whose
        recognition sites should be avoided - useful when the optimized
        sequence needs to be cloned into a specific vector/backbone. Also
        DNAChisel functionality this library didn't previously expose.

    Returns
    -------
    (final_sequence, objectives_summary, num_edits)
    """
    # Nothing downstream validates this - confirmed directly that an IUPAC
    # ambiguity code (e.g. "N", "R") crashes deep inside DNAChisel with an
    # unhandled TranslationError or KeyError depending on which constraint
    # reaches it first, with no eso-level message at all. Checked here, at
    # the lowest-level entry point, so both the file-based CLI/pipeline path
    # (which also checks earlier, in eso.io_utils.test_input, for a
    # friendlier pre-flight message) and any direct API caller are covered.
    validate_dna_alphabet(seq)

    if df_recombination is None:
        df_recombination = pd.DataFrame()
    if df_slippage is None:
        df_slippage = pd.DataFrame()
    if df_motifs is None:
        df_motifs = pd.DataFrame()

    if len(orf_regions) == 0:
        new_last_index = (len(seq) // 3) * 3
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
        obj = _codon_optimization_objectives(organism_name, orf_regions, method, codon_usage_table=codon_usage_table)

    cnst = [dnachisel.EnforceGCContent(mini=mini_gc, maxi=maxi_gc, window=window_size_gc)]
    for orf in orf_regions:
        cnst.append(dnachisel.EnforceTranslation(location=orf))

    cnst.extend(dnachisel.AvoidChanges(location=region) for region in exclusion_regions)

    if avoid_hairpins:
        cnst.append(dnachisel.AvoidHairpins(stem_size=hairpin_stem_size, hairpin_window=hairpin_window))

    for enzyme_name in avoid_enzymes:
        try:
            enzyme_pattern = dnachisel.EnzymeSitePattern(enzyme_name)
        except KeyError:
            raise ValueError(
                f"'{enzyme_name}' isn't a recognized restriction enzyme name. See "
                "dnachisel.list_common_enzymes() for the full list of recognized names "
                "(case-sensitive, e.g. 'BsaI', 'EcoRI')."
            ) from None
        cnst.append(dnachisel.AvoidPattern(enzyme_pattern))

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
    # 60 was chosen empirically back when constraint counts were always small
    # (num_sites capped how many hotspot-avoidance constraints existed at
    # all). Since num_sites became report-only and constraint-building
    # started using every raw, uncollapsed candidate (see
    # docs/detector-comparisons.md's overlap-collapse coverage-gap entry), a
    # single dense/repetitive sequence can legitimately produce well over 60
    # individual AvoidPattern constraints - a fixed 60-round budget could then
    # raise NoSolutionError purely from running out of retries, not from any
    # actual unresolvable conflict. Scale the budget with how many
    # constraints there actually are (computed once, before any get dropped,
    # so it doesn't shrink alongside `cnst` as rounds progress), keeping 60 as
    # the floor for the common, small-constraint-count case this was
    # originally tuned for.
    retry_budget = max(60, len(cnst))
    dropped = []  # every constraint removed below, across both branches - see the shrink pass after this loop
    while flag < retry_budget:  # retry, dropping unsatisfiable/crashing constraints
        problem = dnachisel.DnaOptimizationProblem(sequence=str(seq), constraints=cnst, objectives=obj)
        try:
            problem.resolve_constraints()
            break
        except NoSolutionError as e:
            named_culprit = e.constraint
        except AttributeError as e:
            # See _LOCALIZED_NONE_CRASH_MESSAGE above: a genuine DNAChisel-internal
            # crash on non-codon-aligned AvoidPattern locations. Only swallow this
            # exact crash - anything else with this type is a real bug, re-raise it.
            if _LOCALIZED_NONE_CRASH_MESSAGE not in str(e):
                raise
            named_culprit = None

        # The common, cheap case: DNAChisel itself named exactly one failing,
        # droppable (not hard/protected) constraint - just drop that one, no
        # need to re-evaluate everything else ourselves.
        if named_culprit is not None and not isinstance(named_culprit, _PROTECTED_CONSTRAINT_TYPES):
            cnst.remove(named_culprit)
            dropped.append(named_culprit)
            flag += 1
            continue

        # Either DNAChisel named a *protected* hard constraint (most commonly
        # EnforceGCContent) as the one it couldn't satisfy, or its failure was
        # ambiguous (perform_final_constraints_check failed with no single
        # culprit - NoSolutionError.constraint defaults to None - or the
        # localized-None crash above). Either way, re-evaluate every
        # constraint ourselves and let _constraints_to_drop decide what to
        # remove: if a protected constraint isn't actually in conflict here,
        # this is the same batch-drop-everything-failing behavior as before
        # (cheap, and fine - none of what's failing is a hard requirement);
        # if one is, only the single least-severe non-protected constraint
        # still actually failing is removed per round, so nothing gets
        # dropped unless it's still part of the conflict once re-evaluated.
        # An empty result means a genuine, irreducible conflict: a hard,
        # user-specified requirement (currently just GC content) can't be
        # satisfied together with the sequence's other hard requirements even
        # after dropping everything else - eso will not silently drop it and
        # ship a sequence outside the range requested.
        # Constraints like EnforceGCContent have location=None until
        # initialized_on_problem() fills it in (as a copy, not in-place) -
        # evaluate that initialized copy, not the raw constraint, which would
        # crash the same way on its own unset location.
        to_drop = _constraints_to_drop(cnst, problem)
        if not to_drop:
            raise NoSolutionError(
                "Could not satisfy EnforceGCContent together with the sequence's other requested "
                "constraints, and there was no other, less-critical constraint (hotspot avoidance, "
                "restriction-enzyme sites, etc.) left to drop instead that would have helped. Unlike "
                "those, this is a hard, user-specified requirement - eso will not silently drop it "
                "and ship a sequence outside the range you asked for. Try widening mini_gc/maxi_gc, "
                "or relaxing whichever other constraint this one is actually in conflict with.",
                problem=None,
            )
        for constraint in to_drop:
            cnst.remove(constraint)
        dropped.extend(to_drop)
        flag += 1
    else:
        raise NoSolutionError(
            f"More than {retry_budget} hard constraints were not satisfied ({flag}).", problem=problem)

    # Verify minimality: the greedy walk above drops constraints as it goes,
    # in ascending-severity order, but never re-checks whether an EARLIER
    # drop is still needed once LATER drops have also happened - so a
    # constraint can end up dropped even though, given everything else that
    # ended up dropped too, it was never actually necessary on its own. This
    # is exactly the "only remove what actually helps resolve the clash"
    # guarantee that matters, not just an ordering preference - so it's
    # verified directly here, not just aimed for during the walk above,
    # UNLESS there are too many dropped constraints to make that practical
    # (see _SHRINK_VERIFICATION_CAP - a real, compounding runtime cost on a
    # genuinely dense conflict, not just a flat overhead).
    if len(dropped) <= _SHRINK_VERIFICATION_CAP:
        # Try restoring each dropped constraint, most severe (most worth
        # keeping) first; whatever restores cleanly (a single
        # resolve_constraints() with it added back - see _try_restore) stays
        # restored, and only a constraint that genuinely still breaks the
        # resolution when restored is left dropped and warned about.
        for constraint in sorted(dropped, key=lambda c: getattr(c, 'eso_severity', float('inf')), reverse=True):
            restored_problem = _try_restore(cnst, constraint, obj, seq)
            if restored_problem is not None:
                cnst.append(constraint)
                problem = restored_problem
            else:
                _warn_dropped_constraint(constraint)
    else:
        for constraint in dropped:
            _warn_dropped_constraint(constraint)

    problem.optimize()
    obj_description = problem.objectives_text_summary()
    num_edits = problem.number_of_edits()
    final_sequence = str(problem.sequence).upper()

    return final_sequence, obj_description, num_edits
