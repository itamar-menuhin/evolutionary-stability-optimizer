"""Mode-based routing across independently-developed detector implementations,
so callers choose a tradeoff without needing to know which module originally
implemented which algorithm.

Covers recombination and slippage detection; methylation has only one
implementation (eso.detection.methylation) - a second was built and compared
here, but removed after it was found to disagree in accuracy (not just
speed) with the first - see docs/detector-comparisons.md.
"""

import numpy as np

from eso.detection import recombination, slippage, staubility_variant

# Each category (recombination/slippage) is dispatched at three stages -
# find the raw candidates, collapse them to a report view, and reduce them
# to a constraint-building view (see eso.detection.recombination's module
# docstring split for what each stage means) - across the same two
# implementations ("thorough"/"fast" or "default"/"fast"). That's 6 mode
# tables in total; `_dispatch` below is the one place that turns a mode
# table + a chosen mode into either a result or a friendly ValueError, so a
# bug fix or message tweak to that logic only has to happen once.

RECOMBINATION_MODES = {
    "thorough": recombination.find_recombination_sites,
    "fast": staubility_variant.find_recombination_sites,
}
RECOMBINATION_CANDIDATE_MODES = {
    "thorough": recombination.find_recombination_candidates,
    "fast": staubility_variant.find_recombination_candidates,
}
RECOMBINATION_COLLAPSE_MODES = {
    "thorough": recombination.collapse_recombination_sites,
    "fast": staubility_variant.collapse_recombination_sites,
}
RECOMBINATION_FOR_CONSTRAINTS_MODES = {
    "thorough": recombination.recombination_sites_for_constraints,
    "fast": staubility_variant.recombination_sites_for_constraints,
}

SLIPPAGE_MODES = {
    "default": slippage.find_slippage_sites,
    "fast": staubility_variant.find_slippage_sites,
}
SLIPPAGE_CANDIDATE_MODES = {
    "default": slippage.find_slippage_candidates,
    "fast": staubility_variant.find_slippage_candidates,
}
SLIPPAGE_COLLAPSE_MODES = {
    "default": slippage.collapse_slippage_sites,
    "fast": staubility_variant.collapse_slippage_sites,
}
SLIPPAGE_FOR_CONSTRAINTS_MODES = {
    "default": slippage.slippage_sites_for_constraints,
    "fast": staubility_variant.slippage_sites_for_constraints,
}


def _dispatch(modes, mode, category_label, *args):
    """Look up `mode` in `modes` and call it with `args`, or raise a
    ValueError naming `category_label` ("recombination"/"slippage") and
    listing the modes actually available - the one error-message format
    every public function below shares.
    """
    try:
        implementation = modes[mode]
    except KeyError:
        raise ValueError(f"Unknown {category_label} mode {mode!r}; choose from {sorted(modes)}")
    return implementation(*args)


def find_recombination_sites(seq, num_sites=np.inf, mode="thorough"):
    """Detect recombination (RMD) hotspots, routed to one of two independently
    developed implementations.

    mode="thorough" (default) - eso.detection.recombination: Levenshtein-tolerant,
        catches pairs of sites within edit distance 1 of each other, not just
        exact duplicates. Benchmarked as roughly linear from 51,400nt through
        1,000,000nt (~126.7s at 1,000,000nt, in local benchmarks after fixing
        a pandas-overhead bottleneck - see docs/detector-comparisons.md), with
        no breakdown point found at any tested scale - recommended by default
        at essentially any realistic sequence length.

    mode="fast" - eso.detection.staubility_variant: exact 16-mer match only,
        via vectorized n-gram counting. Will miss a near-duplicate whenever
        its point of divergence sits centrally enough that no 16-consecutive-nt
        exact window survives on either side (verified: catches a duplicate
        with a 1nt substitution near either edge, since 16+nt of exact match
        remains; misses the same case when the substitution is centered).
        19-34x faster than "thorough" at every length tested - reach for this
        only when that speed gap itself matters (e.g. many-sequence batch
        workloads), not because "thorough" becomes intractable.
    """
    return _dispatch(RECOMBINATION_MODES, mode, "recombination", seq, num_sites)


def find_slippage_sites(seq, num_sites=np.inf, mode="default"):
    """Detect slippage (SSR) hotspots, routed to one of two independently
    developed implementations.

    Unlike recombination, both modes detect exactly the same hotspots -
    verified via 300 randomized trials with zero sensitivity or row-count
    mismatches after fixing bugs in both implementations (see
    docs/detector-comparisons.md). This is purely a speed choice, and
    "default" wins it outright:

    mode="default" - eso.detection.slippage: after fixing an O(n^2) candidate-
        scan (see docs/detector-comparisons.md), this is faster than "fast"
        at every length tested, from a few hundred nt through 300,000nt,
        with the gap widening as length grows (10x faster at 300kb).

    mode="fast" - eso.detection.staubility_variant: kept as an independent
        second implementation (useful as a cross-check, and it's a distinct
        algorithm, not just a slower copy) - but there is no longer a length
        range where it's actually faster than "default".
    """
    return _dispatch(SLIPPAGE_MODES, mode, "slippage", seq, num_sites)


def find_recombination_candidates(seq, mode="thorough"):
    """Every candidate recombination (RMD) site-pair, WITHOUT collapsing
    overlapping pairs down to one representative per real hotspot - routed to
    one of the two implementations documented on find_recombination_sites.

    This is the shared base find_recombination_sites (report view, via
    collapse_recombination_sites) and recombination_sites_for_constraints
    (constraint-building view) are both derived from - use whichever of
    those two fits, not this function directly, unless you specifically want
    every raw candidate with nothing reduced at all.
    """
    return _dispatch(RECOMBINATION_CANDIDATE_MODES, mode, "recombination", seq)


def find_slippage_candidates(seq, mode="default"):
    """Every candidate slippage (SSR) hotspot, WITHOUT collapsing overlapping
    candidates down to one representative per physical site - routed to one
    of the two implementations documented on find_slippage_sites.

    This is the shared base find_slippage_sites (report view, via
    collapse_slippage_sites) and slippage_sites_for_constraints
    (constraint-building view) are both derived from - use whichever of
    those two fits, not this function directly, unless you specifically want
    every raw candidate with nothing reduced at all.
    """
    return _dispatch(SLIPPAGE_CANDIDATE_MODES, mode, "slippage", seq)


def collapse_recombination_sites(df_pairs, num_sites=np.inf, mode="thorough"):
    """Collapse a raw recombination-candidates dataframe (from
    find_recombination_candidates, same `mode`) down to the human-facing
    "distinct sites" view - one representative pair per real hotspot, limited
    to `num_sites`. Lets a caller run detection once and derive both the raw
    (constraint-building) and collapsed (report) views without re-detecting.
    """
    return _dispatch(RECOMBINATION_COLLAPSE_MODES, mode, "recombination", df_pairs, num_sites)


def collapse_slippage_sites(df_slippage, num_sites=np.inf, mode="default"):
    """Collapse a raw slippage-candidates dataframe (from
    find_slippage_candidates, same `mode`) down to the human-facing "distinct
    sites" view - one representative per physical site, limited to
    `num_sites`. Lets a caller run detection once and derive both the raw
    (constraint-building) and collapsed (report) views without re-detecting.
    """
    return _dispatch(SLIPPAGE_COLLAPSE_MODES, mode, "slippage", df_slippage, num_sites)


def recombination_sites_for_constraints(df_pairs, mode="thorough"):
    """Reduce a raw recombination-candidates dataframe (from
    find_recombination_candidates, same `mode`) for feeding into
    eso.constraints.recombination_to_multiple_avoidance_sites, WITHOUT the
    coverage-loss risk collapse_recombination_sites has - see
    eso.detection.recombination.recombination_sites_for_constraints.
    """
    return _dispatch(RECOMBINATION_FOR_CONSTRAINTS_MODES, mode, "recombination", df_pairs)


def slippage_sites_for_constraints(df_slippage, mode="default"):
    """Reduce a raw slippage-candidates dataframe (from
    find_slippage_candidates, same `mode`) for feeding into
    eso.detection.slippage.modify_df_slippage, WITHOUT the coverage-loss risk
    collapse_slippage_sites has - see
    eso.detection.slippage.slippage_sites_for_constraints.
    """
    return _dispatch(SLIPPAGE_FOR_CONSTRAINTS_MODES, mode, "slippage", df_slippage)
