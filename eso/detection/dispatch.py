"""Mode-based routing across independently-developed detector implementations,
so callers choose a tradeoff without needing to know which module originally
implemented which algorithm.

Covers recombination and slippage detection; methylation has its own
primary/alternate split too (see eso.detection.methylation vs.
eso.detection.staubility_variant) but isn't wired into dispatch yet.
"""

import numpy as np

from eso.detection import recombination, slippage, staubility_variant

RECOMBINATION_MODES = {
    "thorough": recombination.find_recombination_sites,
    "fast": staubility_variant.find_recombination_sites,
}

SLIPPAGE_MODES = {
    "default": slippage.find_slippage_sites,
    "fast": staubility_variant.find_slippage_sites,
}


def find_recombination_sites(seq, num_sites=np.inf, mode="thorough"):
    """Detect recombination (RMD) hotspots, routed to one of two independently
    developed implementations.

    mode="thorough" (default) - eso.detection.recombination: Levenshtein-tolerant,
        catches pairs of sites within edit distance 1 of each other, not just
        exact duplicates. Recommended for essentially all sequence lengths up
        through tens of kb (~6.5s at 51,300nt in local benchmarks, after
        fixing a pandas-overhead bottleneck - see docs/detector-comparisons.md)
        - there's no reason to pay "fast"'s sensitivity gap until this
        actually becomes inconvenient.

    mode="fast" - eso.detection.staubility_variant: exact 16-mer match only,
        via vectorized n-gram counting. Will miss a near-duplicate whenever
        its point of divergence sits centrally enough that no 16-consecutive-nt
        exact window survives on either side (verified: catches a duplicate
        with a 1nt substitution near either edge, since 16+nt of exact match
        remains; misses the same case when the substitution is centered).
        Reach for this only once "thorough"'s runtime actually matters (very
        large multi-kb+ constructs, or many-sequence batch workloads).
    """
    try:
        detector = RECOMBINATION_MODES[mode]
    except KeyError:
        raise ValueError(f"Unknown recombination mode {mode!r}; choose from {sorted(RECOMBINATION_MODES)}")
    return detector(seq, num_sites)


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
    try:
        detector = SLIPPAGE_MODES[mode]
    except KeyError:
        raise ValueError(f"Unknown slippage mode {mode!r}; choose from {sorted(SLIPPAGE_MODES)}")
    return detector(seq, num_sites)
