"""Mode-based routing across independently-developed detector implementations,
so callers choose a thoroughness/speed tradeoff without needing to know which
module originally implemented which algorithm.

Currently covers recombination detection only; slippage and methylation each
have their own primary/alternate split (see eso.detection.slippage vs.
eso.detection.staubility_variant) but aren't wired into dispatch yet.
"""

import numpy as np

from eso.detection import recombination, staubility_variant

RECOMBINATION_MODES = {
    "thorough": recombination.find_recombination_sites,
    "fast": staubility_variant.find_recombination_sites,
}


def find_recombination_sites(seq, num_sites=np.inf, mode="thorough"):
    """Detect recombination (RMD) hotspots, routed to one of two independently
    developed implementations.

    mode="thorough" (default) - eso.detection.recombination: Levenshtein-tolerant,
        catches pairs of sites within edit distance 1 of each other, not just
        exact duplicates. Recommended for gene-length sequences (hundreds to
        low-thousands of nt - ESO's typical case) where near-duplicate hotspots
        matter biologically. Cost scales with sequence length in a way that
        gets expensive for multi-kb sequences (~3.5s for an ~830nt sequence
        with one real repeat, in local benchmarks).

    mode="fast" - eso.detection.staubility_variant: exact 16-mer match only,
        via vectorized n-gram counting. ~100x faster at gene scale. Will miss
        a near-duplicate whenever its point of divergence sits centrally
        enough that no 16-consecutive-nt exact window survives on either side
        (verified: catches a duplicate with a 1nt substitution near either
        edge, since 16+nt of exact match remains; misses the same case when
        the substitution is centered). Recommended for much longer sequences
        (multi-kb constructs, whole plasmids) where the thorough mode's
        per-site cost becomes impractical.
    """
    try:
        detector = RECOMBINATION_MODES[mode]
    except KeyError:
        raise ValueError(f"Unknown recombination mode {mode!r}; choose from {sorted(RECOMBINATION_MODES)}")
    return detector(seq, num_sites)
