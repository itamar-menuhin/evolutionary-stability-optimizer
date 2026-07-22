"""Replication-slippage (SSR - simple sequence repeat) hotspot detection.

Finds short tandem repeats (base units of length 1-15 repeated 3+ times, or a
single nucleotide repeated 12+ times) and scores them with the empirical
mutation-rate formula from the EFM Calculator paper (Jack et al. 2015, ACS
Synthetic Biology, DOI: 10.1021/acssynbio.5b00068).

The length-1 (homopolymer) minimum is 12, not 3, for two stacked reasons:

1. log10_prob = -12.9 + 0.729*n crosses the -9 filter cutoff (applied below)
   at n=6 - runs of n=4/5 always fail it and are never worth detecting.
2. Every homopolymer run of n>=6 is *also* always detected as a length-2 site
   (any 6+ identical characters trivially contain "XX" repeated 3+ times),
   scored as -4.749 + 0.063*floor(n/2) - and collapse_overlapping_intervals
   below only keeps the single highest-scoring representation per physical
   site. Comparing the two formulas: the length-1 score grows ~0.729/nt while
   the length-2 score grows only ~0.0315/nt, so length-2 always wins for
   n=6..11 (e.g. n=8: length-1 -7.068 vs length-2 -4.497) - the length-1
   detection is real but its result is always discarded. The crossover is
   n=12 exactly (length-1 -4.152 vs length-2 -4.371) - only from there does
   the length-1 reading start being the one that survives, correctly
   reflecting a homopolymer run's much higher risk than a linear-in-repeats
   length-2 score would suggest.

So detecting length-1 sites below n=12 is provably pure waste - the result
never changes the final output regardless, since length-2 detection (run
independently, over the same sequence) always produces a same-or-better-scoring
overlapping candidate for that region. This doesn't apply to length>1 units:
log10_prob = -4.749 + 0.063*n is already > -9 at the smallest detectable n=3,
so every length>1 candidate always passes the -9 filter on its own.
"""

import numpy as np
import pandas as pd

from eso.detection._overlap import collapse_overlapping_intervals

SLIPPAGE_COLUMNS = ['start', 'end', 'length_base_unit', 'sequence']


def _find_all(a_str, sub):
    """Non-overlapping occurrences of `sub` in `a_str`."""
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1:
            return
        yield start
        start += len(sub)


def _find_longest_match(seq, subunit, start_index):
    """Given a start location and a repeat unit, find the full extent of the run."""
    curr_seq = seq[start_index:]
    len_sub = len(subunit)
    end_shift = 0
    while curr_seq[:len_sub] == subunit:
        end_shift += len_sub
        curr_seq = curr_seq[len_sub:]
    end_index = start_index + end_shift
    return start_index, end_index, seq[start_index:end_index], len_sub


def _generate_slippage_sites_current_subunit(seq, subunit):
    curr_slippage_sites = []
    # 12 for length-1 units (below that, length-2 detection of the same
    # region always outscores it - see module docstring), 3 for longer units.
    subseq = subunit * (12 if len(subunit) == 1 else 3)

    indexes_curr = list(_find_all(seq, subseq))
    if not indexes_curr:
        return curr_slippage_sites

    curr_slippage_sites.append(_find_longest_match(seq, subunit, indexes_curr[0]))

    for index in indexes_curr[1:]:
        end_previous_site = curr_slippage_sites[-1][1]
        if index <= end_previous_site:
            continue
        curr_slippage_sites.append(_find_longest_match(seq, subunit, index))

    return curr_slippage_sites


def _find_slippage_len1(seq):
    slippage_sites = []
    for nt in ['A', 'C', 'G', 'T']:
        slippage_sites.extend(_generate_slippage_sites_current_subunit(seq, nt))
    return pd.DataFrame.from_records(
        data=slippage_sites, columns=['start', 'end', 'sequence', 'length_base_unit'],
    ).drop_duplicates().reset_index(drop=True)


def _find_relevant_subunits_len_l(seq, length):
    """Unique substrings of `length` that appear 3+ times in a row, as repeat-unit candidates.

    Profiling showed the original approach - build every unique substring,
    then call `seq.find(subunit*3)` once per candidate - spends ~93% of total
    find_slippage_sites runtime in `str.find` (88,151 calls on a 9.6kb random
    sequence), since each call rescans the whole sequence: O(unique
    candidates) x O(n) = O(n^2) for a sequence with mostly-unique substrings.
    Building a position dict in one pass and checking adjacency via O(1) set
    membership (positions p, p+length, p+2*length all present means
    `subunit*3` occurs at p) is equivalent but O(n) instead.
    """
    positions_by_subunit = {}
    for ii in range(len(seq) - length + 1):
        positions_by_subunit.setdefault(seq[ii:ii + length], []).append(ii)

    relevant = []
    for subunit, positions in positions_by_subunit.items():
        position_set = set(positions)
        if any((p + length) in position_set and (p + 2 * length) in position_set for p in positions):
            relevant.append(subunit)

    return sorted(relevant)


def _find_slippage_len_l(seq, length):
    """Note: 1 < length < 16."""
    relevant_subunits = _find_relevant_subunits_len_l(seq, length)
    slippage_sites = []
    for subunit in relevant_subunits:
        slippage_sites.extend(_generate_slippage_sites_current_subunit(seq, subunit))
    return pd.DataFrame.from_records(data=slippage_sites, columns=['start', 'end', 'sequence', 'length_base_unit'])


def find_slippage_sites(seq, num_sites=np.inf):
    """Find candidate slippage (SSR) hotspots in `seq` across repeat-unit lengths 1-15.

    Returns a dataframe sorted by descending mutation risk
    (`log10_prob_slippage_ecoli`), limited to `num_sites` sites if given.
    """
    slippage_sites_list = [_find_slippage_len1(seq)]
    for length in range(2, 16):
        slippage_sites_list.append(_find_slippage_len_l(seq, length))
    df_slippage = pd.concat(slippage_sites_list, ignore_index=True)[SLIPPAGE_COLUMNS]

    df_slippage.loc[:, 'num_base_units'] = (
        df_slippage.sequence.apply(len) / df_slippage.length_base_unit
    ).astype(int)

    df_slippage.loc[:, 'log10_prob_slippage_ecoli'] = -4.749 + 0.063 * df_slippage['num_base_units']
    df_slippage.loc[df_slippage.length_base_unit == 1, 'log10_prob_slippage_ecoli'] = (
        -12.9 + 0.729 * df_slippage['num_base_units']
    )

    df_slippage = df_slippage[df_slippage.log10_prob_slippage_ecoli > -9]

    # different base-unit lengths (and different phase offsets within the same
    # length) can each detect the same physical repeat as a separate row -
    # e.g. "GCGCGCGC" is a valid length-2 run starting at position N, AND its
    # 1-shifted substring "CGCGCG" is a separate valid length-2 run starting
    # at N+1. Collapsing by exact 'start' alone (the original approach) misses
    # this, since the rows don't share a start position. Keep one
    # representative (highest scoring) per group of overlapping ranges instead.
    df_slippage = collapse_overlapping_intervals(df_slippage, score_col='log10_prob_slippage_ecoli')
    df_slippage = df_slippage.sort_values(['log10_prob_slippage_ecoli', 'length_base_unit'], ascending=[False, False])

    if num_sites < np.inf:
        df_slippage = df_slippage.head(int(num_sites))

    return df_slippage


def modify_df_slippage(df_slippage):
    """Convert each slippage-site row (N repeated base units) into ~N/2 rows of
    individual base units to avoid (skipping every other one) - enough to
    disrupt the repeat without necessarily eliminating it entirely.

    Example: a row for "TGTGTGTGTG" (base unit "TG", num_base_units=5) becomes
    3 rows, each a single "TG" occurrence and its coordinates.
    """
    df_list = []

    for idx in df_slippage.index:
        num_base_units = df_slippage.loc[idx].num_base_units
        length = df_slippage.loc[idx].length_base_unit
        for i in range(0, (num_base_units - 1), 2):
            df_list.append({
                'sequence': df_slippage.loc[idx].sequence[int(i * length):int((i + 1) * length)],
                'start': int(df_slippage.loc[idx].start + i * length),
                'end': int(df_slippage.loc[idx].start + (i + 1) * length),
            })

    return pd.DataFrame(df_list)
