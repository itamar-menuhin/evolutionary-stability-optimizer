"""Replication-slippage (SSR - simple sequence repeat) hotspot detection.

Finds short tandem repeats (base units of length 1-15 repeated 3+ times, or a
single nucleotide repeated 4+ times) and scores them with the empirical
mutation-rate formula from the EFM Calculator paper (Jack et al. 2015, ACS
Synthetic Biology, DOI: 10.1021/acssynbio.5b00068).
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
    subseq = subunit * (4 if len(subunit) == 1 else 3)

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
    """Unique substrings of `length` that appear 3+ times in a row, as repeat-unit candidates."""
    list_subunits = [seq[ii:ii + length] for ii in range(len(seq))]
    list_subunits = [subunit for subunit in list_subunits if len(subunit) == length]
    list_subunits = sorted(set(list_subunits))
    return [subunit for subunit in list_subunits if seq.find(subunit * 3) > -1]


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
        df_slippage = df_slippage.head(num_sites)

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
