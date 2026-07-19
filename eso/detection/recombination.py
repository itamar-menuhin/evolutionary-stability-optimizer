"""Repeat-mediated deletion (RMD) / recombination hotspot detection.

Finds pairs of near-identical sites (length >= 16, Levenshtein distance <= 1,
non-overlapping) that are candidates for recombination-mediated deletion, and
scores them with the empirical mutation-rate formula from the EFM Calculator
paper (Jack et al. 2015, ACS Synthetic Biology, DOI: 10.1021/acssynbio.5b00068).
"""

import numpy as np
import pandas as pd
from Levenshtein import distance as levenshtein_distance

from eso.detection._overlap import ranges_overlap
from eso.sequence_utils import add_backward_sites, shorten_sequences

# matches the columns actually produced by the non-empty path below
# (sequence_1/sequence_2, not a single 'sequence' column)
RECOMBINATION_COLUMNS = [
    'sequence_1', 'start_1', 'end_1', 'sequence_2', 'start_2', 'end_2',
    'location_delta', 'site_length', 'log10_prob_recombination_ecoli',
]


def _generate_all_recombination_sites(seq):
    """Generate candidate 16-17mers (forward + reverse complement) at every offset,
    plus their single-insertion/deletion/substitution-shortened variants, so that
    near-duplicate site pairs can be found by exact match instead of all-pairs comparison.
    """
    forward_insertions = [(seq[ii:ii + 17], ii, ii + 16) for ii in range(len(seq))]
    forward_insertions = [x for x in forward_insertions if len(x[0]) == 17]

    df_forward_insertions = pd.DataFrame.from_records(
        data=forward_insertions, columns=['sequence', 'start', 'end'])

    df_insertions = add_backward_sites(df_forward_insertions)
    df_first_sites = shorten_sequences(df_forward_insertions)
    df_substitutions = add_backward_sites(df_first_sites)

    df_deletions = shorten_sequences(df_first_sites)
    df_deletions = add_backward_sites(df_deletions)

    return df_first_sites, df_insertions, df_substitutions, df_deletions


def _generate_neighbors(curr_seq):
    """All sequences up to one insertion, deletion, or substitution away from curr_seq."""
    neighbors_substitutions = set()
    neighbors_deletions = set()
    neighbors_insertions = set()

    for ii in range(len(curr_seq)):
        prefix = curr_seq[:ii]
        suffix = curr_seq[ii:]
        for nt in ['A', 'C', 'G', 'T']:
            neighbors_insertions.add(prefix + nt + suffix)
        suffix = suffix[1:]
        for nt in ['A', 'C', 'G', 'T']:
            neighbors_substitutions.add(prefix + nt + suffix)
        neighbors_deletions.add(prefix + suffix)

    return neighbors_substitutions, neighbors_deletions, neighbors_insertions


def _order_site12(df):
    lower_second = (df.start_1 > df.start_2)
    df_output = df.copy()
    for col in ['sequence', 'start', 'end']:
        df_output.loc[lower_second, f'{col}_1'] = df[f'{col}_2']
        df_output.loc[lower_second, f'{col}_2'] = df[f'{col}_1']
    return df_output


def _build_sequence_index(df):
    """Map sequence -> list of (start, end) tuples, for O(1) dict lookup instead
    of a `df.sequence.isin(...)` filter (which rescans/reindexes the whole
    frame) called once per candidate site.
    """
    index = {}
    for sequence, start, end in zip(df.sequence, df.start, df.end):
        index.setdefault(sequence, []).append((start, end))
    return index


def _generate_relevant_pairs(df_first_sites, df_insertions, df_substitutions, df_deletions):
    """Profiling showed this function dominates find_recombination_sites'
    runtime (5.1s of 7.0s on an ~830nt sequence), almost entirely in pandas
    overhead: `.isin()` plus the surrounding boolean-mask filtering and
    `.iloc`/`.loc` indexing, called once per candidate site (potentially
    O(n) sites). Building a plain-dict index once up front, and iterating
    the input frames as plain Python tuples via zip() rather than
    `.iloc[ii, col]`, avoids nearly all of that per-call pandas overhead
    while producing identical pairs.
    """
    pairs = []

    substitutions_index = _build_sequence_index(df_substitutions)
    deletions_index = _build_sequence_index(df_deletions)
    insertions_index = _build_sequence_index(df_insertions)

    for sequence_1, start_1, end_1 in zip(df_first_sites.sequence, df_first_sites.start, df_first_sites.end):
        neighbors_substitutions, neighbors_deletions, neighbors_insertions = _generate_neighbors(sequence_1)

        for neighbors, index in (
            (neighbors_substitutions, substitutions_index),
            (neighbors_deletions, deletions_index),
            (neighbors_insertions, insertions_index),
        ):
            for neighbor_seq in neighbors:
                for start_2, end_2 in index.get(neighbor_seq, ()):
                    pairs.append((sequence_1, start_1, end_1, neighbor_seq, start_2, end_2))

    df_pairs = pd.DataFrame.from_records(
        pairs, columns=['sequence_1', 'start_1', 'end_1', 'sequence_2', 'start_2', 'end_2'])

    df_pairs = _order_site12(df_pairs)

    # drop overlapping pairs and duplicates
    df_pairs = df_pairs[df_pairs.end_1 < df_pairs.start_2].drop_duplicates().reset_index(drop=True)

    return df_pairs


def _elongate_sites(row, full_seq):
    """Extend each site pair outward while they remain within Levenshtein distance 1
    and non-overlapping, so the reported site is the full recombination hotspot
    rather than just its detected 16-mer core.
    """
    sequence_1, start_1, end_1 = row.sequence_1, row.start_1, row.end_1 + 1
    sequence_2, start_2, end_2 = row.sequence_2, row.start_2, row.end_2 + 1

    while (levenshtein_distance(full_seq[start_1:end_1 + 1], full_seq[start_2:end_2 + 1], score_cutoff=1) < 2) \
            and (end_1 + 1 < start_2):
        end_1 += 1
        end_2 += 1

    while (levenshtein_distance(full_seq[start_1 - 1:end_1], full_seq[start_2 - 1:end_2], score_cutoff=1) < 2) \
            and (end_1 + 1 < start_2):
        start_1 -= 1
        start_2 -= 1

    row.sequence_1, row.start_1, row.end_1 = full_seq[start_1:end_1], start_1, end_1
    row.sequence_2, row.start_2, row.end_2 = full_seq[start_2:end_2], start_2, end_2
    return row


def calc_recombination_score(location_delta, site_length):
    """Empirical log10 probability of recombination-mediated deletion (EFM Calculator paper)."""
    a, b, c, alpha = 5.8, 1465.6, 0, 29

    first_component = a + location_delta
    second_component = -1 * (alpha / site_length)
    third_component = site_length / (1 + b * site_length + c * location_delta)

    recombination_probability = (first_component ** second_component) * third_component
    return np.log10(recombination_probability)


def _collapse_overlapping_pairs(df_pairs):
    """Different seed 16-mers for the same real hotspot converge, via
    elongation, to slightly different (start, end) extents rather than one
    canonical pair - so even exact-coordinate dedup leaves several
    near-identical rows per real site. Collapse them via non-max suppression:
    walk pairs in descending score order, keeping a pair only if BOTH its
    site_1 and site_2 ranges are still free of an already-kept pair's
    corresponding range. Requiring overlap on both sides (not just one) avoids
    merging two genuinely distinct hotspots that happen to share one site.
    """
    kept_rows = []
    kept_ranges = []  # list of ((start_1,end_1), (start_2,end_2))

    for _, row in df_pairs.sort_values('log10_prob_recombination_ecoli', ascending=False).iterrows():
        range_1 = (row.start_1, row.end_1)
        range_2 = (row.start_2, row.end_2)
        if any(ranges_overlap(range_1, r1) and ranges_overlap(range_2, r2) for r1, r2 in kept_ranges):
            continue
        kept_rows.append(row)
        kept_ranges.append((range_1, range_2))

    return pd.DataFrame(kept_rows, columns=df_pairs.columns) if kept_rows else df_pairs.iloc[0:0]


def find_recombination_sites(seq, num_sites=np.inf):
    """Find candidate recombination (RMD) hotspots in `seq`.

    A pair of sites is a candidate if both are >=16nt, within Levenshtein
    distance 1 of each other, and non-overlapping. Returns a dataframe sorted
    by descending mutation risk (`log10_prob_recombination_ecoli`), limited to
    `num_sites` distinct site-pairs if given.
    """
    empty_df = pd.DataFrame(columns=RECOMBINATION_COLUMNS)

    df_first_sites, df_insertions, df_substitutions, df_deletions = _generate_all_recombination_sites(seq)
    df_pairs = _generate_relevant_pairs(df_first_sites, df_insertions, df_substitutions, df_deletions)

    if df_pairs.shape[0] == 0:
        return empty_df

    df_pairs = df_pairs.apply(lambda row: _elongate_sites(row, seq), axis=1)

    df_pairs.loc[:, 'location_delta'] = df_pairs.start_2 - df_pairs.end_1
    df_pairs.loc[:, 'site_length'] = df_pairs.apply(
        lambda row: max(row.end_1 - row.start_1, row.end_2 - row.start_2), axis=1)

    df_pairs.loc[:, 'log10_prob_recombination_ecoli'] = df_pairs.apply(
        lambda x: calc_recombination_score(x.location_delta, x.site_length), axis=1)

    df_pairs = df_pairs[df_pairs.log10_prob_recombination_ecoli > -9]

    # collapse the many candidate-seed rows that converge (via elongation) on
    # the same or overlapping real hotspot down to one representative row each
    df_pairs = _collapse_overlapping_pairs(df_pairs)
    df_pairs = df_pairs.sort_values('log10_prob_recombination_ecoli', ascending=False)

    if num_sites < np.inf:
        df_pairs = df_pairs.head(num_sites)

    for col in ['start_1', 'end_1', 'start_2', 'end_2', 'location_delta', 'site_length']:
        df_pairs.loc[:, col] = df_pairs[col].astype(int)

    return df_pairs
