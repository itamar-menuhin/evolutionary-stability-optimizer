"""Convert detected hotspot dataframes into DNAChisel AvoidPattern constraints,
respecting user-specified exclusion (locked) regions.
"""

from os import path

import dnachisel
import pandas as pd

from eso.detection.recombination import _generate_neighbors


def has_overlap_exclusion(start, end, exclusions):
    """True if the (start, end) region overlaps any exclusion region.

    (start, end) is exclusive-end throughout this codebase (matches Python
    slicing, eso.sequence_utils.parse_region's output, and
    eso.detection.recombination's output after _elongate_sites) - a region
    ending exactly where another begins shares no actual nucleotide with it,
    so `<=`/`>=` (not `<`/`>`) is required here to avoid treating merely
    touching regions as overlapping (the same bug class fixed in
    eso.detection._overlap.ranges_overlap, found independently in this
    separate module).
    """
    for ex in exclusions:
        if end <= ex[0] or start >= ex[1]:
            continue
        return True
    return False


def _indel_recombinations(row, exclusions):
    """For an indel-type recombination pair: enforce a change in the smaller
    region if it doesn't overlap an exclusion; otherwise enforce a change in
    the larger region's inserted nucleotide (the only edit that doesn't
    increase the Levenshtein distance further).
    """
    start_small, end_small, sequence_small = row.start_1, row.end_1, row.sequence_1
    start_large, end_large, sequence_large = row.start_2, row.end_2, row.sequence_2
    if row.len_1 > row.len_2:
        start_small, end_small, sequence_small = row.start_2, row.end_2, row.sequence_2
        start_large, end_large, sequence_large = row.start_1, row.end_1, row.sequence_1

    if not has_overlap_exclusion(start_small, end_small, exclusions):
        return [(start_small, end_small, sequence_small)]

    if has_overlap_exclusion(start_large, end_large, exclusions):
        return [(start_small, end_small, sequence_small)]

    prefix = path.commonprefix([sequence_small, sequence_large])
    suffix = sequence_small[len(prefix):] if len(prefix) < len(sequence_small) else ''
    return [(start_large, end_large, prefix + nt + suffix) for nt in ['A', 'C', 'G', 'T']]


def _substitution_recombinations(row, exclusions):
    """For a substitution-type recombination pair: enforce a change in
    whichever region doesn't overlap an exclusion, at all its single-substitution neighbors.
    """
    start_1, end_1, sequence_1 = row.start_1, row.end_1, row.sequence_1
    start_2, end_2 = row.start_2, row.end_2

    if has_overlap_exclusion(start_2, end_2, exclusions) and not has_overlap_exclusion(start_1, end_1, exclusions):
        start_2, end_2 = start_1, end_1

    substitution_neighbours, _, _ = _generate_neighbors(sequence_1)
    return [(start_2, end_2, sub) for sub in substitution_neighbours]


def recombination_to_multiple_avoidance_sites(df, exclusion_regions):
    """Translate detected recombination pairs into individual sites to mutate,
    such that mutating them breaks the Levenshtein-distance-1 relationship
    between the pair (see _indel_recombinations / _substitution_recombinations).
    """
    df_copy = df[['start_1', 'end_1', 'sequence_1', 'start_2', 'end_2', 'sequence_2']].copy()
    df_copy.loc[:, 'len_1'] = df_copy.sequence_1.apply(len)
    df_copy.loc[:, 'len_2'] = df_copy.sequence_2.apply(len)

    recombination_sites = []
    for ii in range(df_copy.shape[0]):
        row = df_copy.iloc[ii]
        if row.len_1 != row.len_2:
            recombination_sites.extend(_indel_recombinations(row, exclusion_regions))
        else:
            recombination_sites.extend(_substitution_recombinations(row, exclusion_regions))

    recombination_sites = sorted(set(recombination_sites))
    return pd.DataFrame.from_records(data=recombination_sites, columns=['start', 'end', 'sequence'])


def exclusion_site_correcter(df, exclusion_regions):
    """Trim detected sites so they don't overlap any exclusion (locked) region.

    (start, end) is exclusive-end (see has_overlap_exclusion) - a region's
    exclusive end is itself the first *excluded* position, so trimming a site
    to stop before an exclusion should set end = region[0] (not
    region[0] - 1), and trimming to start after one should set
    start = region[1] (not region[1] + 1). The old `-1`/`+1` didn't corrupt
    the `sequence` field (it always matched the true substring at the
    resulting start/end - just a shorter one than necessary), but it did
    needlessly discard one usable, still-modifiable nucleotide per exclusion
    boundary.
    """
    if exclusion_regions == 'error':
        return df

    for region in exclusion_regions:
        df_before = df[df.start < region[0]]
        df_before.loc[:, 'end'] = df_before.end.apply(lambda x: int(min(x, region[0])))
        df_before.loc[:, 'sequence'] = df_before.apply(lambda row: row.sequence[:(row.end - row.start)], axis=1)

        df_after = df[df.end > region[1]]
        df_after.loc[:, 'start'] = df_after.start.apply(lambda x: int(max(x, region[1])))
        df_after.loc[:, 'sequence'] = df_after.apply(lambda row: row.sequence[-(row.end - row.start):], axis=1)

        df = pd.concat([df_before, df_after], ignore_index=True)

        # keep at least one codon, so there's still something to modify
        df = df[df.start < df.end - 2]

    return df.sort_values('start').drop_duplicates().reset_index(drop=True)


def convert_df_to_constraints(df):
    """Convert a {sequence, start, end} dataframe of patterns-to-avoid into
    DNAChisel AvoidPattern constraints.
    """
    if df.shape[0] == 0:
        return []

    df = df[['sequence', 'start', 'end']].drop_duplicates()
    return [
        dnachisel.AvoidPattern(df.loc[idx, 'sequence'], location=(int(df.loc[idx, 'start']), int(df.loc[idx, 'end'])))
        for idx in df.index
    ]
