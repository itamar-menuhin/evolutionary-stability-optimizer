"""Shared non-max-suppression helper: several detectors (independently, in both
the primary and STABLES-derived implementations) generate multiple candidate
rows - at different phases, frames, or seed positions - that describe the same
real hotspot rather than genuinely distinct ones. This collapses a dataframe
of scored, ranged candidates down to one representative per group of
mutually-overlapping candidates.
"""

import pandas as pd


def ranges_overlap(a, b):
    """True if ranges (start, end) overlap, using the same exclusive-end
    convention as Python slicing (seq[start:end]) - so two ranges that merely
    touch at a shared boundary (e.g. (0, 950) and (950, 1250)) do NOT count
    as overlapping, since they share no actual character position. A
    previous `<=` version treated touching-but-adjacent ranges as
    overlapping, which silently discarded a genuinely distinct, adjacent
    hotspot whenever it happened to sit immediately next to a higher-scoring
    one - found via a chunk-boundary stress test for an unrelated prototype,
    reproduced with two homopolymer runs separated by a third with no gap.
    """
    return a[0] < b[1] and b[0] < a[1]


def collapse_overlapping_intervals(df, score_col, start_col='start', end_col='end'):
    """Walk rows in descending `score_col` order, keeping a row only if its
    [start_col, end_col] range doesn't overlap an already-kept row's range.
    """
    kept_rows = []
    kept_ranges = []

    for _, row in df.sort_values(score_col, ascending=False).iterrows():
        current_range = (row[start_col], row[end_col])
        if any(ranges_overlap(current_range, kept_range) for kept_range in kept_ranges):
            continue
        kept_rows.append(row)
        kept_ranges.append(current_range)

    return pd.DataFrame(kept_rows, columns=df.columns) if kept_rows else df.iloc[0:0]
