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


def range_contains(outer, inner):
    """True if `inner` (start, end) is fully contained in `outer` (start, end),
    exclusive-end throughout (matches ranges_overlap's convention).
    """
    return outer[0] <= inner[0] and inner[1] <= outer[1]


def _collapse_by_predicate(df, score_col, start_col, end_col, should_drop):
    """Shared non-max-suppression walk: process rows in descending `score_col`
    order, dropping a row exactly when `should_drop(current_range, kept_range)`
    is True for some already-kept range, otherwise keeping it. The two public
    collapse functions below differ only in which predicate they pass in -
    everything else (the sort, the kept_rows/kept_ranges bookkeeping, the
    empty-dataframe fallback) is identical between them, so it lives here once.
    """
    kept_rows = []
    kept_ranges = []

    for _, row in df.sort_values(score_col, ascending=False).iterrows():
        current_range = (row[start_col], row[end_col])
        if any(should_drop(current_range, kept_range) for kept_range in kept_ranges):
            continue
        kept_rows.append(row)
        kept_ranges.append(current_range)

    return pd.DataFrame(kept_rows, columns=df.columns) if kept_rows else df.iloc[0:0]


def collapse_overlapping_intervals(df, score_col, start_col='start', end_col='end'):
    """Walk rows in descending `score_col` order, keeping a row only if its
    [start_col, end_col] range doesn't overlap an already-kept row's range.

    This is non-max suppression: it only requires ranges to *overlap*, not
    that one *contain* the other, to drop the lower-scoring one. For a
    "distinct sites" report or count that's the right call - but it means a
    row can be dropped even though part of its range isn't actually covered
    by anything kept, silently losing coverage of that part entirely. Do NOT
    use this to build correction constraints for that reason - use
    collapse_overlapping_intervals_no_coverage_loss instead. See
    docs/detector-comparisons.md for the concrete failure this caused.
    """
    return _collapse_by_predicate(df, score_col, start_col, end_col, should_drop=ranges_overlap)


def collapse_overlapping_intervals_no_coverage_loss(df, score_col, start_col='start', end_col='end'):
    """Like collapse_overlapping_intervals, but only drops a row when its
    range is FULLY CONTAINED in an already-kept row's range - not merely
    overlapping it.

    Walking in descending `score_col` order, a redundant candidate (one whose
    entire range is already covered by a higher-scoring kept row - the common
    case: several detections of the same physical site converging on
    essentially the same extent) is still dropped, exactly as before. But a
    candidate that only *partially* overlaps a kept row - some real,
    independently-detected hotspot whose span isn't a subset of the kept
    row's - is kept too, so nothing downstream ever has to build a correction
    constraint for a region no surviving row actually covers.

    Rows returned by this function can legitimately still overlap each other
    (that's the point) - use it to feed a correction/constraint-building step,
    never for a "how many distinct sites" count or report (use
    collapse_overlapping_intervals for that instead).
    """
    def _kept_fully_covers_current(current_range, kept_range):
        return range_contains(outer=kept_range, inner=current_range)

    return _collapse_by_predicate(df, score_col, start_col, end_col, should_drop=_kept_fully_covers_current)
