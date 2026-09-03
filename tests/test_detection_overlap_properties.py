"""Property-based (fuzz) tests for eso.detection._overlap.

The example-based tests in tests/test_detection_overlap.py pin specific,
hand-picked scenarios (including the exact coverage-gap bug this module was
fixed for - see docs/detector-comparisons.md). This file instead generates
many random interval sets and checks invariants that must hold for ANY input,
not just the cases someone thought to write by hand - the overlap-collapse
logic in this codebase has a real history of edge cases (touching-but-not-
overlapping ranges, partially-overlapping-but-not-nested ranges) that hand-
picked examples missed the first time around, so this is exactly the kind of
module worth fuzzing.
"""

import pandas as pd
from hypothesis import given, strategies as st

from eso.detection._overlap import (
    ranges_overlap,
    collapse_overlapping_intervals,
    collapse_overlapping_intervals_no_coverage_loss,
)

# A single (start, end, score) row: end > start always (a zero/negative-length
# range isn't a real interval any detector in this codebase would produce).
_interval = st.tuples(
    st.integers(min_value=0, max_value=200),
    st.integers(min_value=1, max_value=50),
    st.floats(min_value=-100, max_value=100, allow_nan=False, allow_infinity=False),
).map(lambda t: (t[0], t[0] + t[1], t[2]))  # (start, start+positive_length, score)

_intervals_df = st.lists(_interval, min_size=0, max_size=25).map(
    lambda rows: pd.DataFrame(rows, columns=['start', 'end', 'score'])
)


def _fully_contains(outer_start, outer_end, inner_start, inner_end):
    return outer_start <= inner_start and inner_end <= outer_end


@given(a=_interval, b=_interval)
def test_ranges_overlap_is_symmetric(a, b):
    range_a, range_b = (a[0], a[1]), (b[0], b[1])
    assert ranges_overlap(range_a, range_b) == ranges_overlap(range_b, range_a)


@given(df=_intervals_df)
def test_collapsed_output_never_mutually_overlaps(df):
    # the defining guarantee of non-max suppression: no two surviving rows
    # can overlap each other, or collapsing didn't actually happen.
    collapsed = collapse_overlapping_intervals(df, score_col='score')
    ranges = list(zip(collapsed.start, collapsed.end))
    for i in range(len(ranges)):
        for j in range(i + 1, len(ranges)):
            assert not ranges_overlap(ranges[i], ranges[j])


@given(df=_intervals_df)
def test_collapsed_output_is_a_subset_of_the_input_rows(df):
    collapsed = collapse_overlapping_intervals(df, score_col='score')
    input_rows = set(zip(df.start, df.end, df.score))
    output_rows = set(zip(collapsed.start, collapsed.end, collapsed.score))
    assert output_rows <= input_rows


@given(df=_intervals_df)
def test_no_coverage_loss_output_is_a_subset_of_the_input_rows(df):
    reduced = collapse_overlapping_intervals_no_coverage_loss(df, score_col='score')
    input_rows = set(zip(df.start, df.end, df.score))
    output_rows = set(zip(reduced.start, reduced.end, reduced.score))
    assert output_rows <= input_rows


@given(df=_intervals_df)
def test_no_coverage_loss_never_drops_a_row_without_something_covering_it(df):
    # THE core property this whole module was fixed to guarantee: every
    # input row either survives, or its entire range is covered by some row
    # that DID survive. A row silently vanishing with nothing covering its
    # range is exactly the coverage-gap bug documented in
    # docs/detector-comparisons.md - this test would have caught it directly,
    # for any input, not just the one hand-constructed scenario that found it.
    reduced = collapse_overlapping_intervals_no_coverage_loss(df, score_col='score')
    kept_ranges = list(zip(reduced.start, reduced.end))

    for start, end, _score in zip(df.start, df.end, df.score):
        covered = any(_fully_contains(k_start, k_end, start, end) for k_start, k_end in kept_ranges)
        assert covered, f"input row ({start}, {end}) has no covering row in the output"
