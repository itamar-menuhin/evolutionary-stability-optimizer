"""Tests for eso.detection._overlap - the shared collapse helpers used by both
recombination and slippage detection.
"""

from eso.detection._overlap import (
    ranges_overlap,
    collapse_overlapping_intervals,
    collapse_overlapping_intervals_no_coverage_loss,
)
import pandas as pd


def test_touching_ranges_do_not_overlap():
    # (start, end) uses Python-slicing exclusive-end semantics (seq[start:end]),
    # so a range ending at 950 and one starting at 950 share no character
    # position and must not count as overlapping.
    assert not ranges_overlap((0, 950), (950, 1250))


def test_genuinely_overlapping_ranges_do_overlap():
    assert ranges_overlap((0, 950), (940, 1250))
    assert ranges_overlap((100, 200), (150, 160))  # fully contained


def test_disjoint_ranges_do_not_overlap():
    assert not ranges_overlap((0, 100), (200, 300))


def test_collapse_keeps_adjacent_non_overlapping_regions():
    # Regression test: a `<=`-based overlap check previously treated two
    # ranges that merely touch (e.g. (0, 950) and (950, 1250)) as
    # overlapping, silently discarding the lower-scoring one even though
    # they share no actual character position. Reproduced with
    # eso.detection.slippage.find_slippage_sites on "A"*950 + "T"*300 + "A"*950 -
    # the middle T-run used to vanish entirely.
    df = pd.DataFrame([
        {'start': 0, 'end': 950, 'score': 679.0},
        {'start': 950, 'end': 1250, 'score': 205.8},
        {'start': 1250, 'end': 2200, 'score': 679.0},
    ])
    collapsed = collapse_overlapping_intervals(df, score_col='score')
    assert collapsed.shape[0] == 3


def test_no_coverage_loss_still_collapses_fully_contained_candidates():
    # A lower-scoring candidate whose entire range is already inside a
    # higher-scoring kept one is genuinely redundant - the common case
    # (several detections converging on the same physical site) - and must
    # still collapse away exactly like plain NMS does.
    df = pd.DataFrame([
        {'start': 0, 'end': 20, 'score': 10.0},
        {'start': 2, 'end': 18, 'score': 5.0},  # fully inside [0, 20)
    ])
    collapsed = collapse_overlapping_intervals_no_coverage_loss(df, score_col='score')
    assert collapsed.shape[0] == 1
    assert (collapsed.iloc[0].start, collapsed.iloc[0].end) == (0, 20)


def test_no_coverage_loss_keeps_a_partially_overlapping_distinct_candidate():
    # Regression test for the real coverage-gap bug: two genuinely distinct
    # candidates whose ranges overlap but neither contains the other. Plain
    # NMS (collapse_overlapping_intervals) drops the lower-scoring one
    # entirely, silently leaving positions [10, 20) - the part only the
    # dropped candidate covered - with no representative at all. The
    # no-coverage-loss variant must keep both, since the higher-scoring one
    # doesn't cover the dropped one's full extent.
    df = pd.DataFrame([
        {'start': 0, 'end': 10, 'score': 10.0},   # scores higher, kept either way
        {'start': 8, 'end': 20, 'score': 5.0},    # overlaps [0, 10) only at [8, 10)
    ])

    plain = collapse_overlapping_intervals(df, score_col='score')
    assert plain.shape[0] == 1  # the real bug: [10, 20) loses all coverage

    no_loss = collapse_overlapping_intervals_no_coverage_loss(df, score_col='score')
    assert no_loss.shape[0] == 2
    assert set(zip(no_loss.start, no_loss.end)) == {(0, 10), (8, 20)}
