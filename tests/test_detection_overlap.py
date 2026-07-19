"""Tests for eso.detection._overlap - the shared collapse helper used by both
recombination and slippage detection.
"""

from eso.detection._overlap import ranges_overlap, collapse_overlapping_intervals
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
