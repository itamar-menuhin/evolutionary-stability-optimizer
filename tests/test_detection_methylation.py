"""Tests for eso.detection.methylation.find_motif_sites - previously untested.

find_motif_sites was rewritten from a per-position/per-motif Python loop that
built a 2*num_motifs*len(seq)-row intermediate dataframe (profiling showed this,
plus a df.apply(axis=1) to extract each site's sequence, dominated runtime -
~14x speedup measured after vectorizing with numpy, see docs/detector-comparisons.md)
into a vectorized numpy score-matrix reduction. These tests pin its behavior
directly, and a couple explicitly document the one known, pre-existing
ambiguity: when two *different* motifs score exactly equally at the same
position, which one is reported was never a deliberate, documented choice
(the original code's answer was an accident of pandas' unstable sort) - this
rewrite picks the lowest motif_number deterministically instead.
"""

import io

import numpy as np
import pytest
from Bio import motifs

from eso.detection.methylation import find_motif_sites


def _parse_meme(matrix_rows, name="test_motif", background=(0.25, 0.25, 0.25, 0.25)):
    bg_line = " ".join(f"{letter} {freq}" for letter, freq in zip("ACGT", background))
    matrix_lines = "\n".join(" ".join(f"{p:.3f}" for p in row) for row in matrix_rows)
    meme_text = f"""MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
{bg_line}

MOTIF {name}
letter-probability matrix: alength= 4 w= {len(matrix_rows)} nsites= 20 E= 0
{matrix_lines}
"""
    return list(motifs.parse(io.StringIO(meme_text), "minimal"))[0]


# strongly prefers A, C, G, T at positions 0-3 respectively
_ACGT_MOTIF_ROWS = [
    [0.85, 0.05, 0.05, 0.05],
    [0.05, 0.85, 0.05, 0.05],
    [0.05, 0.05, 0.85, 0.05],
    [0.05, 0.05, 0.05, 0.85],
]


def test_finds_exact_match():
    motif = _parse_meme(_ACGT_MOTIF_ROWS)
    seq = "GGGGGG" + "ACGT" + "GGGGGG"

    df = find_motif_sites(seq, np.inf, [motif])

    assert len(df) == 1
    row = df.iloc[0]
    assert row.start_index == 6
    assert row.end_index == 9
    assert row.actual_site == "ACGT"
    assert row.matching_motif == "test_motif"


def test_no_match_in_unrelated_sequence():
    motif = _parse_meme(_ACGT_MOTIF_ROWS)
    seq = "TTTTTTTTTTTTTTTTTTTT"  # no A/C/G/T-in-order 4-mer at all

    df = find_motif_sites(seq, np.inf, [motif])

    assert df.empty


def test_finds_match_on_reverse_complement_strand():
    # prefers A,A,C,C at positions 0-3; reverse complement of "AACC" is "GGTT",
    # so a motif that only matches "AACC" on the forward strand should still
    # flag a "GGTT" window via its reverse-complement scan.
    rows = [
        [0.85, 0.05, 0.05, 0.05],
        [0.85, 0.05, 0.05, 0.05],
        [0.05, 0.85, 0.05, 0.05],
        [0.05, 0.85, 0.05, 0.05],
    ]
    motif = _parse_meme(rows, name="aacc_motif")
    seq = "TCTCTC" + "GGTT" + "TCTCTC"  # neutral, non-repetitive flanks

    df = find_motif_sites(seq, np.inf, [motif])

    match = df[df.start_index == 6]
    assert len(match) == 1
    row = match.iloc[0]
    assert row.actual_site == "GGTT"
    assert row.actual_site_reverse_conjugate == "AACC"


def test_keeps_only_best_scoring_motif_per_position():
    strong = _parse_meme(_ACGT_MOTIF_ROWS, name="strong")
    weak_rows = [[0.4, 0.2, 0.2, 0.2]] * 4
    weak = _parse_meme(weak_rows, name="weak")
    seq = "GGGGGG" + "ACGT" + "GGGGGG"

    df = find_motif_sites(seq, np.inf, [weak, strong])

    matches_at_6 = df[df.start_index == 6]
    assert len(matches_at_6) == 1
    assert matches_at_6.iloc[0].matching_motif == "strong"


def test_num_sites_truncates_to_highest_scoring():
    motif = _parse_meme(_ACGT_MOTIF_ROWS)
    seq = "ACGT" * 20  # many tandem exact matches, identical top score

    df_all = find_motif_sites(seq, np.inf, [motif])
    df_top3 = find_motif_sites(seq, 3, [motif])

    assert len(df_top3) == 3
    assert set(df_top3.start_index).issubset(set(df_all.start_index))


def test_output_sorted_by_score_descending():
    strong = _parse_meme(_ACGT_MOTIF_ROWS, name="strong")
    weak_rows = [[0.4, 0.2, 0.2, 0.2]] * 4
    weak = _parse_meme(weak_rows, name="weak")
    seq = "ACGT" + "GGGG" + "ACGT"  # strong match at 0 and 8, weak-only region in between

    df = find_motif_sites(seq, np.inf, [weak, strong])

    assert list(df.PSSM_score) == sorted(df.PSSM_score, reverse=True)


def test_empty_motif_list_returns_empty_dataframe():
    df = find_motif_sites("ACGTACGT", np.inf, [])
    assert df.empty
    assert list(df.columns) == [
        'start_index', 'end_index', 'matching_motif', 'PSSM_score',
        'actual_site', 'actual_site_reverse_conjugate',
    ]


def test_motif_longer_than_sequence_is_skipped_not_crashed():
    long_motif = _parse_meme([[0.7, 0.1, 0.1, 0.1]] * 20, name="long")
    df = find_motif_sites("ACGT", np.inf, [long_motif])
    assert df.empty


def test_tie_between_distinct_motifs_resolves_to_lower_motif_number():
    # Two motifs scoring identically at the same position is a genuine,
    # previously-undocumented ambiguity (see module docstring) - this pins the
    # deterministic tie-break this rewrite chose (lowest motif_number first),
    # so the behavior is at least reproducible and explicit going forward.
    rows = [[0.7, 0.1, 0.1, 0.1]] * 4
    motif_a = _parse_meme(rows, name="a")
    motif_b = _parse_meme(rows, name="b")
    seq = "GGGGGG" + "AAAA" + "GGGGGG"

    df = find_motif_sites(seq, np.inf, [motif_a, motif_b])
    at_pos = df[df.start_index == 6]
    assert len(at_pos) == 1
    assert at_pos.iloc[0].matching_motif == "a"

    df_reordered = find_motif_sites(seq, np.inf, [motif_b, motif_a])
    at_pos_reordered = df_reordered[df_reordered.start_index == 6]
    assert at_pos_reordered.iloc[0].matching_motif == "b"
