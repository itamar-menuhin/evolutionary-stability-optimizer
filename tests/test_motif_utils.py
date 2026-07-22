"""Tests for eso.detection.motif_utils (IUPAC-consensus motif building) and
eso.detection.common_motifs (the bundled Dam/Dcm library) - the easy,
no-MEME-file path for defining or importing motifs.
"""

import numpy as np
import pytest

from eso.detection.common_motifs import COMMON_MOTIFS, load_common_motifs
from eso.detection.methylation import find_motif_sites
from eso.detection.motif_utils import motif_from_consensus, motifs_from_consensus_dict


def test_exact_consensus_scores_a_perfect_match_highest():
    motif = motif_from_consensus("gatc", "GATC")
    seq = "TTTTTT" + "GATC" + "TTTTTT"

    df = find_motif_sites(seq, np.inf, [motif])
    match = df[df.start_index == 6]

    assert len(match) == 1
    assert match.iloc[0].actual_site == "GATC"


def test_ambiguity_code_matches_either_represented_base():
    motif = motif_from_consensus("dcm_like", "CCWGG")  # W = A or T

    for variant in ("CCAGG", "CCTGG"):
        seq = "TTTTTT" + variant + "TTTTTT"
        df = find_motif_sites(seq, np.inf, [motif])
        match = df[df.start_index == 6]
        assert len(match) == 1, f"expected a match for {variant}"
        assert match.iloc[0].actual_site == variant


def test_unrecognized_letter_raises_with_position():
    with pytest.raises(ValueError, match="position 3"):
        motif_from_consensus("bad", "GAXTC")


def test_empty_consensus_raises():
    with pytest.raises(ValueError, match="empty"):
        motif_from_consensus("empty", "")


def test_motifs_from_consensus_dict_builds_all_entries():
    motifs_list = motifs_from_consensus_dict({"a": "GATC", "b": "CCWGG"})
    assert [m.name for m in motifs_list] == ["a", "b"]
    assert len(motifs_list[0]) == 4
    assert len(motifs_list[1]) == 5


def test_common_motifs_covers_methylation_rbs_and_promoter_elements():
    assert set(COMMON_MOTIFS) == {
        "dam", "dcm", "shine_dalgarno", "sigma70_minus35", "sigma70_minus10",
    }
    assert COMMON_MOTIFS["dam"] == "GATC"
    assert COMMON_MOTIFS["dcm"] == "CCWGG"
    assert COMMON_MOTIFS["shine_dalgarno"] == "AGGAGG"
    assert COMMON_MOTIFS["sigma70_minus35"] == "TTGACA"
    assert COMMON_MOTIFS["sigma70_minus10"] == "TATAAT"


def test_load_common_motifs_default_loads_all_five():
    motifs_list = load_common_motifs()
    assert {m.name for m in motifs_list} == set(COMMON_MOTIFS)


def test_shine_dalgarno_finds_a_cryptic_rbs_inside_a_sequence():
    seq = "ATGCATGCAT" + "AGGAGG" + "ATGCATGCAT"
    df = find_motif_sites(seq, np.inf, load_common_motifs(["shine_dalgarno"]))
    match = df[df.start_index == 10]
    assert len(match) == 1
    assert match.iloc[0].actual_site == "AGGAGG"


def test_load_common_motifs_is_case_insensitive_and_selectable():
    motifs_list = load_common_motifs(["DAM"])
    assert [m.name for m in motifs_list] == ["dam"]


def test_load_common_motifs_unknown_name_gives_friendly_error():
    with pytest.raises(ValueError, match="Unknown common motif"):
        load_common_motifs(["not_a_real_motif"])


def test_common_motifs_find_real_sites_in_a_sequence():
    seq = "ATGCATGC" + "GATC" + "ATGCATGC" + "CCAGG" + "ATGCATGC" + "CCTGG" + "ATGC"

    df = find_motif_sites(seq, np.inf, load_common_motifs())

    dam_hits = set(df[df.matching_motif == "dam"].start_index) & {8}
    dcm_hits = set(df[df.matching_motif == "dcm"].start_index)
    assert dam_hits == {8}
    assert dcm_hits == {20, 33}
