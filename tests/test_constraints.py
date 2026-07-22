"""Tests for eso.constraints - previously untested, despite sitting directly
between hotspot detection and DNAChisel optimization (a bug here silently
changes what actually gets optimized, with no visible symptom).

Two real bugs found and fixed here, both at range boundaries under this
codebase's exclusive-end convention (matching Python slicing and
eso.sequence_utils.parse_region's output, and the same convention
eso.detection.recombination's output uses after _elongate_sites):

1. has_overlap_exclusion treated a site that merely *touches* an exclusion
   region (shares no actual nucleotide with it) as overlapping it - the same
   "touching ranges" bug class already found and fixed in
   eso.detection._overlap.ranges_overlap, but in this separate, never-touched
   module.
2. exclusion_site_correcter trimmed one nucleotide more than necessary on
   both sides of an exclusion region (using region[0]-1 / region[1]+1
   instead of region[0] / region[1]) - not a data-corruption bug (the
   trimmed `sequence` field always matched the true substring at its new
   `start`/`end`), but a real, needless loss of one nucleotide of usable,
   still-modifiable sequence per exclusion boundary.
"""

import pandas as pd
import pytest

from eso.constraints import (
    convert_df_to_constraints,
    exclusion_site_correcter,
    has_overlap_exclusion,
    recombination_to_multiple_avoidance_sites,
)


# --- has_overlap_exclusion ---------------------------------------------

def test_touching_ranges_do_not_overlap():
    # site [5,10) and exclusion [10,15) share no actual nucleotide - they merely touch.
    assert has_overlap_exclusion(5, 10, [(10, 15)]) is False


def test_touching_ranges_do_not_overlap_other_side():
    assert has_overlap_exclusion(15, 20, [(10, 15)]) is False


def test_genuinely_overlapping_ranges_do_overlap():
    assert has_overlap_exclusion(5, 11, [(10, 15)]) is True


def test_disjoint_ranges_do_not_overlap():
    assert has_overlap_exclusion(0, 3, [(10, 15)]) is False


def test_checks_against_every_exclusion_in_the_list():
    assert has_overlap_exclusion(20, 25, [(0, 5), (20, 30)]) is True
    assert has_overlap_exclusion(20, 25, [(0, 5), (40, 50)]) is False


# --- exclusion_site_correcter -------------------------------------------

FULL_SEQ = "ACGTACGTACGTACGTACGTACGT"  # 24nt


def test_trims_to_the_maximal_safe_region_before_an_exclusion():
    # exclusion [8,15) excludes positions 8-14; position 7 is still free and
    # should be kept, not discarded.
    site = pd.DataFrame([{"start": 3, "end": 10, "sequence": FULL_SEQ[3:10]}])

    result = exclusion_site_correcter(site, [(8, 15)])

    assert len(result) == 1
    row = result.iloc[0]
    assert (row.start, row.end) == (3, 8)
    assert row.sequence == FULL_SEQ[row.start:row.end]


def test_trims_to_the_maximal_safe_region_after_an_exclusion():
    # exclusion [1,5) excludes positions 1-4; position 5 is still free and
    # should be kept, not discarded.
    site = pd.DataFrame([{"start": 3, "end": 10, "sequence": FULL_SEQ[3:10]}])

    result = exclusion_site_correcter(site, [(1, 5)])

    assert len(result) == 1
    row = result.iloc[0]
    assert (row.start, row.end) == (5, 10)
    assert row.sequence == FULL_SEQ[row.start:row.end]


def test_trimmed_sequence_always_matches_the_true_substring():
    # general correctness check across several boundary placements, not just
    # the two cases above.
    site = pd.DataFrame([{"start": 2, "end": 16, "sequence": FULL_SEQ[2:16]}])

    for exclusion in [(0, 3), (5, 9), (13, 16), (15, 20)]:
        result = exclusion_site_correcter(site, [exclusion])
        for _, row in result.iterrows():
            assert row.sequence == FULL_SEQ[row.start:row.end]


def test_site_entirely_inside_an_exclusion_is_dropped():
    site = pd.DataFrame([{"start": 5, "end": 8, "sequence": FULL_SEQ[5:8]}])
    result = exclusion_site_correcter(site, [(0, 20)])
    assert result.empty


def test_site_untouched_by_any_exclusion_is_unchanged():
    site = pd.DataFrame([{"start": 5, "end": 12, "sequence": FULL_SEQ[5:12]}])
    result = exclusion_site_correcter(site, [(15, 20)])
    assert len(result) == 1
    row = result.iloc[0]
    assert (row.start, row.end, row.sequence) == (5, 12, FULL_SEQ[5:12])


def test_error_sentinel_passes_through_unchanged():
    site = pd.DataFrame([{"start": 5, "end": 12, "sequence": FULL_SEQ[5:12]}])
    result = exclusion_site_correcter(site, "error")
    assert result is site


# --- recombination_to_multiple_avoidance_sites ---------------------------

def test_substitution_pair_targets_region_2_with_all_single_substitution_neighbors():
    # equal-length pair -> substitution-type: forces region_2 away from every
    # sequence within Levenshtein distance 1 of sequence_1 (including an
    # exact match), so the resulting sequence must end up at distance >= 2.
    df = pd.DataFrame([{
        "start_1": 0, "end_1": 6, "sequence_1": "ACGTAC",
        "start_2": 20, "end_2": 26, "sequence_2": "ACGTAG",
    }])

    result = recombination_to_multiple_avoidance_sites(df, ())

    assert set(result.start) == {20}
    assert set(result.end) == {26}
    # exact match + 3 substitutions/position * 6 positions = 19 neighbors
    assert len(result) == 19
    assert "ACGTAC" in set(result.sequence)  # the exact original is included


def test_substitution_pair_targets_region_1_if_region_2_is_excluded():
    df = pd.DataFrame([{
        "start_1": 0, "end_1": 6, "sequence_1": "ACGTAC",
        "start_2": 20, "end_2": 26, "sequence_2": "ACGTAG",
    }])

    result = recombination_to_multiple_avoidance_sites(df, [(20, 26)])

    assert set(result.start) == {0}
    assert set(result.end) == {6}


def test_indel_pair_targets_the_smaller_region_when_not_excluded():
    # unequal lengths -> indel-type: mutate the smaller region directly.
    df = pd.DataFrame([{
        "start_1": 0, "end_1": 6, "sequence_1": "ACGTAC",       # smaller (6nt)
        "start_2": 20, "end_2": 27, "sequence_2": "ACGTACG",    # larger (7nt)
    }])

    result = recombination_to_multiple_avoidance_sites(df, ())

    assert len(result) == 1
    assert (result.iloc[0].start, result.iloc[0].end) == (0, 6)


def test_indel_pair_targets_the_larger_regions_insertion_when_smaller_is_excluded():
    df = pd.DataFrame([{
        "start_1": 0, "end_1": 6, "sequence_1": "ACGTAC",
        "start_2": 20, "end_2": 27, "sequence_2": "ACGTACG",
    }])

    result = recombination_to_multiple_avoidance_sites(df, [(0, 6)])

    # targets the larger region (7nt window), one AvoidPattern per possible
    # inserted nucleotide (4 candidates)
    assert set(result.start) == {20}
    assert set(result.end) == {27}
    assert len(result) == 4


# --- convert_df_to_constraints -------------------------------------------

def test_empty_dataframe_returns_no_constraints():
    df = pd.DataFrame(columns=["sequence", "start", "end"])
    assert convert_df_to_constraints(df) == []


def test_builds_one_avoidpattern_per_row():
    df = pd.DataFrame([
        {"sequence": "ACGT", "start": 0, "end": 4},
        {"sequence": "TTTT", "start": 10, "end": 14},
    ])

    constraints = convert_df_to_constraints(df)

    assert len(constraints) == 2
    locations = sorted((c.location.start, c.location.end) for c in constraints)
    assert locations == [(0, 4), (10, 14)]


def test_duplicate_rows_produce_one_constraint():
    df = pd.DataFrame([
        {"sequence": "ACGT", "start": 0, "end": 4},
        {"sequence": "ACGT", "start": 0, "end": 4},
    ])
    assert len(convert_df_to_constraints(df)) == 1
