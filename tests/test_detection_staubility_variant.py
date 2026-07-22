"""Direct tests for the STABLES-derived detector variant - separate from
tests/test_detection_dispatch.py, which only exercises it through the mode
router for recombination. Slippage isn't wired into dispatch (both
implementations agree after the fixes below), so it needs its own coverage.
"""

import pytest

from eso.detection.staubility_variant import find_recombination_sites, find_slippage_sites


def test_find_slippage_sites_no_longer_crashes_on_repeat_longer_than_scan_window():
    # Regression test: `df.loc[:, 'end'] = values` silently kept the column's
    # existing float64 dtype (from an earlier `= None` assignment) even after
    # `.astype(int)`, so slicing with a float index raised TypeError - but
    # only for repeats long enough to trigger the back-to-back merge logic.
    seq = "ATGCTAAT" + "GC" * 8 + "TTAGGCATGCCTAGC"
    df = find_slippage_sites(seq)
    assert not df.empty


def test_find_recombination_sites_no_longer_crashes_on_repeat_longer_than_16nt():
    # Same underlying dtype bug, in the sibling recombination function.
    site = "ACGTGGCATTAGCTAGCCTA"  # 20nt, >16
    seq = "ATGCATGCAT" + site + "ATCGGATCCAAGCTTGGATCCAAGCTTGGA" + site + "TTGGCCAATT"
    df = find_recombination_sites(seq)
    assert not df.empty


def test_borderline_homopolymer_run_is_filtered_out():
    # A 4-5nt homopolymer run scores below the -9 risk cutoff (log10 ~= -9.3 to
    # -9.98) and must be filtered, matching eso.detection.slippage's behavior
    # and the filter STABLES itself applied at its call site (not inside this
    # function as originally extracted - see docs/detector-comparisons.md).
    seq = "CCAGCGCGGTCAGTTCCAAAAACACCCTAAGTAACCGAATA"
    df = find_slippage_sites(seq)
    assert df.empty


def test_no_indexerror_when_repeat_sits_at_a_length_boundary():
    # Regression test: is_followed1's `curr_seq_split[ii + 3]` access ran
    # unconditionally whenever the first two equality checks passed, even for
    # length > 1 (where end_of_range only reserves room up to ii+2) - a real
    # bug in the original STABLES source, found via randomized fuzz testing,
    # that crashed with IndexError whenever a length>1 repeat sat at the
    # boundary of its scan range. This exact sequence reproduces it.
    seq = 'CACGCATTTCCCCCCTACATCACCAGAGAG'
    df = find_slippage_sites(seq)
    assert not df.empty


def test_homopolymer_run_of_11_found_via_length_2_run_of_12_via_length_1():
    # Same crossover as eso.detection.slippage (see its module docstring):
    # n=11 is still detected, but only via an independently-run length-2
    # reading (which outscores length-1's); n=12 is the first case where the
    # length-1 reading itself wins - confirms the detection seed here (12)
    # matches the primary implementation's boundary exactly.
    seq_11 = "ATGCTAGCCATTAGGC" + "A" * 11 + "TTGCCTAGCATGC"
    df_11 = find_slippage_sites(seq_11)
    hit_11 = df_11[(df_11.start <= 27) & (df_11.end >= 16)]
    assert not hit_11.empty
    assert hit_11.iloc[0].length_base_unit == 2

    seq_12 = "ATGCTAGCCATTAGGC" + "A" * 12 + "TTGCCTAGCATGC"
    df_12 = find_slippage_sites(seq_12)
    hit_12 = df_12[(df_12.start <= 28) & (df_12.end >= 16)]
    assert not hit_12.empty
    assert hit_12.iloc[0].length_base_unit == 1
    assert hit_12.iloc[0].num_base_units == 12


def test_no_redundant_rows_for_same_start_competing_lengths():
    # A run of 6 T's is valid as both a length-1 and length-2 site at the same
    # start; only the highest-scoring representation should survive.
    seq = "ATGCTAGCCATTAGGC" + "TTTTTT" + "ATGCCTAGCATGC"
    df = find_slippage_sites(seq)
    assert df.shape[0] == 1


def test_recombination_score_matches_efm_calculator_reference():
    # Same fix and same reference cross-check as
    # eso.detection.recombination.test_calc_recombination_score_matches_efm_calculator_reference -
    # this module duplicates the formula inline rather than calling the
    # shared function, so it needs its own pin against the reference.
    import math

    site = "ACGTGGCATTAGCTAGCCTA"  # 20nt
    seq = "ATGCATGCAT" + site + "ATCGGATCCAAGCTTGGATCCAAGCTTGGA" + site + "TTGGCCAATT"
    df = find_recombination_sites(seq)
    assert not df.empty

    row = df.iloc[0]
    location_delta = row.start_2 - row.end_1
    site_length = row.end_1 - row.start_1
    reference = math.log10(
        ((5.8 + location_delta) ** (-29.0 / site_length))
        * (site_length / (1 + 1465.6 * site_length))
    )
    assert row.log10_prob_recombination_ecoli == pytest.approx(reference, abs=1e-6)
