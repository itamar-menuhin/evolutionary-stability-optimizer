"""Direct tests for the STABLES-derived detector variant - separate from
tests/test_detection_dispatch.py, which only exercises it through the mode
router for recombination. Slippage isn't wired into dispatch (both
implementations agree after the fixes below), so it needs its own coverage.
"""

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


def test_six_nucleotide_homopolymer_run_is_still_detected():
    # n=6 is the exact minimum that survives the -9 filter; confirms the
    # detection seed raised from 4 to 6 (see eso.detection.slippage's module
    # docstring) didn't accidentally exclude the boundary case here too.
    seq = "ATGCTAGCCATTAGGC" + "AAAAAA" + "ATGCCTAGCATGC"
    df = find_slippage_sites(seq)
    assert (df.sequence.str.contains("AAAAAA")).any()


def test_no_redundant_rows_for_same_start_competing_lengths():
    # A run of 6 T's is valid as both a length-1 and length-2 site at the same
    # start; only the highest-scoring representation should survive.
    seq = "ATGCTAGCCATTAGGC" + "TTTTTT" + "ATGCCTAGCATGC"
    df = find_slippage_sites(seq)
    assert df.shape[0] == 1
