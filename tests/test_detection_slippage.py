from eso.detection.slippage import find_slippage_sites, modify_df_slippage


def test_finds_dinucleotide_repeat():
    # "GC" repeated 4x is a length-2 base-unit slippage site. Flanking bases are
    # chosen so they don't accidentally extend the "GC" run in either direction.
    seq = "ATGCTAAT" + "GCGCGCGC" + "TTAGGCATGCCTAGC"
    df = find_slippage_sites(seq)

    # collapsed to exactly one row: "GCGCGCGC" starting at position 8 and its
    # 1-shifted reading "CGCGCG" starting at position 9 are the same physical
    # repeat and must not be reported as two separate sites.
    assert df.shape[0] == 1
    row = df.iloc[0]
    assert row.length_base_unit == 2
    assert row.num_base_units == 4


def test_finds_single_nucleotide_run():
    # A run of 6 T's is simultaneously a valid length-1 site (TxTxTxTxTxT) and a
    # length-2 site (TT x3); find_slippage_sites keeps only the single
    # highest-scoring representation per start position, which here is the
    # length-2 one (log10 prob -4.56 vs. -8.526 for the length-1 reading).
    seq = "ATGCTAGCCATTAGGC" + "TTTTTT" + "ATGCCTAGCATGC"
    df = find_slippage_sites(seq)

    matches = df[df.sequence.str.contains("TTTTTT")]
    assert not matches.empty
    assert matches.iloc[0].length_base_unit == 2
    assert matches.iloc[0].num_base_units == 3


def test_no_repeats_returns_empty():
    # Hand-checked to contain no repeated base unit of length 1-15.
    seq = "ACGTAGCTTGACCTGAAGCTAGCA"
    df = find_slippage_sites(seq)
    assert df.empty


def test_short_homopolymer_runs_below_filter_threshold_are_empty():
    # A run of 4 or 5 identical nucleotides always fails the -9 log10-prob
    # filter outright (n=4 -> -9.984, n=5 -> -9.255) - no length-2 fallback
    # rescues these (a length-2 "XX" site needs 6nt minimum), so they
    # correctly produce no detection at all, at any length_base_unit.
    seq_4 = "ATGCTAGCCATTAGGC" + "AAAA" + "TTGCCTAGCATGC"
    seq_5 = "ATGCTAGCCATTAGGC" + "AAAAA" + "TTGCCTAGCATGC"
    assert find_slippage_sites(seq_4).empty
    assert find_slippage_sites(seq_5).empty


def test_homopolymer_runs_6_to_11_are_found_via_length_2_not_length_1():
    # n=6..11 all pass the -9 filter at length_base_unit=1, but an
    # independently-run length-2 detection of the same region always scores
    # higher (e.g. n=8: length-1 -7.068 vs length-2 -4.497) and wins the
    # overlap collapse - so the region is still found, just reported as a
    # length-2 site, never as length-1. This is why the length-1 detection
    # seed is 12, not 6: below 12 the length-1 reading never survives anyway.
    for n in (6, 8, 11):
        seq = "ATGCTAGCCATTAGGC" + "A" * n + "TTGCCTAGCATGC"
        df = find_slippage_sites(seq)
        hit = df[(df.start <= 16 + n) & (df.end >= 16)]
        assert not hit.empty, f"n={n} should still be detected (via length-2)"
        assert hit.iloc[0].length_base_unit == 2, f"n={n} should be reported as length-2, not length-1"


def test_homopolymer_run_of_12_is_the_first_to_win_as_length_1():
    # n=12 is the exact crossover: length-1 score -4.152 first exceeds
    # length-2's -4.371, so from here the length-1 reading survives collapse
    # instead - confirming the detection seed (12) doesn't cut off the first
    # case where length-1 actually matters.
    seq = "ATGCTAGCCATTAGGC" + "A" * 12 + "TTGCCTAGCATGC"
    df = find_slippage_sites(seq)
    hit = df[(df.start <= 28) & (df.end >= 16)]
    assert not hit.empty
    assert hit.iloc[0].length_base_unit == 1
    assert hit.iloc[0].num_base_units == 12


def test_adjacent_distinct_repeats_are_both_reported():
    # Regression test: two genuinely distinct repeats sitting immediately
    # next to each other (no gap) used to have the lower-scoring one
    # silently discarded, because the overlap check treated touching-but-
    # non-overlapping ranges as overlapping (see test_detection_overlap.py).
    seq = "A" * 950 + "T" * 300 + "A" * 950
    df = find_slippage_sites(seq)
    assert df.shape[0] == 3

    t_run = df[df.sequence.str.contains("TTT")]
    assert not t_run.empty
    assert t_run.iloc[0].num_base_units == 300


def test_modify_df_slippage_splits_into_alternating_units():
    seq = "ATGCTAAT" + "GCGCGCGC" + "TTAGGCATGCCTAGC"
    df = find_slippage_sites(seq)
    matches = df[df.sequence.str.contains("GCGCGCGC")]

    df_modified = modify_df_slippage(matches)
    assert not df_modified.empty
    assert set(df_modified.sequence.unique()) == {"GC"}
