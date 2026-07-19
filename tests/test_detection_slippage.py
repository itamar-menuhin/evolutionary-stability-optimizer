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


def test_modify_df_slippage_splits_into_alternating_units():
    seq = "ATGCTAAT" + "GCGCGCGC" + "TTAGGCATGCCTAGC"
    df = find_slippage_sites(seq)
    matches = df[df.sequence.str.contains("GCGCGCGC")]

    df_modified = modify_df_slippage(matches)
    assert not df_modified.empty
    assert set(df_modified.sequence.unique()) == {"GC"}
