from eso.detection.junction_linker import (
    find_slippage_1,
    find_slippage_l,
    find_recombination,
    find_suspect,
)


def test_find_slippage_1_detects_homopolymer_run_at_junction():
    # target tail "TGA" + linker "AAA" + host head "TGC" = "TGAAAATGC", contains "AAAA"
    result = find_slippage_1("ATGA", "AAA", "TGCA")
    assert result == "AAAA"


def test_find_slippage_1_no_run():
    result = find_slippage_1("ATGC", "CT", "GATC")
    assert result == ""


def test_find_slippage_l_detects_dinucleotide_repeat_at_junction():
    # junction forms "...GCGCGC..." (GC x3) spanning target tail + linker + host head
    result = find_slippage_l("AAAAAGC", "GC", "GCAAAAA", length=2)
    assert result == "GCGCGC"


def test_find_recombination_detects_junction_kmer_reused_in_target():
    # junction forms "GATC", which already appears at the start of target_seq
    result = find_recombination("GATCAAAAAA", "G", "ATCTTTTTTT", length=4, mode="linker")
    assert result == "GATC"


def test_find_recombination_no_shared_kmer():
    result = find_recombination("ATGCATGC", "TT", "CCGGCCGG", length=4, mode="linker")
    assert result == ""


def test_find_suspect_aggregates_all_checks():
    # combine a homopolymer-run junction with an otherwise clean target/host
    result = find_suspect("ATGA", "AAA", "TGCA", max_l=4)
    assert "AAAA" in result


def test_find_suspect_empty_for_clean_junction():
    result = find_suspect("ATGCATGCATGCATGC", "TT", "CCGGCCGGCCGGCCGG", max_l=4)
    assert result == ""
