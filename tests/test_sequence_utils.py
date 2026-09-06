import pytest

from eso.sequence_utils import (
    InvalidSequenceError,
    first_invalid_base,
    parse_region,
    reverse_complement_seq,
    validate_dna_alphabet,
)


def test_reverse_complement_seq():
    assert reverse_complement_seq("ACGT") == "ACGT"
    assert reverse_complement_seq("AATTCCGG") == "CCGGAATT"


def test_parse_region_empty():
    assert parse_region('') == ()
    assert parse_region('None') == ()


def test_parse_region_single():
    assert parse_region('1-9') == [(0, 9)]


def test_parse_region_multiple():
    assert parse_region('1-9, 21-29') == [(0, 9), (20, 29)]


def test_parse_region_malformed():
    assert parse_region('not-a-region') == 'error'


# --- first_invalid_base / validate_dna_alphabet ---------------------------

def test_first_invalid_base_finds_the_position_and_character():
    assert first_invalid_base("ATGNTT") == (3, "N")


def test_first_invalid_base_none_for_valid_sequence():
    assert first_invalid_base("ATGCGTACGT") is None


def test_first_invalid_base_is_case_insensitive():
    assert first_invalid_base("atgcgt") is None
    assert first_invalid_base("atgRgt") == (3, "R")


def test_validate_dna_alphabet_raises_with_a_plain_english_message():
    with pytest.raises(InvalidSequenceError, match=r"position 4 \('N'\)"):
        validate_dna_alphabet("ATGNTT")


def test_validate_dna_alphabet_accepts_valid_sequence():
    validate_dna_alphabet("ATGCGTACGT")  # must not raise


def test_validate_dna_alphabet_uses_the_given_label():
    with pytest.raises(InvalidSequenceError, match="^'my_gene' sequence 0 contains"):
        validate_dna_alphabet("ATGNTT", sequence_label="'my_gene' sequence 0")


def test_validate_dna_alphabet_rejects_empty_sequence():
    # Regression test for a real bug: an empty sequence reached DNAChisel's
    # EnforceTranslation with a zero-length ORF region and crashed with a
    # bare `IndexError: string index out of range` deep inside DNAChisel -
    # confirmed directly before this fix, via
    # eso.optimize.optimization_engine("").
    with pytest.raises(InvalidSequenceError, match="is empty"):
        validate_dna_alphabet("")
