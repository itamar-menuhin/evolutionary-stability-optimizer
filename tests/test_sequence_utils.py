import pandas as pd
import pytest

from eso.sequence_utils import (
    InvalidSequenceError,
    add_backward_sites,
    first_invalid_base,
    parse_region,
    reverse_complement_seq,
    shorten_sequences,
    validate_dna_alphabet,
)


def test_reverse_complement_seq():
    assert reverse_complement_seq("ACGT") == "ACGT"
    assert reverse_complement_seq("AATTCCGG") == "CCGGAATT"


def test_reverse_complement_seq_is_case_insensitive_and_uppercases_output():
    assert reverse_complement_seq("aattccgg") == "CCGGAATT"


def test_reverse_complement_seq_raises_a_clear_error_for_a_non_dna_character():
    # Regression test for a real gap: this used to leave every caller exposed
    # to a raw, unhelpful `KeyError` for even a single unexpected character -
    # confirmed directly before this fix (KeyError: 'N', no eso-level message
    # at all) - the same failure mode validate_dna_alphabet was written to
    # close off for a whole sequence, but this function had no equivalent.
    with pytest.raises(InvalidSequenceError, match=r"'N'"):
        reverse_complement_seq("ACGNT")


def test_parse_region_empty():
    assert parse_region('') == ()
    assert parse_region('None') == ()


def test_parse_region_single():
    assert parse_region('1-9') == [(0, 9)]


def test_parse_region_multiple():
    assert parse_region('1-9, 21-29') == [(0, 9), (20, 29)]


def test_parse_region_malformed():
    assert parse_region('not-a-region') == 'error'


def test_add_backward_sites_appends_reverse_complement_rows_unchanged_coordinates():
    df = pd.DataFrame([{"sequence": "AATTCCGG", "start": 3, "end": 11}])

    result = add_backward_sites(df)

    assert len(result) == 2
    assert set(result.sequence) == {"AATTCCGG", "CCGGAATT"}
    # coordinates are untouched by the reverse-complement - only the letters change
    assert (result.start == 3).all() and (result.end == 11).all()


def test_add_backward_sites_does_not_mutate_the_input():
    df = pd.DataFrame([{"sequence": "AATTCCGG", "start": 3, "end": 11}])
    add_backward_sites(df)
    assert df.sequence.iloc[0] == "AATTCCGG"


def test_add_backward_sites_result_has_a_unique_index():
    # Regression test for a real, hardened-regardless-of-reachability gap:
    # concatenating with ignore_index=False left the result with each
    # original row's index repeated twice (once per half) - harmless for
    # every current caller (none does label-based .loc[idx] access), but a
    # real trap for one that did, silently getting two rows back instead of
    # one. Confirmed directly before this fix.
    df = pd.DataFrame([{"sequence": "AATT", "start": 0, "end": 4}, {"sequence": "CCGG", "start": 10, "end": 14}])
    result = add_backward_sites(df)
    assert result.index.is_unique


def test_shorten_sequences_drops_the_last_base_and_decrements_end():
    df = pd.DataFrame([{"sequence": "ACGTAC", "start": 5, "end": 11}])

    result = shorten_sequences(df)

    assert result.sequence.iloc[0] == "ACGTA"
    assert result.start.iloc[0] == 5
    assert result.end.iloc[0] == 10


def test_shorten_sequences_does_not_mutate_the_input():
    df = pd.DataFrame([{"sequence": "ACGTAC", "start": 5, "end": 11}])
    shorten_sequences(df)
    assert df.sequence.iloc[0] == "ACGTAC" and df.end.iloc[0] == 11


def test_shorten_sequences_rejects_a_length_1_sequence():
    # Regression test for a real, hardened-regardless-of-reachability gap:
    # shortening a 1nt sequence used to silently produce an empty-string
    # "sequence" (x[:-1] on a 1-character string) rather than raising -
    # every current caller only ever passes 16-17nt windows, so this isn't
    # reachable today, but nothing enforced that, and a degenerate empty
    # pattern silently reaching detection/DNAChisel code downstream is worse
    # than failing loudly here instead.
    df = pd.DataFrame([{"sequence": "A", "start": 5, "end": 6}])
    with pytest.raises(ValueError, match="length >= 2"):
        shorten_sequences(df)


def test_shorten_sequences_handles_an_empty_dataframe():
    df = pd.DataFrame(columns=["sequence", "start", "end"])
    result = shorten_sequences(df)
    assert result.empty


def test_parse_region_rejects_an_extra_dash_instead_of_silently_truncating():
    # Regression test for a real, previously-silent bug: an extra dash (e.g.
    # a typo like "10-20-30") used to parse as (9, 20), silently dropping
    # "-30" instead of being flagged as malformed. Confirmed directly before
    # this fix.
    assert parse_region('10-20-30') == 'error'
    assert parse_region('1-9,20-30-40') == 'error'
    assert parse_region('10--5') == 'error'


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
