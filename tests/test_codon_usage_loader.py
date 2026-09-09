"""Tests for eso.codon_usage.load_custom_codon_table_from_file - the
file-based loader behind `eso-optimize --codon-usage-table-file`, aimed at
someone bringing in a project-specific codon usage table rather than a
bioinformatics/pandas expert.
"""

import pytest

from eso.codon_usage import CustomCodonTableFileError, load_custom_codon_table_from_file

_VALID_CSV = (
    'codon,aa,freq_within_aa\n'
    'GCT,A,0.5\nGCC,A,0.5\n'
    'AAA,K,0.4\nAAG,K,0.6\n'
    'TAA,*,0.5\nTAG,*,0.5\n'
)


def _write(tmp_path, name, content):
    file_path = tmp_path / name
    file_path.write_text(content)
    return str(file_path)


def test_loads_a_well_formed_table(tmp_path):
    file_path = _write(tmp_path, "good.csv", _VALID_CSV)
    table = load_custom_codon_table_from_file(file_path)
    assert table['A'] == {'GCT': 0.5, 'GCC': 0.5}
    assert table['K'] == {'AAA': 0.4, 'AAG': 0.6}


def test_missing_file_gives_friendly_message():
    with pytest.raises(CustomCodonTableFileError, match="Can't find"):
        load_custom_codon_table_from_file("no_such_file_here.csv")


def test_not_actually_a_csv_gives_friendly_message(tmp_path):
    # A single unterminated quote is one of the few things pandas' C parser
    # itself refuses outright, rather than silently absorbing into a field.
    file_path = _write(tmp_path, "bad.csv", 'codon,aa,freq_within_aa\n"ATG,M,1\n')
    with pytest.raises(CustomCodonTableFileError, match="could not be read as a CSV"):
        load_custom_codon_table_from_file(file_path)


def test_missing_columns_gives_friendly_message(tmp_path):
    file_path = _write(tmp_path, "missing_col.csv", 'codon,aa\nATG,M\n')
    with pytest.raises(CustomCodonTableFileError, match="missing the column"):
        load_custom_codon_table_from_file(file_path)


def test_empty_table_gives_friendly_message(tmp_path):
    file_path = _write(tmp_path, "empty.csv", 'codon,aa,freq_within_aa\n')
    with pytest.raises(CustomCodonTableFileError, match="no data rows"):
        load_custom_codon_table_from_file(file_path)


def test_invalid_codon_gives_friendly_message(tmp_path):
    file_path = _write(tmp_path, "bad_codon.csv", 'codon,aa,freq_within_aa\nATGX,M,1\n')
    with pytest.raises(CustomCodonTableFileError, match="invalid entries in its 'codon' column"):
        load_custom_codon_table_from_file(file_path)


def test_non_numeric_frequency_gives_friendly_message(tmp_path):
    file_path = _write(tmp_path, "bad_freq.csv", 'codon,aa,freq_within_aa\nATG,M,not_a_number\n')
    with pytest.raises(CustomCodonTableFileError, match="non-numeric value"):
        load_custom_codon_table_from_file(file_path)


def test_missing_stop_codon_falls_back_to_default_distribution(tmp_path):
    # Matches eso's own bundled CSV-backed tables (see test_codon_usage.py's
    # test_csv_backed_tables_fall_back_to_a_default_stop_codon_distribution) -
    # a table with no stop-codon row shouldn't fail, just fall back. Also
    # confirms a table need NOT cover every amino acid to be considered
    # well-formed (only the amino acids a given real sequence actually uses
    # matter, which this loader has no way to know in advance).
    file_path = _write(tmp_path, "no_stop.csv", 'codon,aa,freq_within_aa\nGCT,A,0.5\nGCC,A,0.5\n')
    table = load_custom_codon_table_from_file(file_path)
    assert table['*'] == {'TAA': 0.33, 'TAG': 0.33, 'TGA': 0.34}


def test_codon_filed_under_the_wrong_amino_acid_gives_friendly_message(tmp_path):
    # GCT translates to Alanine ('A'), not Lysine ('K') - a realistic
    # copy-paste/column-order mistake, checked directly rather than only
    # surfacing as a confusing failure deep inside DNAChisel later.
    file_path = _write(tmp_path, "mislabeled.csv", 'codon,aa,freq_within_aa\nGCT,K,0.5\nAAA,K,0.5\n')
    with pytest.raises(CustomCodonTableFileError, match="filed under the wrong amino acid"):
        load_custom_codon_table_from_file(file_path)
