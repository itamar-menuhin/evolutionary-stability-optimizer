"""Tests for eso.io_utils - previously untested.

Includes a regression test for a real bug found while reviewing: exclusion
region validation in test_input looped over `orf_indexes` (already validated
above it) instead of the just-parsed `region_indexes`, so a malformed
exclusion region (e.g. start >= end, or an unparseable string) was never
actually checked - it silently passed validation and only surfaced later as
a confusing, unrelated GC-limit error message.
"""

import gzip
import os

import pytest

from eso.io_utils import IndexesFileError, file_opener, load_indexes_from_file, relevant_file_paths
from eso.io_utils import test_input as validate_input


# --- test_input: GC bounds -----------------------------------------------

def test_negative_mini_gc_rejected():
    assert validate_input(-0.1, 0.7, {}, []) == 'The minimal GC content must be at least 0!'


def test_maxi_gc_above_one_rejected():
    assert validate_input(0.3, 1.1, {}, []) == 'The maximal GC content must be no more than 1!'


def test_mini_gc_not_less_than_maxi_gc_rejected():
    assert validate_input(0.7, 0.7, {}, []) == 'The minimal GC content must be less than the maximum!'


def test_valid_gc_bounds_with_no_indexes_succeeds():
    assert validate_input(0.3, 0.7, {}, []) == 'Success!'


# --- test_input: ORF region validation -----------------------------------

def test_malformed_orf_region_rejected():
    result = validate_input(0.3, 0.7, {('f', '0'): ('not-a-region', '')}, [])
    assert 'ORF regions must be formatted' in result


def test_orf_region_not_divisible_by_three_rejected():
    result = validate_input(0.3, 0.7, {('f', '0'): ('1-9,10-14', '')}, [])  # 4nt second region
    assert 'divisible by 3' in result


def test_orf_start_not_before_end_rejected():
    result = validate_input(0.3, 0.7, {('f', '0'): ('9-3', '')}, [])
    assert 'Start index must be smaller than end index!' == result


def test_valid_orf_region_with_no_exclusion_succeeds(tmp_path):
    # indexes being non-empty always triggers exclusion_gc_tester over `files`,
    # even with no exclusion region given, so a real file matching the index
    # key's filename is required for the success path to be reachable at all.
    file_path = tmp_path / 'f.fasta'
    file_path.write_text('>x\n' + 'ACGT' * 15 + '\n')

    result = validate_input(0.3, 0.7, {('f', '0'): ('1-9', '')}, [(str(file_path), 'fasta')])
    assert result == 'Success!'


# --- test_input: exclusion region validation (regression) ---------------

def test_malformed_exclusion_region_is_actually_rejected():
    # regression test: this used to fall through to a confusing GC-limit
    # error instead of being caught here, because the validation loop
    # checked orf_indexes instead of region_indexes.
    result = validate_input(0.3, 0.7, {('f', '0'): ('1-9', 'not-a-region')}, [])
    assert 'Exclusion regions must be formatted' in result


def test_exclusion_region_start_not_before_end_is_actually_rejected():
    result = validate_input(0.3, 0.7, {('f', '0'): ('1-9', '10-5')}, [])
    assert result == 'Start index must be smaller than end index, also for exclusion sites!'


def test_valid_orf_and_exclusion_region_together_succeeds(tmp_path):
    file_path = tmp_path / 'f.fasta'
    file_path.write_text('>x\n' + 'ACGT' * 15 + '\n')

    result = validate_input(0.3, 0.7, {('f', '0'): ('1-9', '1-6')}, [(str(file_path), 'fasta')])
    assert result == 'Success!'


# --- relevant_file_paths ---------------------------------------------------

def test_finds_files_directly_in_and_one_level_under_input_folder(tmp_path):
    (tmp_path / 'a.fasta').write_text('>x\nACGT\n')
    (tmp_path / 'b.gb').write_text('')
    sub = tmp_path / 'sub'
    sub.mkdir()
    (sub / 'c.fna').write_text('>y\nTTTT\n')
    with gzip.open(tmp_path / 'd.fasta.gz', 'wt') as f:
        f.write('>z\nGGGG\n')

    files = relevant_file_paths(str(tmp_path))
    found = {(os.path.basename(f), t) for f, t in files}

    assert ('a.fasta', 'fasta') in found
    assert ('b.gb', 'genbank') in found
    assert ('c.fna', 'fasta') in found
    assert ('d.fasta.gz', 'fasta') in found


def test_ignores_unrelated_file_extensions(tmp_path):
    (tmp_path / 'notes.txt').write_text('hello')
    files = relevant_file_paths(str(tmp_path))
    assert files == []


def test_empty_folder_returns_no_files(tmp_path):
    assert relevant_file_paths(str(tmp_path)) == []


# --- file_opener -----------------------------------------------------------

def test_opens_plain_fasta_file(tmp_path):
    file_path = tmp_path / 'a.fasta'
    file_path.write_text('>x\nACGT\n')

    records = file_opener((str(file_path), 'fasta'))

    assert len(records) == 1
    assert str(records[0].seq) == 'ACGT'


def test_opens_gzipped_fasta_file(tmp_path):
    file_path = tmp_path / 'a.fasta.gz'
    with gzip.open(file_path, 'wt') as f:
        f.write('>x\nTTTT\n')

    records = file_opener((str(file_path), 'fasta'))

    assert len(records) == 1
    assert str(records[0].seq) == 'TTTT'


def test_opens_multi_record_fasta_file(tmp_path):
    file_path = tmp_path / 'a.fasta'
    file_path.write_text('>x\nACGT\n>y\nTTTT\n')

    records = file_opener((str(file_path), 'fasta'))

    assert [str(r.seq) for r in records] == ['ACGT', 'TTTT']


# --- load_indexes_from_file ------------------------------------------------

def test_load_indexes_from_file_builds_the_expected_dict(tmp_path):
    file_path = tmp_path / 'indexes.json'
    file_path.write_text(
        '[{"file": "my_gene", "seq_index": "0", "orf_regions": "1-6, 51-68", '
        '"exclusion_regions": "1-6, 50-68"}]'
    )

    indexes = load_indexes_from_file(str(file_path))

    assert indexes == {("my_gene", "0"): ("1-6, 51-68", "1-6, 50-68")}


def test_load_indexes_from_file_defaults_missing_exclusion_regions_to_empty(tmp_path):
    file_path = tmp_path / 'indexes.json'
    file_path.write_text('[{"file": "my_gene", "seq_index": "0", "orf_regions": "1-6"}]')

    indexes = load_indexes_from_file(str(file_path))

    assert indexes == {("my_gene", "0"): ("1-6", "")}


def test_load_indexes_from_file_missing_file_raises_friendly_error(tmp_path):
    with pytest.raises(IndexesFileError, match="Can't find the indexes file"):
        load_indexes_from_file(str(tmp_path / 'does_not_exist.json'))


def test_load_indexes_from_file_invalid_json_raises_friendly_error(tmp_path):
    file_path = tmp_path / 'indexes.json'
    file_path.write_text('{not valid json')

    with pytest.raises(IndexesFileError, match="isn't valid JSON"):
        load_indexes_from_file(str(file_path))


def test_load_indexes_from_file_non_list_raises_friendly_error(tmp_path):
    file_path = tmp_path / 'indexes.json'
    file_path.write_text('{"file": "my_gene"}')

    with pytest.raises(IndexesFileError, match="must contain a JSON list"):
        load_indexes_from_file(str(file_path))


def test_load_indexes_from_file_entry_missing_required_key_raises_friendly_error(tmp_path):
    file_path = tmp_path / 'indexes.json'
    file_path.write_text('[{"file": "my_gene"}]')

    with pytest.raises(IndexesFileError, match='"file" and "seq_index" keys'):
        load_indexes_from_file(str(file_path))
