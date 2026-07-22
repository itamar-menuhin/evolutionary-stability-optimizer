"""FASTA/GenBank file discovery and input validation."""

import gzip
import json
import os
from os import path
from pathlib import Path

from Bio import SeqIO
from dnachisel import biotools

from eso.sequence_utils import parse_region

FILE_ENDINGS = {
    'fasta': ['fasta', 'fna', 'ffn', 'faa', 'frn', 'fa'],
    'genbank': ['genbank', 'gb', 'gbk'],
}


def file_stem(filepath):
    """Basename with its extension(s) removed, preserving any other dots in
    the stem itself. `path.basename(filepath).split('.')[0]` (used to
    determine this) truncated at the FIRST dot, so a filename like
    "sample.2024.fasta" silently became the stem "sample" instead of
    "sample.2024" - meaning an `indexes` entry keyed on the intended full
    stem would never match, silently falling back to no ORF/exclusion
    regions instead of erroring.
    """
    name = path.basename(filepath)
    if name.endswith('.gz'):
        name = name[:-3]
    return path.splitext(name)[0]


class IndexesFileError(Exception):
    """An --indexes-file failed to load or was malformed.

    Raised with a plain-English message aimed at someone filling in ORF and
    exclusion regions for their sequences, not a Python or JSON expert - the
    message alone should be enough to fix the problem.
    """


def load_indexes_from_file(file_path):
    """Load the `indexes` argument to eso.pipeline.main from a JSON file, for
    use as the CLI's `--indexes-file`.

    The file must be a JSON list of objects, each with:
    - "file": the FASTA/GenBank file's stem (no extension - matches
      eso.io_utils.file_stem), e.g. "my_gene" for "my_gene.fasta".
    - "seq_index": which record within that file, as a string, 0-indexed in
      file order, e.g. "0" for the first sequence.
    - "orf_regions": 1-indexed, inclusive region string, e.g. "1-6, 51-68".
    - "exclusion_regions": same format; omit or use "" for no exclusions.

    Region strings themselves are validated later, the same way as when
    `indexes` is passed directly to eso.pipeline.main - this only checks the
    file's own structure (valid JSON, a list, well-formed entries).

    Returns
    -------
    dict mapping (file_stem, seq_index) -> (orf_regions, exclusion_regions),
    matching eso.pipeline.main's `indexes` parameter.
    """
    if not path.isfile(file_path):
        raise IndexesFileError(
            f"Can't find the indexes file '{file_path}'. Check the path is correct.")

    with open(file_path, "r", encoding="utf-8") as handle:
        try:
            entries = json.load(handle)
        except json.JSONDecodeError as e:
            raise IndexesFileError(
                f"'{file_path}' isn't valid JSON: {e}. Check for missing commas, quotes, or brackets."
            ) from e

    if not isinstance(entries, list):
        raise IndexesFileError(
            f"'{file_path}' must contain a JSON list of entries (got a {type(entries).__name__} instead). "
            'Each entry looks like {"file": "my_gene", "seq_index": "0", "orf_regions": "1-6, 51-68", '
            '"exclusion_regions": "1-6, 50-68"}.')

    indexes = {}
    for ii, entry in enumerate(entries):
        if not isinstance(entry, dict) or 'file' not in entry or 'seq_index' not in entry:
            raise IndexesFileError(
                f"Entry {ii} of '{file_path}' must be an object with at least \"file\" and "
                f'"seq_index" keys - got {entry!r}.')

        key = (str(entry['file']), str(entry['seq_index']))
        orf_regions = str(entry.get('orf_regions', ''))
        exclusion_regions = str(entry.get('exclusion_regions', ''))
        indexes[key] = (orf_regions, exclusion_regions)

    return indexes


def file_opener(file):
    """`file` is a (filepath, filetype) tuple, filetype in {'fasta', 'genbank'}."""
    filepath, filetype = file
    if filepath.endswith('.gz'):
        with gzip.open(filepath, "rt", encoding="utf-8") as handle:
            return list(SeqIO.parse(handle, filetype))
    with open(filepath, "rt", encoding="utf-8") as handle:
        return list(SeqIO.parse(handle, filetype))


def _matching_filetype(filename):
    """Return the FILE_ENDINGS key `filename` belongs to (case-insensitively,
    ignoring a trailing .gz), or None if it doesn't match any known extension.

    Matching case-insensitively (rather than relying on glob('*.fasta')) matters
    because glob's case sensitivity is filesystem-dependent - Windows filesystems
    are case-insensitive so `*.fasta` there also matches `GENE.FASTA`, but the
    same glob call on a (case-sensitive) Mac/Linux filesystem would silently skip
    it. Confirmed this is a real, reachable discrepancy, not hypothetical -
    without this, a file discovery result depends on which OS the tool happens
    to run on for the exact same input folder.
    """
    name = filename.lower()
    if name.endswith('.gz'):
        name = name[:-3]
    if '.' not in name:
        return None
    ext = name.rsplit('.', 1)[-1]
    for filetype, endings in FILE_ENDINGS.items():
        if ext in endings:
            return filetype
    return None


def relevant_file_paths(input_folder=None):
    """Find all FASTA/GenBank files (optionally gzipped) directly in or one level
    under `input_folder`, returned as a list of (filepath, filetype) tuples.
    """
    if input_folder is None:
        input_folder = os.getcwd()
    input_as_path = Path(input_folder)

    candidates = [p for p in input_as_path.glob('*') if p.is_file()]
    for entry in input_as_path.glob('*'):
        if entry.is_dir():
            candidates.extend(p for p in entry.glob('*') if p.is_file())

    files = []
    for candidate in candidates:
        filetype = _matching_filetype(candidate.name)
        if filetype is not None:
            files.append((str(candidate), filetype))
    return files


def exclusion_gc_tester(file, indexes):
    """Given a file's exclusion (locked) regions, compute the range of GC-content
    limits that can still legally be enforced across a sliding window that
    partially overlaps those regions.
    """
    data = file_opener(file)
    filename_indexes = file_stem(file[0])

    legal_mini_gc = 1.0
    legal_maxi_gc = 0.0

    for ii, record in enumerate(data):
        example_seq = str(record.seq).upper()
        windowed_gc_content = biotools.gc_content(example_seq, window_size=50)
        seq_indexes = str(ii)

        if (filename_indexes, seq_indexes) not in indexes:
            exclusion_regions = ()
        else:
            relevant_index_data = indexes[(filename_indexes, seq_indexes)]
            exclusion_regions = parse_region(relevant_index_data[1])

        for region in exclusion_regions:
            start_ind, end_ind = region
            curr_min_gc, curr_max_gc = 1.0, 0.0

            if start_ind + 50 <= end_ind:
                curr_max_gc = max(windowed_gc_content[start_ind:end_ind - 50])
                curr_min_gc = min(windowed_gc_content[start_ind:end_ind - 50])

            for jj in range(50):
                if start_ind - jj < 0:
                    break
                curr_region = (start_ind - jj, start_ind + 49 - jj)
                overlap = ''
                for ex_reg in exclusion_regions:
                    overlap_reg = biotools.windows_overlap(curr_region, ex_reg)
                    if overlap_reg is not None:
                        overlap += example_seq[overlap_reg[0]:overlap_reg[1]]
                # biotools.gc_content('') is 0/0 (NaN, with a RuntimeWarning) -
                # `overlap` is empty whenever none of exclusion_regions actually
                # overlaps this particular sliding window, which does happen (e.g.
                # a window entirely before the region under test). The GC content
                # of an empty overlap contributes nothing regardless, so skip the
                # call rather than asking for a GC content that doesn't exist.
                curr_gc = biotools.gc_content(overlap) if overlap else 0.0
                curr_max_gc = max(curr_max_gc, curr_gc * len(overlap) / 50.0)
                curr_min_gc = min(curr_min_gc, curr_gc * len(overlap) / 50.0 + (1.0 - len(overlap) / 50.0))

            for jj in range(50):
                if end_ind + jj >= len(example_seq):
                    break
                curr_region = (end_ind - 49 + jj, end_ind + jj)
                overlap = ''
                for ex_reg in exclusion_regions:
                    overlap_reg = biotools.windows_overlap(curr_region, ex_reg)
                    if overlap_reg is not None:
                        overlap += example_seq[overlap_reg[0]:overlap_reg[1]]
                curr_gc = biotools.gc_content(overlap) if overlap else 0.0
                curr_max_gc = max(curr_max_gc, curr_gc * len(overlap) / 50.0)
                curr_min_gc = min(curr_min_gc, curr_gc * len(overlap) / 50.0 + (1.0 - len(overlap) / 50.0))

            legal_mini_gc = min(curr_min_gc, legal_mini_gc)
            legal_maxi_gc = max(curr_max_gc, legal_maxi_gc)

    return legal_mini_gc, legal_maxi_gc


def test_input(mini_gc, maxi_gc, indexes, files):
    """Validate GC-content bounds and ORF/exclusion region formatting before running.

    `indexes` maps (filename, seq_index) -> (orf_region_string, exclusion_region_string),
    e.g. {"myfile": ("1-9, 21-29", "")} - 1-indexed, inclusive, as entered by a biologist.
    """
    if mini_gc < 0.0:
        return 'The minimal GC content must be at least 0!'
    if maxi_gc > 1.0:
        return 'The maximal GC content must be no more than 1!'
    if mini_gc >= maxi_gc:
        return 'The minimal GC content must be less than the maximum!'

    if len(indexes) > 0:
        for index_value in indexes.values():
            orf_indexes = parse_region(index_value[0])
            if orf_indexes == 'error':
                return 'ORF regions must be formatted as "start_1-end_1,start_2-end_2,...", for example "1-9, 21-29"'
            for curr_ind in orf_indexes:
                if curr_ind[0] < 0:
                    return 'Start index must be greater than 0!'
                if curr_ind[0] >= curr_ind[1]:
                    return 'Start index must be smaller than end index!'
                if (curr_ind[1] - curr_ind[0]) % 3 != 0:
                    return 'Indexes must describe sequence length divisible by 3!'

            if index_value[1] not in ('', 'None'):
                region_indexes = parse_region(index_value[1])
                if region_indexes == 'error':
                    return 'Exclusion regions must be formatted as "start_1-end_1,start_2-end_2,...", for example "1-8, 50-103"'
                # was `for curr_ind in orf_indexes` - validated the already-checked
                # ORF indexes a second time and never actually checked the
                # exclusion regions themselves, so a malformed exclusion region
                # (e.g. start >= end) slipped through validation entirely and
                # only surfaced later as a confusing, unrelated GC-limit error.
                for curr_ind in region_indexes:
                    if curr_ind[0] < 0:
                        return 'Start index must be greater than 0, also for exclusion sites!!'
                    if curr_ind[0] >= curr_ind[1]:
                        return 'Start index must be smaller than end index, also for exclusion sites!'

        legal_mini_gc, legal_maxi_gc = 0.0, 1.0
        for file in files:
            curr_mini_gc, curr_maxi_gc = exclusion_gc_tester(file, indexes)
            legal_mini_gc = max(legal_mini_gc, curr_mini_gc) + 0.01
            legal_maxi_gc = min(legal_maxi_gc, curr_maxi_gc) - 0.01

        if legal_mini_gc < mini_gc or legal_maxi_gc > maxi_gc:
            return (
                f'Given the exclusion regions defined by the user, the maximal GC limit must be at least '
                f'{legal_maxi_gc} and the minimal GC limit must be at most {legal_mini_gc}. Note that the '
                f'maximal/minimal GC limit required might be a bit higher/lower to account also for ORF regions.'
            )

    return 'Success!'
