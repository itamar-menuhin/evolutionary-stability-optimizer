"""FASTA/GenBank file discovery and input validation."""

import gzip
import os
from glob import glob
from os import path
from pathlib import Path

from Bio import SeqIO
from dnachisel import biotools

from eso.sequence_utils import parse_region

FILE_ENDINGS = {
    'fasta': ['fasta', 'fna', 'ffn', 'faa', 'frn', 'fa'],
    'genbank': ['genbank', 'gb', 'gbk'],
}


def file_opener(file):
    """`file` is a (filepath, filetype) tuple, filetype in {'fasta', 'genbank'}."""
    filepath, filetype = file
    if filepath.endswith('.gz'):
        with gzip.open(filepath, "rt") as handle:
            return list(SeqIO.parse(handle, filetype))
    with open(filepath, "rt") as handle:
        return list(SeqIO.parse(handle, filetype))


def relevant_file_paths(input_folder=None):
    """Find all FASTA/GenBank files (optionally gzipped) directly in or one level
    under `input_folder`, returned as a list of (filepath, filetype) tuples.
    """
    if input_folder is None:
        input_folder = os.getcwd()
    input_as_path = Path(input_folder)

    files = []
    for filetype, endings in FILE_ENDINGS.items():
        for ending in endings:
            for gzipping in ['', '.gz']:
                pattern = ending + gzipping
                files.extend(
                    (filepath, filetype) for filepath in glob(path.join(input_as_path, '*', f'*.{pattern}'))
                )
                files.extend(
                    (filepath, filetype) for filepath in glob(path.join(input_as_path, f'*.{pattern}'))
                )
    return files


def exclusion_gc_tester(file, indexes):
    """Given a file's exclusion (locked) regions, compute the range of GC-content
    limits that can still legally be enforced across a sliding window that
    partially overlaps those regions.
    """
    data = file_opener(file)
    filename_indexes = path.basename(file[0]).split('.')[0]

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
                curr_gc = biotools.gc_content(overlap)
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
                curr_gc = biotools.gc_content(overlap)
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
                for curr_ind in orf_indexes:
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
