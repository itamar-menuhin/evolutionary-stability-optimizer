"""Tests for eso.codon_usage - previously untested.

Includes a regression test for a real bug found while reviewing:
cub_kompas() converted its 3-letter amino acid keys (including "END" for
stop codons) via Bio.SeqUtils.seq1, but seq1("END") returns "X" (undefined
amino acid), not "*" (stop) - so the kompas table's stop-codon frequencies
were silently filed under the wrong key and never used for stop-codon
scoring during CodonOptimize.
"""

import pytest
from Bio.Seq import Seq

from eso.codon_usage import CODON_USAGE_TABLES, cub_kompas

STANDARD_AMINO_ACIDS = set('ACDEFGHIKLMNPQRSTVWY') | {'*'}


def test_kompas_stop_codons_are_filed_under_star_not_x():
    table = cub_kompas()
    assert '*' in table
    assert table['*'] == {'TAA': 0.399841, 'TAG': 0.339428, 'TGA': 0.260731}
    assert 'X' not in table


@pytest.mark.parametrize('name', list(CODON_USAGE_TABLES))
def test_table_has_exactly_the_20_amino_acids_plus_stop(name):
    table = CODON_USAGE_TABLES[name]()
    assert set(table.keys()) == STANDARD_AMINO_ACIDS


@pytest.mark.parametrize('name', list(CODON_USAGE_TABLES))
def test_every_codon_translates_to_the_amino_acid_it_is_filed_under(name):
    table = CODON_USAGE_TABLES[name]()
    for aa, codons in table.items():
        for codon in codons:
            assert str(Seq(codon).translate()) == aa, f"{name}: {codon} filed under {aa}"


@pytest.mark.parametrize('name', list(CODON_USAGE_TABLES))
def test_frequencies_within_each_amino_acid_sum_to_one(name):
    table = CODON_USAGE_TABLES[name]()
    for aa, codons in table.items():
        assert sum(codons.values()) == pytest.approx(1.0, abs=1e-2), f"{name}: {aa} sums to {sum(codons.values())}"


def test_csv_backed_tables_fall_back_to_a_default_stop_codon_distribution():
    # the antibody CSVs don't include stop-codon rows at all - confirm the
    # documented fallback default is what actually gets used.
    from eso.codon_usage import cub_human_antibody_heavy_chain, cub_human_antibody_light_chain

    for table in (cub_human_antibody_heavy_chain(), cub_human_antibody_light_chain()):
        assert table['*'] == {'TAA': 0.33, 'TAG': 0.33, 'TGA': 0.34}
