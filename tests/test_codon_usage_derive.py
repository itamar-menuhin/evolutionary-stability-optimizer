"""Tests for eso.codon_usage.derive_table_from_genome and
detect_genetic_code_num_from_gff - the ENc-based, real Sharp & Li (1987)
codon-adaptation-table derivation described in eso/codon_usage.py's own
docstring. All fixtures are small, synthetic, local FASTA/GFF content - no
network access.
"""

import pytest

from eso.codon_usage import (
    CustomCodonTableFileError,
    derive_table_from_genome,
    detect_genetic_code_num_from_gff,
)

_GENETIC_CODE_NUM = 11  # bacteria/archaea/plant plastid


def _write_fasta(tmp_path, records):
    """records: list of (header, sequence) pairs."""
    path = tmp_path / "cds.fna"
    with open(path, "w") as handle:
        for header, seq in records:
            handle.write(f">{header}\n{seq}\n")
    return str(path)


def _biased_gene(n_units):
    # Lys (AAA only) + Glu (GAA only), both twofold-degenerate amino acids
    # with the OTHER synonymous codon (AAG/GAG) never used - maximally
    # biased for these two degeneracy-2 groups.
    return ("AAAGAA" * n_units) + "TAA"


def _unbiased_gene(n_units):
    # Cycles Lys's two codons and Glu's two codons evenly - unbiased for
    # the same two amino acids, so this should rank with a much higher
    # (less biased) ENc than _biased_gene.
    return ("AAAAAGGAAGAG" * n_units) + "TAA"


def test_reference_set_picks_the_biased_gene_over_the_unbiased_one(tmp_path):
    fasta = _write_fasta(tmp_path, [
        ("biased", _biased_gene(60)),      # 120 codons + stop, well above the length floor
        ("unbiased", _unbiased_gene(60)),  # 240 codons + stop
    ])
    table = derive_table_from_genome(
        fasta, genetic_code_num=_GENETIC_CODE_NUM,
        min_len_codons=100, top_perc=0.5, min_gene_count=1,
    )
    assert table['K'] == {'AAA': 1.0, 'AAG': 0.0}
    assert table['E'] == {'GAA': 1.0, 'GAG': 0.0}


def test_short_gene_below_the_length_floor_is_excluded_regardless_of_bias(tmp_path):
    # Regression test for the real contamination this session found in
    # STABLES' own create_he.py (no length floor at all): a short,
    # extremely-biased gene (mimicking a leader peptide/toxin fragment)
    # must not end up in the reference set just because it's the most
    # biased sequence present - only min_len_codons=100+ genes are
    # eligible, however biased a shorter one is.
    long_unbiased = _unbiased_gene(60)  # 240 codons + stop - the only eligible gene
    short_biased = _biased_gene(10)     # 20 codons + stop - below the length floor

    fasta_without_short = _write_fasta(tmp_path, [("long", long_unbiased)])
    table_without_short = derive_table_from_genome(
        fasta_without_short, genetic_code_num=_GENETIC_CODE_NUM,
        min_len_codons=100, top_perc=1.0, min_gene_count=1,
    )

    fasta_with_short = tmp_path / "with_short.fna"
    fasta_with_short.write_text(f">long\n{long_unbiased}\n>short\n{short_biased}\n")
    table_with_short = derive_table_from_genome(
        str(fasta_with_short), genetic_code_num=_GENETIC_CODE_NUM,
        min_len_codons=100, top_perc=1.0, min_gene_count=1,
    )
    assert table_with_short == table_without_short


def test_min_gene_count_floor_widens_a_too_small_top_perc_selection(tmp_path):
    fasta = _write_fasta(tmp_path, [
        (f"gene{i}", _unbiased_gene(40)) for i in range(10)
    ])
    # top_perc=0.05 of 10 genes rounds to 0 - min_gene_count should still
    # guarantee at least 1 gene is used, not an empty/degenerate table.
    table = derive_table_from_genome(
        fasta, genetic_code_num=_GENETIC_CODE_NUM, min_len_codons=100, top_perc=0.05, min_gene_count=1,
    )
    assert table['K']['AAA'] > 0 or table['K']['AAG'] > 0


def test_no_gene_meets_the_length_floor_gives_a_friendly_message(tmp_path):
    fasta = _write_fasta(tmp_path, [("short", _biased_gene(5))])
    with pytest.raises(CustomCodonTableFileError, match="No gene"):
        derive_table_from_genome(fasta, genetic_code_num=_GENETIC_CODE_NUM, min_len_codons=100)


def test_genetic_code_num_is_required():
    with pytest.raises(CustomCodonTableFileError, match="genetic_code_num"):
        derive_table_from_genome("irrelevant.fna", genetic_code_num=None)


def test_stop_codon_default_is_present_when_not_derivable_from_the_reference_set(tmp_path):
    fasta = _write_fasta(tmp_path, [("gene", _unbiased_gene(60))])
    table = derive_table_from_genome(fasta, genetic_code_num=_GENETIC_CODE_NUM, min_len_codons=100, top_perc=1.0, min_gene_count=1)
    assert table['*'] == {'TAA': 0.33, 'TAG': 0.33, 'TGA': 0.34}


def test_detect_genetic_code_num_from_gff_majority_votes_transl_table(tmp_path):
    gff_path = tmp_path / "genomic.gff"
    gff_path.write_text(
        "##gff-version 3\n"
        "NC_1\tRefSeq\tCDS\t1\t100\t.\t+\t0\tID=cds-1;transl_table=11\n"
        "NC_1\tRefSeq\tCDS\t200\t300\t.\t+\t0\tID=cds-2;transl_table=11\n"
        "NC_1\tRefSeq\tCDS\t400\t500\t.\t+\t0\tID=cds-3;transl_table=4\n"
    )
    assert detect_genetic_code_num_from_gff(str(gff_path)) == 11


def test_detect_genetic_code_num_from_gff_raises_a_friendly_message_when_absent(tmp_path):
    gff_path = tmp_path / "genomic.gff"
    gff_path.write_text("##gff-version 3\nNC_1\tRefSeq\tgene\t1\t100\t.\t+\t.\tID=gene-1\n")
    with pytest.raises(CustomCodonTableFileError, match="transl_table"):
        detect_genetic_code_num_from_gff(str(gff_path))
