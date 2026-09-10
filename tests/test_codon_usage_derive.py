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


def test_genetic_code_with_a_different_degeneracy_structure_still_works_correctly(tmp_path):
    # Regression test for a real gap, confirmed directly (not assumed): the
    # ENc formula _calc_enc uses had fixed coefficients (2, 9, 1, 5, 3) that
    # are the standard genetic code's own count of amino-acid families at
    # each degeneracy level - true for tables 1 and 11, but genuinely
    # different for every other NCBI genetic code table (e.g. table 2,
    # vertebrate mitochondrial, has no 3-fold-degenerate amino acid at all -
    # AGA/AGG are stop codons there, not Arg; ATA is Met, not Ile). Fixed by
    # computing those coefficients directly from whichever genetic code is
    # actually given (_degeneracy_family_counts), rather than assuming the
    # standard code's own numbers or rejecting anything else outright.
    # Table 2's reassignments don't touch Lys/Glu's own codons, so the same
    # biased-vs-unbiased fixture used for the standard code (table 11)
    # elsewhere in this file still exercises the real selection logic here.
    fasta = _write_fasta(tmp_path, [
        ("biased", _biased_gene(60)),
        ("unbiased", _unbiased_gene(60)),
    ])
    table = derive_table_from_genome(
        fasta, genetic_code_num=2,  # vertebrate mitochondrial
        min_len_codons=100, top_perc=0.5, min_gene_count=1,
    )
    assert table['K'] == {'AAA': 1.0, 'AAG': 0.0}
    assert table['E'] == {'GAA': 1.0, 'GAG': 0.0}


def test_calc_enc_matches_the_standard_formula_for_the_standard_genetic_code(tmp_path):
    # Direct check that the generalized, any-genetic-code formula
    # (_degeneracy_family_counts + _calc_enc) reproduces the exact original
    # hardcoded standard-code formula (2 + 9/d2 + 1/d3 + 5/d4 + 3/d6) for
    # tables 1/11 - not just "doesn't crash", the actual number.
    from Bio.Data.CodonTable import unambiguous_dna_by_id
    from collections import defaultdict
    from eso.codon_usage import _calc_enc, _degeneracy_family_counts

    aa_to_codons = defaultdict(list)
    for codon, aa in unambiguous_dna_by_id[11].forward_table.items():
        aa_to_codons[aa].append(codon)
    degeneracy_family_counts = _degeneracy_family_counts(aa_to_codons)
    assert degeneracy_family_counts == {1: 2, 2: 9, 3: 1, 4: 5, 6: 3}

    codons = (["AAA"] * 8 + ["AAG"] * 2       # Lys, 2-fold, biased
              + ["ATA"] * 5 + ["ATC"] * 3 + ["ATT"] * 2)  # Ile, 3-fold

    d2_group = (8 / 10) ** 2 + (2 / 10) ** 2
    d3_group = (5 / 10) ** 2 + (3 / 10) ** 2 + (2 / 10) ** 2
    expected = 2 + 9 / d2_group + 1 / d3_group + 5 / (1 / 4) + 3 / (1 / 6)

    assert _calc_enc(codons, aa_to_codons, degeneracy_family_counts) == pytest.approx(expected)


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
