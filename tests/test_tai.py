"""Tests for eso.tai - the tRNA Adaptation Index (tAI) implementation ported
from a sibling project's already real-R-verified `compute_ecoli_tai_weights.py`
(github.com/mariodosreis/tai's R/tAI.R get.ws(), dos Reis et al. 2004).
"""

import math

import pytest

from eso.tai import (
    CODON_ORDER,
    _get_ws,
    _parse_trna_gene_counts,
    build_tai_score_fn,
    derive_tai_weights_from_gff,
)

# Real tRNA gene copy number vector, indexed by CODON_ORDER, counted from
# E. coli K-12 MG1655's real NCBI RefSeq annotation (GCF_000005845.2,
# genomic.gff) - 86 tRNA genes total, the same exact count and vector this
# session's sibling project (STABLES) independently verified against a real
# `Rscript` run of the reference implementation.
_REAL_ECOLI_TRNA_COUNTS = [
    0, 2, 1, 1, 0, 2, 1, 1, 0, 3, 0, 0, 0, 1, 1, 1, 0, 1, 1, 4, 0, 1, 1, 1, 0, 1, 2, 2, 4, 0, 0, 1,
    0, 3, 0, 8, 0, 2, 1, 1, 0, 4, 6, 0, 0, 1, 1, 1, 0, 2, 5, 0, 0, 2, 3, 0, 0, 3, 4, 0, 0, 4, 1, 1,
]

# Golden reference, reused verbatim from STABLES' own
# qa/tests/test_compute_ecoli_tai_weights.py - the real output of get.ws()
# run through an actual Rscript invocation of the unmodified R source
# (github.com/mariodosreis/tai), not re-derived from this port.
_GOLDEN_R_OUTPUT = {
    'TTT': 0.196666666666667, 'TTC': 0.333333333333333, 'TTA': 0.166666666666667, 'TTG': 0.22,
    'TCT': 0.196666666666667, 'TCC': 0.333333333333333, 'TCA': 0.166666666666667, 'TCG': 0.22,
    'TAT': 0.295, 'TAC': 0.5,
    'TGT': 0.0983333333333333, 'TGC': 0.166666666666667, 'TGG': 0.22,
    'CTT': 0.0983333333333333, 'CTC': 0.166666666666667, 'CTA': 0.166666666666667, 'CTG': 0.72,
    'CCT': 0.0983333333333333, 'CCC': 0.166666666666667, 'CCA': 0.166666666666667, 'CCG': 0.22,
    'CAT': 0.0983333333333333, 'CAC': 0.166666666666667, 'CAA': 0.333333333333333, 'CAG': 0.44,
    'CGT': 0.666666666666667, 'CGC': 0.48, 'CGA': 6.66666666666593e-05, 'CGG': 0.166666666666667,
    'ATT': 0.295, 'ATC': 0.5, 'ATA': 0.0183333333333333,
    'ACT': 0.196666666666667, 'ACC': 0.333333333333333, 'ACA': 0.166666666666667, 'ACG': 0.22,
    'AAT': 0.393333333333333, 'AAC': 0.666666666666667, 'AAA': 1.0, 'AAG': 0.32,
    'AGT': 0.0983333333333333, 'AGC': 0.166666666666667, 'AGA': 0.166666666666667, 'AGG': 0.22,
    'GTT': 0.196666666666667, 'GTC': 0.333333333333333, 'GTA': 0.833333333333333, 'GTG': 0.266666666666667,
    'GCT': 0.196666666666667, 'GCC': 0.333333333333333, 'GCA': 0.5, 'GCG': 0.16,
    'GAT': 0.295, 'GAC': 0.5, 'GAA': 0.666666666666667, 'GAG': 0.213333333333333,
    'GGT': 0.393333333333333, 'GGC': 0.666666666666667, 'GGA': 0.166666666666667, 'GGG': 0.22,
}


def test_codon_order_matches_get_ws_own_special_case_indices():
    assert CODON_ORDER[10] == 'TAA'
    assert CODON_ORDER[11] == 'TAG'
    assert CODON_ORDER[14] == 'TGA'
    assert CODON_ORDER[35] == 'ATG'


def test_get_ws_matches_real_r_reference_implementation():
    weights = _get_ws(_REAL_ECOLI_TRNA_COUNTS, sking=1)
    assert set(weights) == set(_GOLDEN_R_OUTPUT)
    for codon, golden_value in _GOLDEN_R_OUTPUT.items():
        assert abs(weights[codon] - golden_value) < 1e-9, f'{codon}: {weights[codon]} != {golden_value}'


def _write_gff(tmp_path, trna_lines):
    path = tmp_path / "genomic.gff"
    header = "##gff-version 3\n"
    path.write_text(header + "".join(trna_lines))
    return str(path)


def test_parse_trna_gene_counts_reads_note_style_anticodons(tmp_path):
    gff_path = _write_gff(tmp_path, [
        "NC_1\tRefSeq\ttRNA\t1\t76\t.\t+\t.\tID=trna-1;Note=tRNA-Lys(UUU)\n",
        "NC_1\tRefSeq\ttRNA\t100\t176\t.\t+\t.\tID=trna-2;Note=tRNA-Lys(UUU)\n",
    ])
    counts = _parse_trna_gene_counts(gff_path)
    # real NCBI annotations give the anticodon in RNA notation (U, not T) -
    # UUU reverse-complements (as DNA) to codon AAA.
    assert counts[CODON_ORDER.index('AAA')] == 2
    assert sum(counts) == 2


def test_parse_trna_gene_counts_reads_product_style_anticodons(tmp_path):
    gff_path = _write_gff(tmp_path, [
        "NC_1\tRefSeq\ttRNA\t1\t76\t.\t+\t.\tID=trna-1;product=tRNA-Glu(UUC)\n",
    ])
    counts = _parse_trna_gene_counts(gff_path)
    assert counts[CODON_ORDER.index('GAA')] == 1
    assert sum(counts) == 1


def test_parse_trna_gene_counts_resolves_a_position_only_anticodon(tmp_path):
    # Regression test: tRNAscan-SE-sourced annotations (confirmed this
    # session on a real archaeon, GCF_000091665.1) give the anticodon only
    # as a genomic position, not embedded in the feature's own attributes -
    # the original text-only parsing silently found zero tRNA genes for
    # every one of these, which would have shipped a completely broken
    # (empty) tAI weight table for any such organism.
    #
    # NC_1 sequence, 1-indexed: position 10-12 = "TTT" (forward strand).
    genome_seqs = {"NC_1": "AAAAAAAAATTTAAAAAAAAAA"}
    gff_path = _write_gff(tmp_path, [
        "NC_1\tRefSeq\ttRNA\t10\t12\t.\t+\t.\tID=trna-1;anticodon=(pos:10..12);product=tRNA-Lys\n",
    ])
    counts = _parse_trna_gene_counts(gff_path, genome_seqs=genome_seqs)
    # anticodon TTT (read directly off the genome, already DNA notation)
    # reverse-complements to codon AAA.
    assert counts[CODON_ORDER.index('AAA')] == 1
    assert sum(counts) == 1


def test_parse_trna_gene_counts_resolves_a_complement_strand_position(tmp_path):
    # Same real annotation style, on the complement strand -
    # `anticodon=(pos:complement(...))`.
    genome_seqs = {"NC_1": "AAAAAAAAAAAAGAAAAAAAAA"}  # positions 13-15 (1-based) = "GAA"
    gff_path = _write_gff(tmp_path, [
        "NC_1\tRefSeq\ttRNA\t13\t15\t.\t-\t.\tID=trna-1;anticodon=(pos:complement(13..15));product=tRNA-Ser\n",
    ])
    counts = _parse_trna_gene_counts(gff_path, genome_seqs=genome_seqs)
    # Two revcomps applied back to back (once for the complement strand,
    # once to go from anticodon to codon) cancel out - for a
    # complement-strand tRNA gene, the resulting codon equals the raw
    # forward-strand-read sequence itself. Verified against 35/36 real
    # tRNA features on the real archaeon assembly this was found on
    # (GCF_000091665.1), cross-checked against each gene's own annotated
    # amino acid.
    assert counts[CODON_ORDER.index('GAA')] == 1
    assert sum(counts) == 1


def test_parse_trna_gene_counts_raises_when_position_only_and_no_genome_given(tmp_path):
    gff_path = _write_gff(tmp_path, [
        "NC_1\tRefSeq\ttRNA\t10\t12\t.\t+\t.\tID=trna-1;anticodon=(pos:10..12);product=tRNA-Lys\n",
    ])
    with pytest.raises(ValueError, match="genomic position"):
        _parse_trna_gene_counts(gff_path, genome_seqs=None)


def test_derive_tai_weights_from_gff_requires_a_valid_kingdom(tmp_path):
    gff_path = _write_gff(tmp_path, ["NC_1\tRefSeq\ttRNA\t1\t76\t.\t+\t.\tNote=tRNA-Lys(UUU)\n"])
    with pytest.raises(ValueError, match="prokaryote' or 'eukaryote'"):
        derive_tai_weights_from_gff(gff_path, kingdom="bacteria")


def test_derive_tai_weights_from_gff_rejects_no_trna_genes(tmp_path):
    gff_path = _write_gff(tmp_path, ["NC_1\tRefSeq\tgene\t1\t76\t.\t+\t.\tID=gene-1\n"])
    with pytest.raises(ValueError, match="No tRNA gene features"):
        derive_tai_weights_from_gff(gff_path, kingdom="prokaryote")


def test_derive_tai_weights_from_gff_rejects_a_mismatched_genetic_code(tmp_path):
    gff_path = _write_gff(tmp_path, ["NC_1\tRefSeq\ttRNA\t1\t76\t.\t+\t.\tNote=tRNA-Lys(UUU)\n"])
    # Genetic code 4 (Mold/Protozoan Mitochondrial) reassigns TGA from a
    # stop codon to Trp - this tAI implementation assumes the standard 3
    # stop codons and should refuse rather than silently mis-score.
    with pytest.raises(ValueError, match="isn't supported"):
        derive_tai_weights_from_gff(gff_path, kingdom="prokaryote", genetic_code_num=4)


def test_build_tai_score_fn_is_the_geometric_mean_of_per_codon_weights():
    weights = {"AAA": 1.0, "GGG": 0.25}
    score_fn = build_tai_score_fn(weights)
    # AAA GGG AAA GGG: geometric mean of [1.0, 0.25, 1.0, 0.25] = sqrt(0.25) = 0.5
    assert math.isclose(score_fn("AAAGGGAAAGGG"), 0.5, rel_tol=1e-9)


def test_build_tai_score_fn_skips_codons_with_no_weight():
    # Stop codons/Met carry no weight (see derive_tai_weights_from_gff's
    # docstring) - a real ORF's trailing stop codon must not crash or
    # silently zero out the score.
    weights = {"AAA": 1.0}
    score_fn = build_tai_score_fn(weights)
    assert math.isclose(score_fn("AAATAA"), 1.0, rel_tol=1e-9)


def test_build_tai_score_fn_on_all_unweighted_codons_returns_zero():
    score_fn = build_tai_score_fn({"AAA": 1.0})
    assert score_fn("TAATAG") == 0.0
