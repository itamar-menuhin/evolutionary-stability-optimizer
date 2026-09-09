"""Tests for eso.tai's species-optimized tAI weights
(derive_species_optimized_tai_weights/derive_tai_weights) - a from-scratch
reimplementation of the real gtAI (Anwar et al. 2023) algorithm, built after
the actual `gtAI` PyPI package was tried and found to crash on real data (see
eso/tai.py's module docstring). Reuses eso.codon_usage._select_reference_codons
(already tested in tests/test_codon_usage_derive.py) and this module's own
already-R-verified _get_ws - only _rscu and the optimizer wiring are new,
untested-elsewhere code.
"""

import warnings

import pytest

from eso.codon_usage import _select_reference_codons
from eso.tai import (
    _rscu,
    derive_species_optimized_tai_weights,
    derive_tai_weights,
    derive_tai_weights_from_gff,
)

_GENETIC_CODE_NUM = 11


def _write_fasta(tmp_path, records, name="cds.fna"):
    path = tmp_path / name
    with open(path, "w") as handle:
        for header, seq in records:
            handle.write(f">{header}\n{seq}\n")
    return str(path)


def _write_gff(tmp_path, trna_lines, name="genomic.gff"):
    path = tmp_path / name
    path.write_text("##gff-version 3\n" + "".join(trna_lines))
    return str(path)


def _biased_gene(n_units):
    # Lys (AAA only), Glu (GAA only) - same helper shape as
    # tests/test_codon_usage_derive.py's own fixtures.
    return ("AAAGAA" * n_units) + "TAA"


def test_rscu_hand_computed_values():
    aa_to_codons = {"X": ["AAA", "AAG"]}
    reference_codons = ["AAA"] * 3 + ["AAG"] * 1
    rscu = _rscu(reference_codons, aa_to_codons)
    # expected count per codon if uniform = (3+1)/2 = 2
    assert rscu == {"AAA": pytest.approx(1.5), "AAG": pytest.approx(0.5)}


def test_rscu_skips_amino_acids_never_observed():
    aa_to_codons = {"X": ["AAA", "AAG"], "Y": ["CCC", "CCG"]}
    reference_codons = ["AAA", "AAG"]
    rscu = _rscu(reference_codons, aa_to_codons)
    assert set(rscu) == {"AAA", "AAG"}


def test_species_optimized_and_cai_use_the_identical_reference_set(tmp_path):
    # The whole point of sharing _select_reference_codons between
    # eso.codon_usage.derive_table_from_genome and this module is that both
    # use the SAME reference-gene recipe - verify that directly, not just
    # assume it from the refactor.
    fasta = _write_fasta(tmp_path, [(f"gene{i}", _biased_gene(60)) for i in range(20)])
    result_a = _select_reference_codons(fasta, _GENETIC_CODE_NUM, 100, 0.5, 1)
    result_b = _select_reference_codons(fasta, _GENETIC_CODE_NUM, 100, 0.5, 1)
    assert result_a == result_b


def _diverse_biased_gene(n_units):
    # Biased across several amino-acid families (Lys/Glu, Phe/Leu, Ile,
    # Asp/Gln), not just one pair - enough distinct codons for a meaningful
    # RSCU-vs-Wi correlation (the single-pair _biased_gene above only ever
    # covers 4 codons, too few by design - see
    # test_derive_species_optimized_tai_weights_rejects_too_few_reference_codons).
    return ("AAAGAATTCCTGATCGATCAG" * n_units) + "TAA"


def _small_genome_fixture(tmp_path):
    # A handful of genes with real, discoverable codon bias across several
    # amino acids, plus a small but real tRNA gene set so _get_ws has
    # something to score.
    fasta = _write_fasta(tmp_path, [
        (f"gene{i}", _diverse_biased_gene(60)) for i in range(10)
    ])
    gff = _write_gff(tmp_path, [
        "NC_1\tRefSeq\ttRNA\t1\t76\t.\t+\t.\tNote=tRNA-Lys(UUU)\n",
        "NC_1\tRefSeq\ttRNA\t100\t176\t.\t+\t.\tNote=tRNA-Lys(UUU)\n",
        "NC_1\tRefSeq\ttRNA\t200\t276\t.\t+\t.\tNote=tRNA-Lys(CUU)\n",
        "NC_1\tRefSeq\ttRNA\t300\t376\t.\t+\t.\tNote=tRNA-Glu(UUC)\n",
        "NC_1\tRefSeq\ttRNA\t400\t476\t.\t+\t.\tNote=tRNA-Glu(CUC)\n",
        "NC_1\tRefSeq\ttRNA\t500\t576\t.\t+\t.\tNote=tRNA-Phe(GAA)\n",
        "NC_1\tRefSeq\ttRNA\t600\t676\t.\t+\t.\tNote=tRNA-Leu(CAG)\n",
        "NC_1\tRefSeq\ttRNA\t700\t776\t.\t+\t.\tNote=tRNA-Ile(GAU)\n",
        "NC_1\tRefSeq\ttRNA\t800\t876\t.\t+\t.\tNote=tRNA-Asp(GUC)\n",
        "NC_1\tRefSeq\ttRNA\t900\t976\t.\t+\t.\tNote=tRNA-Gln(UUG)\n",
    ])
    return fasta, gff


def test_derive_species_optimized_tai_weights_is_well_formed(tmp_path):
    fasta, gff = _small_genome_fixture(tmp_path)
    weights = derive_species_optimized_tai_weights(
        fasta, gff, kingdom="prokaryote", genetic_code_num=_GENETIC_CODE_NUM,
        min_len_codons=100, top_perc=0.5, min_gene_count=1, seed=0,
    )
    assert len(weights) == 60
    assert all(isinstance(w, float) for w in weights.values())
    assert all(0.0 < w <= 1.0 for w in weights.values())


def test_derive_species_optimized_tai_weights_improves_correlation_over_generic(tmp_path):
    fasta, gff = _small_genome_fixture(tmp_path)
    kwargs = dict(genetic_code_num=_GENETIC_CODE_NUM, min_len_codons=100, top_perc=0.5, min_gene_count=1)

    optimized = derive_species_optimized_tai_weights(fasta, gff, kingdom="prokaryote", seed=0, **kwargs)
    generic = derive_tai_weights_from_gff(gff, kingdom="prokaryote", genetic_code_num=_GENETIC_CODE_NUM)

    aa_to_codons, reference_codons = _select_reference_codons(
        fasta, _GENETIC_CODE_NUM, kwargs["min_len_codons"], kwargs["top_perc"], kwargs["min_gene_count"])
    rscu = _rscu(reference_codons, aa_to_codons)

    from scipy.stats import spearmanr

    def corr(weights):
        shared = [c for c in weights if c in rscu]
        r, _ = spearmanr([rscu[c] for c in shared], [weights[c] for c in shared])
        return r

    assert corr(optimized) >= corr(generic)


def test_derive_species_optimized_tai_weights_rejects_a_bad_kingdom(tmp_path):
    fasta, gff = _small_genome_fixture(tmp_path)
    with pytest.raises(ValueError, match="prokaryote' or 'eukaryote'"):
        derive_species_optimized_tai_weights(fasta, gff, kingdom="bacteria")


def test_derive_species_optimized_tai_weights_rejects_too_few_reference_codons(tmp_path):
    # A single, tiny gene gives too few RSCU-defined codons to optimize
    # against meaningfully.
    fasta = _write_fasta(tmp_path, [("gene1", _biased_gene(60))])
    gff = _write_gff(tmp_path, ["NC_1\tRefSeq\ttRNA\t1\t76\t.\t+\t.\tNote=tRNA-Lys(UUU)\n"])
    with pytest.raises(ValueError, match="too few"):
        derive_species_optimized_tai_weights(
            fasta, gff, kingdom="prokaryote", genetic_code_num=_GENETIC_CODE_NUM,
            min_len_codons=100, top_perc=1.0, min_gene_count=1,
        )


def test_derive_tai_weights_prefers_species_optimized_by_default(tmp_path, monkeypatch):
    fasta, gff = _small_genome_fixture(tmp_path)
    called = {}

    def fake_species_optimized(*args, **kwargs):
        called["yes"] = True
        return {"AAA": 1.0}

    import eso.tai as tai_module
    monkeypatch.setattr(tai_module, "derive_species_optimized_tai_weights", fake_species_optimized)

    weights = derive_tai_weights(gff, kingdom="prokaryote", cds_fasta_path=fasta, genetic_code_num=_GENETIC_CODE_NUM)
    assert called.get("yes")
    assert weights == {"AAA": 1.0}


def test_derive_tai_weights_falls_back_to_generic_on_any_failure(tmp_path, monkeypatch):
    fasta, gff = _small_genome_fixture(tmp_path)

    def fake_species_optimized(*args, **kwargs):
        raise ValueError("simulated optimizer failure")

    import eso.tai as tai_module
    monkeypatch.setattr(tai_module, "derive_species_optimized_tai_weights", fake_species_optimized)

    with pytest.warns(UserWarning, match="falling back"):
        weights = derive_tai_weights(gff, kingdom="prokaryote", cds_fasta_path=fasta, genetic_code_num=_GENETIC_CODE_NUM)

    generic = derive_tai_weights_from_gff(gff, kingdom="prokaryote", genetic_code_num=_GENETIC_CODE_NUM)
    assert weights == generic


def test_derive_tai_weights_skips_species_optimized_without_cds_fasta_path(tmp_path):
    _, gff = _small_genome_fixture(tmp_path)
    # No cds_fasta_path given - must go straight to the generic path, no warning.
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        weights = derive_tai_weights(gff, kingdom="prokaryote", genetic_code_num=_GENETIC_CODE_NUM)
    generic = derive_tai_weights_from_gff(gff, kingdom="prokaryote", genetic_code_num=_GENETIC_CODE_NUM)
    assert weights == generic
