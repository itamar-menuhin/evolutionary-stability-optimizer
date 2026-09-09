"""One real, opt-in, end-to-end network test against NCBI's live Datasets
API - deliberately NOT part of the routine test suite (no real network calls
anywhere else in this codebase's tests). Skipped unless the
ESO_RUN_NETWORK_TESTS environment variable is set, e.g.:

    ESO_RUN_NETWORK_TESTS=1 pytest tests/test_ncbi_genome_integration.py

Uses the smallest of the three organisms verified by hand this session
(M. jannaschii, GCF_000091665.1, ~1.2MB) - catches real-world regressions
(NCBI API/package-format changes) that a mocked unit test can't.
"""

import os

import pytest

from eso.codon_usage import derive_table_from_genome, detect_genetic_code_num_from_gff
from eso.ncbi_genome import fetch_genome_package
from eso.tai import derive_tai_weights_from_gff

pytestmark = pytest.mark.skipif(
    not os.environ.get("ESO_RUN_NETWORK_TESTS"),
    reason="real-network test - set ESO_RUN_NETWORK_TESTS=1 to run",
)

_ARCHAEON_ACCESSION = "GCF_000091665.1"  # Methanocaldococcus jannaschii DSM 2661


def test_fetch_derive_cai_and_tai_end_to_end_for_a_real_archaeon(tmp_path):
    package = fetch_genome_package(_ARCHAEON_ACCESSION, dest_dir=tmp_path)
    genetic_code_num = detect_genetic_code_num_from_gff(package.gff_path)
    assert genetic_code_num == 11

    cai_table = derive_table_from_genome(package.cds_fasta_path, genetic_code_num)
    assert 'K' in cai_table
    assert all(0.0 <= v <= 1.0 for freqs in cai_table.values() for v in freqs.values())

    tai_weights = derive_tai_weights_from_gff(
        package.gff_path, kingdom="prokaryote", genetic_code_num=genetic_code_num,
        genome_fasta_path=package.genome_fasta_path)
    assert len(tai_weights) == 60
    assert all(0.0 < w <= 1.0 for w in tai_weights.values())
