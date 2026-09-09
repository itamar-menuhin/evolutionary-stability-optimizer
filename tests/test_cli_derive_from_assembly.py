"""Tests that --derive-codon-usage-table-from-assembly and
--derive-tai-score-from-assembly are wired correctly through eso.cli.main -
mocking eso.ncbi_genome.fetch_genome_package_for (no real network access), same
pattern as tests/test_cli_codon_usage_table.py.
"""

import eso.cli as cli
from eso.ncbi_genome import GenomeFetchError, GenomePackage


def _write(tmp_path, name, content):
    file_path = tmp_path / name
    file_path.write_text(content)
    return str(file_path)


def _fake_ecoli_package(tmp_path):
    gff_path = _write(
        tmp_path, "genomic.gff",
        "##gff-version 3\n"
        "NC_1\tRefSeq\tCDS\t1\t100\t.\t+\t0\tID=cds-1;transl_table=11\n"
        "NC_1\tRefSeq\ttRNA\t200\t276\t.\t+\t.\tID=trna-1;Note=tRNA-Lys(UUU)\n",
    )
    cds_path = _write(
        tmp_path, "cds_from_genomic.fna",
        ">gene1\n" + ("AAAGAA" * 40) + "TAA\n" + ">gene2\n" + ("AAAAAGGAAGAG" * 40) + "TAA\n",
    )
    genome_path = _write(tmp_path, "genome_genomic.fna", ">NC_1\n" + "ACGT" * 100 + "\n")
    return GenomePackage(cds_fasta_path=cds_path, gff_path=gff_path, genome_fasta_path=genome_path)


def _mock_run_pipeline(monkeypatch):
    captured = {}

    def fake_run_pipeline(**kwargs):
        captured.update(kwargs)
        return 'Success!', []

    monkeypatch.setattr(cli, 'run_pipeline', fake_run_pipeline)
    return captured


def test_derive_codon_usage_table_from_assembly_reaches_pipeline_main(tmp_path, monkeypatch):
    captured = _mock_run_pipeline(monkeypatch)
    package = _fake_ecoli_package(tmp_path)
    monkeypatch.setattr(cli, 'fetch_genome_package_for', lambda accession: package)

    exit_code = cli.main(['--derive-codon-usage-table-from-assembly', 'GCF_000005845.2'])

    assert exit_code == 0
    assert captured['codon_usage_table'] is not None
    assert 'K' in captured['codon_usage_table']


def test_derive_tai_score_from_assembly_reaches_pipeline_main(tmp_path, monkeypatch):
    captured = _mock_run_pipeline(monkeypatch)
    package = _fake_ecoli_package(tmp_path)
    monkeypatch.setattr(cli, 'fetch_genome_package_for', lambda accession: package)

    exit_code = cli.main([
        '--derive-tai-score-from-assembly', 'GCF_000005845.2',
        '--tai-kingdom', 'prokaryote',
    ])

    assert exit_code == 0
    assert callable(captured['custom_score_fn'])
    assert captured['custom_score_fn']('AAATAA') > 0


def test_derive_tai_score_from_assembly_without_kingdom_fails_fast(capsys):
    exit_code = cli.main(['--derive-tai-score-from-assembly', 'GCF_000005845.2'])

    assert exit_code == 1
    assert "--tai-kingdom" in capsys.readouterr().err


def test_derive_codon_usage_table_conflicts_with_codon_usage_table_file(tmp_path, capsys):
    file_path = _write(tmp_path, "good.csv", 'codon,aa,freq_within_aa\nGCT,A,0.5\nGCC,A,0.5\n')

    exit_code = cli.main([
        '--codon-usage-table-file', file_path,
        '--derive-codon-usage-table-from-assembly', 'GCF_000005845.2',
    ])

    assert exit_code == 1
    assert "can't both be given" in capsys.readouterr().err


def test_a_fetch_failure_fails_fast_with_the_underlying_message(monkeypatch, capsys):
    def fake_fetch(accession):
        raise GenomeFetchError(f"Could not fetch NCBI genome package for {accession!r}.")

    monkeypatch.setattr(cli, 'fetch_genome_package_for', fake_fetch)

    exit_code = cli.main(['--derive-codon-usage-table-from-assembly', 'GCF_bad.1'])

    assert exit_code == 1
    assert "Could not fetch" in capsys.readouterr().err
