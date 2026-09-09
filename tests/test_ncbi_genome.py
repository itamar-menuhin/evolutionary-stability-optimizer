"""Tests for eso.ncbi_genome.fetch_genome_package - the network-fetch layer
behind eso.codon_usage.derive_table_from_genome/eso.tai.derive_tai_weights_from_gff.
No real network access - urllib.request.urlopen is mocked throughout.
"""

import io
import zipfile
from unittest.mock import patch

import pytest

from eso.ncbi_genome import GenomeFetchError, fetch_genome_package


def _fake_zip_bytes(accession, include_cds=True, include_gff=True, include_genome=True):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, 'w') as archive:
        prefix = f'ncbi_dataset/data/{accession}/'
        if include_cds:
            archive.writestr(prefix + 'cds_from_genomic.fna', '>gene1\nATGAAATAA\n')
        if include_gff:
            archive.writestr(prefix + 'genomic.gff', '##gff-version 3\n')
        if include_genome:
            archive.writestr(prefix + f'{accession}_asm_genomic.fna', '>NC_1\nACGT\n')
    return buffer.getvalue()


class _FakeResponse:
    def __init__(self, data):
        self._data = data

    def read(self, size):
        chunk, self._data = self._data[:size], self._data[size:]
        return chunk

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        return False


def test_fetches_and_extracts_a_real_looking_package(tmp_path):
    accession = 'GCF_000005845.2'
    with patch('urllib.request.urlopen', return_value=_FakeResponse(_fake_zip_bytes(accession))):
        package = fetch_genome_package(accession, dest_dir=tmp_path)

    assert package.cds_fasta_path.endswith('cds_from_genomic.fna')
    assert package.gff_path.endswith('genomic.gff')
    assert package.genome_fasta_path.endswith('_genomic.fna')
    assert package.genome_fasta_path != package.cds_fasta_path
    with open(package.cds_fasta_path) as handle:
        assert '>gene1' in handle.read()


def test_network_failure_gives_friendly_message(tmp_path):
    import urllib.error

    with patch('urllib.request.urlopen', side_effect=urllib.error.URLError('boom')):
        with pytest.raises(GenomeFetchError, match='Could not fetch'):
            fetch_genome_package('GCF_000005845.2', dest_dir=tmp_path)


class _DroppedConnectionResponse(_FakeResponse):
    """Returns one real (partial) chunk successfully, then raises on the next
    read - simulating a connection dropping mid-download, a genuine failure
    mode for a real, multi-second, multi-MB transfer, distinct from a
    failure at connection-open time (already covered by URLError)."""

    def __init__(self, data, error):
        super().__init__(data)
        self._error = error
        self._reads = 0

    def read(self, size):
        self._reads += 1
        if self._reads > 1:
            raise self._error
        return super().read(size)


def test_connection_reset_mid_download_gives_friendly_message(tmp_path):
    # Regression test for a real gap: ConnectionResetError isn't a
    # urllib.error.URLError subclass, so a connection dropping mid-transfer
    # (not at connection-open time, which URLError alone does cover) used to
    # crash with a raw, unhandled exception instead of this function's own
    # documented GenomeFetchError. Confirmed directly before this fix.
    response = _DroppedConnectionResponse(_fake_zip_bytes('GCF_000005845.2'), ConnectionResetError('reset'))
    with patch('urllib.request.urlopen', return_value=response):
        with pytest.raises(GenomeFetchError, match='Could not fetch'):
            fetch_genome_package('GCF_000005845.2', dest_dir=tmp_path)


def test_incomplete_read_mid_download_gives_friendly_message(tmp_path):
    # Same real gap, via the OTHER exception http.client can raise for a
    # transfer that ends early (the server itself closing the connection
    # before the declared content length is reached) - also not a URLError
    # subclass.
    import http.client

    response = _DroppedConnectionResponse(
        _fake_zip_bytes('GCF_000005845.2'), http.client.IncompleteRead(b'partial'))
    with patch('urllib.request.urlopen', return_value=response):
        with pytest.raises(GenomeFetchError, match='Could not fetch'):
            fetch_genome_package('GCF_000005845.2', dest_dir=tmp_path)


def test_not_a_zip_file_gives_friendly_message(tmp_path):
    with patch('urllib.request.urlopen', return_value=_FakeResponse(b'not actually a zip file')):
        with pytest.raises(GenomeFetchError, match="isn't a valid genome-package zip"):
            fetch_genome_package('GCF_bad_accession.1', dest_dir=tmp_path)


def test_missing_cds_fasta_gives_friendly_message(tmp_path):
    accession = 'GCF_000005845.2'
    with patch('urllib.request.urlopen', return_value=_FakeResponse(_fake_zip_bytes(accession, include_cds=False))):
        with pytest.raises(GenomeFetchError, match='no CDS FASTA'):
            fetch_genome_package(accession, dest_dir=tmp_path)


def test_missing_gff_gives_friendly_message(tmp_path):
    accession = 'GCF_000005845.2'
    with patch('urllib.request.urlopen', return_value=_FakeResponse(_fake_zip_bytes(accession, include_gff=False))):
        with pytest.raises(GenomeFetchError, match='no GFF3 annotation'):
            fetch_genome_package(accession, dest_dir=tmp_path)


def test_missing_genome_fasta_gives_friendly_message(tmp_path):
    accession = 'GCF_000005845.2'
    with patch('urllib.request.urlopen', return_value=_FakeResponse(_fake_zip_bytes(accession, include_genome=False))):
        with pytest.raises(GenomeFetchError, match='no whole-genome FASTA'):
            fetch_genome_package(accession, dest_dir=tmp_path)
