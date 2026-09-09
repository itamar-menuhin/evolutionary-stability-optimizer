"""Tests for eso.ncbi_genome.resolve_assembly_accession/fetch_genome_package_for -
resolving a bare species name/TaxID to its official NCBI assembly accession.
No real network access - urllib.request.urlopen is mocked throughout (see
tests/test_ncbi_genome_integration.py for the real, opt-in network version).
"""

import json
from unittest.mock import patch

import pytest

from eso.ncbi_genome import GenomeFetchError, fetch_genome_package_for, resolve_assembly_accession


class _FakeJsonResponse:
    def __init__(self, data):
        self._body = json.dumps(data).encode()

    def read(self, size=None):
        if size is None:
            body, self._body = self._body, b""
            return body
        chunk, self._body = self._body[:size], self._body[size:]
        return chunk

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        return False


def _report(accession, organism_name):
    return {"accession": accession, "organism": {"organism_name": organism_name}}


def test_single_report_resolves_directly():
    data = {"reports": [_report("GCF_000091665.1", "Methanocaldococcus jannaschii DSM 2661")]}
    with patch('urllib.request.urlopen', return_value=_FakeJsonResponse(data)):
        assert resolve_assembly_accession("Methanocaldococcus jannaschii") == "GCF_000091665.1"


def test_multi_report_same_species_picks_lowest_accession_and_warns():
    # E. coli-shaped fixture: 2 official assemblies, same species.
    data = {"reports": [
        _report("GCF_000008865.2", "Escherichia coli O157:H7 str. Sakai"),
        _report("GCF_000005845.2", "Escherichia coli str. K-12 substr. MG1655"),
    ]}
    with patch('urllib.request.urlopen', return_value=_FakeJsonResponse(data)):
        with pytest.warns(UserWarning, match="2 official assemblies"):
            accession = resolve_assembly_accession("Escherichia coli")
    assert accession == "GCF_000005845.2"  # the lower of the two


def test_multi_report_different_species_raises_clear_error():
    # Zobellia/Xanthomonas-shaped fixture: multiple reports, genuinely
    # different species (a genus-level query) - must not silently pick one.
    data = {"reports": [
        _report("GCF_000973105.1", "Zobellia galactanivorans"),
        _report("GCF_009725995.2", "Zobellia laminariae"),
        _report("GCF_029323795.1", "Zobellia alginiliquefaciens"),
    ]}
    with patch('urllib.request.urlopen', return_value=_FakeJsonResponse(data)):
        with pytest.raises(GenomeFetchError, match="3 different species"):
            resolve_assembly_accession("Zobellia")


def test_no_reports_raises_clear_error():
    with patch('urllib.request.urlopen', return_value=_FakeJsonResponse({})):
        with pytest.raises(GenomeFetchError, match="doesn't resolve to any NCBI organism"):
            resolve_assembly_accession("ThisIsNotARealSpecies123")


def test_network_failure_gives_friendly_message():
    import urllib.error

    with patch('urllib.request.urlopen', side_effect=urllib.error.URLError('boom')):
        with pytest.raises(GenomeFetchError, match="Could not look up"):
            resolve_assembly_accession("Escherichia coli")


def test_fetch_genome_package_for_skips_resolution_for_an_accession_shaped_input(monkeypatch):
    called = {"resolve": False, "fetch": None}

    def fake_resolve(organism):
        called["resolve"] = True
        return "should-not-be-used"

    def fake_fetch(accession, dest_dir=None):
        called["fetch"] = accession
        return "sentinel-package"

    import eso.ncbi_genome as ncbi_genome_module
    monkeypatch.setattr(ncbi_genome_module, "resolve_assembly_accession", fake_resolve)
    monkeypatch.setattr(ncbi_genome_module, "fetch_genome_package", fake_fetch)

    result = fetch_genome_package_for("GCF_000005845.2")

    assert called["resolve"] is False
    assert called["fetch"] == "GCF_000005845.2"
    assert result == "sentinel-package"


def test_fetch_genome_package_for_resolves_a_name_shaped_input(monkeypatch):
    called = {"resolve_arg": None, "fetch": None}

    def fake_resolve(organism):
        called["resolve_arg"] = organism
        return "GCF_000091665.1"

    def fake_fetch(accession, dest_dir=None):
        called["fetch"] = accession
        return "sentinel-package"

    import eso.ncbi_genome as ncbi_genome_module
    monkeypatch.setattr(ncbi_genome_module, "resolve_assembly_accession", fake_resolve)
    monkeypatch.setattr(ncbi_genome_module, "fetch_genome_package", fake_fetch)

    result = fetch_genome_package_for("Methanocaldococcus jannaschii")

    assert called["resolve_arg"] == "Methanocaldococcus jannaschii"
    assert called["fetch"] == "GCF_000091665.1"
    assert result == "sentinel-package"


def test_fetch_genome_package_for_recognizes_gca_accessions_too(monkeypatch):
    # GCA_ (GenBank), not just GCF_ (RefSeq), must also be treated as an
    # already-resolved accession, not a name to look up.
    called = {"resolve": False}

    def fake_resolve(organism):
        called["resolve"] = True
        return "should-not-be-used"

    def fake_fetch(accession, dest_dir=None):
        return accession

    import eso.ncbi_genome as ncbi_genome_module
    monkeypatch.setattr(ncbi_genome_module, "resolve_assembly_accession", fake_resolve)
    monkeypatch.setattr(ncbi_genome_module, "fetch_genome_package", fake_fetch)

    result = fetch_genome_package_for("GCA_000005845.2")
    assert called["resolve"] is False
    assert result == "GCA_000005845.2"
