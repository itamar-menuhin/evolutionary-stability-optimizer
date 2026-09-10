"""Fetch a single organism's genome package (CDS FASTA + GFF3 annotation)
directly from NCBI's Datasets API, given its RefSeq/GenBank assembly
accession.

This is the network-I/O layer only - see eso.codon_usage.derive_table_from_genome
and eso.tai.derive_tai_weights_from_gff for what to do with the files this
returns. Download cost is negligible in practice (measured: 1.2-9MB,
1.5-8.7s for three real test organisms spanning bacteria/eukaryote/archaea),
so this is meant to be called directly, not cached/bundled ahead of time.

Mirrors the exact download approach already used and proven in
`prepare_ecoli_host_reference_data.py` (a sibling project's real,
already-working NCBI-fetch code) - stdlib-only (urllib + zipfile), no new
dependency.
"""

import http.client
import json
import re
import urllib.error
import urllib.parse
import urllib.request
import warnings
import zipfile
from dataclasses import dataclass
from pathlib import Path
from tempfile import mkdtemp

_CHUNK_SIZE = 1024 * 1024
_NCBI_DOWNLOAD_URL_TEMPLATE = (
    "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}"
    "/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,CDS_FASTA"
)
_NCBI_TAXON_REPORT_URL_TEMPLATE = (
    "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{name_or_taxid}"
    "/dataset_report?filters.reference_only=true&filters.assembly_source=RefSeq"
)
#: A real assembly accession, e.g. "GCF_000005845.2" or "GCA_000005845". If a
#: string doesn't match this, fetch_genome_package_for treats it as a species
#: name/TaxID to resolve instead.
_ACCESSION_RE = re.compile(r'^GC[AF]_\d+(\.\d+)?$')


class GenomeFetchError(Exception):
    """Fetching or unpacking an NCBI genome package failed.

    Raised with a plain-English message - the most common real causes are a
    typo'd/withdrawn assembly accession or a network problem, neither of
    which should require reading a traceback through this module to
    diagnose.
    """


@dataclass(frozen=True)
class GenomePackage:
    """Paths to the files this module's callers actually need, inside the
    extracted NCBI genome package. `genome_fasta_path` (the whole-genome
    sequence, not just CDS regions) is needed by eso.tai.derive_tai_weights_from_gff
    for tRNA annotations that only give an anticodon's genomic position
    (`anticodon=(pos:...)`, common in tRNAscan-SE-sourced annotations) rather
    than embedding the anticodon triplet directly in the feature's own
    attributes (as RefSeq's own curated E. coli annotation does) - confirmed
    this session that both styles occur in real NCBI packages, not just one."""

    cds_fasta_path: str
    gff_path: str
    genome_fasta_path: str


#: Genuine network-layer failures - a connection dropping mid-transfer
#: (ConnectionError, an OSError subclass but NOT a urllib.error.URLError
#: one) or the server closing early (http.client.IncompleteRead, its own
#: separate exception hierarchy) - in addition to urllib's own URLError
#: (raised at connection-open time: DNS failure, refused connection, a
#: non-2xx status). Deliberately narrower than bare OSError: a previous
#: version of this fix caught OSError broadly around the whole download,
#: which also happened to catch a LOCAL disk error (permission denied,
#: disk full) from opening/writing the destination file and misreported it
#: as a network problem - confirmed directly this was a real, reachable
#: regression, not hypothetical. Local I/O errors are now caught separately,
#: right where they occur, with their own honest message.
_NETWORK_ERRORS = (urllib.error.URLError, ConnectionError, TimeoutError, http.client.HTTPException)


def _download(url, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")

    try:
        try:
            handle = open(tmp, "wb")
        except OSError as exc:
            raise GenomeFetchError(
                f"Could not create '{tmp}' to download into - check that the destination "
                f"directory is writable and has enough free disk space. The underlying error "
                f"was: {exc!r}."
            ) from exc

        with handle:
            try:
                with urllib.request.urlopen(url, timeout=120) as response:
                    while True:
                        chunk = response.read(_CHUNK_SIZE)
                        if not chunk:
                            break
                        try:
                            handle.write(chunk)
                        except OSError as exc:
                            raise GenomeFetchError(
                                f"Could not write to '{tmp}' - check that the destination disk "
                                f"has enough free space. The underlying error was: {exc!r}."
                            ) from exc
            except _NETWORK_ERRORS as exc:
                raise GenomeFetchError(
                    f"Could not fetch the NCBI genome package from {url!r}. This is usually either "
                    f"a network problem, or '{url.rsplit('/accession/', 1)[-1].split('/download')[0]}' "
                    f"isn't a real/current NCBI assembly accession - check it at "
                    f"https://www.ncbi.nlm.nih.gov/datasets/genome. The underlying error was: {exc!r}."
                ) from exc
    except GenomeFetchError:
        tmp.unlink(missing_ok=True)
        raise

    tmp.replace(dest)


def fetch_genome_package(assembly_accession, dest_dir=None):
    """Download and extract one organism's genome package from NCBI.

    `assembly_accession` is an explicit RefSeq/GenBank accession, e.g.
    `"GCF_000005845.2"` (E. coli K-12 MG1655) - find one for your organism
    at https://www.ncbi.nlm.nih.gov/datasets/genome, or reuse one you
    already know. This does NOT resolve a bare species name or TaxID to its
    official assembly on its own - use `fetch_genome_package_for` for that.

    `dest_dir`: where to download/extract into (default: a fresh temp
    directory). The returned paths live under here - if you pass your own
    `dest_dir`, you're responsible for cleaning it up afterward.

    Returns
    -------
    GenomePackage

    Raises
    ------
    GenomeFetchError
        On a bad accession, a network failure, or a package that doesn't
        contain what was requested (e.g. an assembly with no CDS annotation).
    """
    work_dir = Path(dest_dir) if dest_dir is not None else Path(mkdtemp(prefix="eso_ncbi_genome_"))
    zip_path = work_dir / f"{assembly_accession}.zip"
    url = _NCBI_DOWNLOAD_URL_TEMPLATE.format(accession=assembly_accession)

    _download(url, zip_path)

    extract_dir = work_dir / "extracted"
    try:
        with zipfile.ZipFile(zip_path) as archive:
            archive.extractall(extract_dir)
    except zipfile.BadZipFile as exc:
        raise GenomeFetchError(
            f"NCBI returned something that isn't a valid genome-package zip file for "
            f"'{assembly_accession}' - it's likely not a real/current assembly accession. "
            f"Check it at https://www.ncbi.nlm.nih.gov/datasets/genome."
        ) from exc

    data_dir = extract_dir / "ncbi_dataset" / "data" / assembly_accession
    if not data_dir.is_dir():
        raise GenomeFetchError(
            f"The downloaded package for '{assembly_accession}' didn't have the expected "
            f"folder layout (looked for {data_dir}) - NCBI's package format may have changed, "
            f"or this accession doesn't exist."
        )

    cds_candidates = list(data_dir.glob("cds_from_genomic.fna"))
    gff_candidates = list(data_dir.glob("genomic.gff"))
    # Matched by accession prefix, not just "*_genomic.fna" - that pattern
    # would also match "cds_from_genomic.fna" itself.
    genome_candidates = list(data_dir.glob(f"{assembly_accession}_*_genomic.fna"))
    if not cds_candidates:
        raise GenomeFetchError(
            f"'{assembly_accession}' has no CDS FASTA in its NCBI package - this assembly "
            f"likely has no protein-coding annotation available, so a codon-usage/tAI table "
            f"can't be derived from it. Try a different (typically RefSeq) accession for the "
            f"same organism at https://www.ncbi.nlm.nih.gov/datasets/genome."
        )
    if not gff_candidates:
        raise GenomeFetchError(
            f"'{assembly_accession}' has no GFF3 annotation in its NCBI package."
        )
    if not genome_candidates:
        raise GenomeFetchError(
            f"'{assembly_accession}' has no whole-genome FASTA in its NCBI package."
        )

    return GenomePackage(
        cds_fasta_path=str(cds_candidates[0]),
        gff_path=str(gff_candidates[0]),
        genome_fasta_path=str(genome_candidates[0]),
    )


def _species_key(organism_name):
    """First two whitespace-separated tokens (genus + species epithet) of an
    NCBI organism_name - used to tell "the same species, multiple official
    strains" (e.g. two E. coli strains) apart from "genuinely different
    species" (e.g. a genus-level query matching many species), since
    individual reports' own tax_id is strain-specific in both cases and
    can't be used for this directly. Verified against ~40 real organisms
    this session, including both a real multi-strain-one-species case
    (E. coli) and two independent real multi-species cases (genus-level
    queries "Zobellia", "Xanthomonas")."""
    return tuple(organism_name.split()[:2])


def resolve_assembly_accession(species_name_or_taxid):
    """Resolve a bare species name (e.g. `"Escherichia coli"`) or NCBI TaxID
    (e.g. `562`, as an int or numeric string) to its official RefSeq
    assembly accession, via NCBI's Datasets API - no bulk-file download
    needed. Works identically for a species name or a TaxID (verified
    directly against ~40 real organisms this session).

    If more than one assembly is officially designated (`"reference
    genome"` in NCBI's own terms - rare in practice, but real: e.g. E. coli
    has two, for K-12 MG1655 and O157:H7 Sakai) and they're all the *same*
    species, the lowest accession is picked and a warning names the others -
    pass that specific assembly's own accession to `fetch_genome_package`
    directly instead of using resolution, if you need a particular one.

    If the name/TaxID is too broad and matches multiple *different* species
    (e.g. a bare genus name), or doesn't resolve to any NCBI organism at
    all, or resolves but has no officially-designated assembly, raises
    `GenomeFetchError` explaining which case it was and how to proceed (a
    more specific name, e.g. a full strain name, or an explicit accession).

    Returns
    -------
    str - an assembly accession, suitable for `fetch_genome_package`.
    """
    url = _NCBI_TAXON_REPORT_URL_TEMPLATE.format(name_or_taxid=urllib.parse.quote(str(species_name_or_taxid)))
    try:
        # See _download's identical fix and _NETWORK_ERRORS' own docstring:
        # URLError/TimeoutError alone misses a connection dropping
        # mid-response. No local file I/O happens in this function, so
        # there's no risk of misattributing a local error here the way a
        # bare `except OSError` would for _download - still using the same
        # narrower _NETWORK_ERRORS tuple for consistency, not because it's
        # required here.
        with urllib.request.urlopen(url, timeout=30) as response:
            raw = response.read()
    except _NETWORK_ERRORS as exc:
        raise GenomeFetchError(
            f"Could not look up '{species_name_or_taxid}' on NCBI - this is usually a network "
            f"problem. The underlying error was: {exc!r}."
        ) from exc

    try:
        data = json.loads(raw)
    except json.JSONDecodeError as exc:
        # Confirmed real: NCBI (or an intermediary proxy) returning an error
        # page - maintenance, rate-limiting, an outage - as HTML/plain text
        # instead of JSON isn't a network-level failure urlopen itself would
        # catch, but json.loads on it raises unhelpfully with no indication
        # this was an NCBI-side response problem rather than a malformed
        # species name.
        raise GenomeFetchError(
            f"NCBI returned something that isn't valid JSON when looking up "
            f"'{species_name_or_taxid}' - this usually means NCBI itself is having a problem "
            f"(maintenance, rate-limiting, an outage), not that the name/TaxID is wrong. Try "
            f"again in a moment, or check https://www.ncbi.nlm.nih.gov/datasets/genome directly."
        ) from exc

    reports = data.get('reports', [])
    if not reports:
        raise GenomeFetchError(
            f"'{species_name_or_taxid}' doesn't resolve to any NCBI organism with an official "
            f"(\"reference genome\") RefSeq assembly - either it's not a real/complete species "
            f"name or TaxID (NCBI's own name matching is exact-ish; check spelling, or search "
            f"at https://www.ncbi.nlm.nih.gov/datasets/genome), or this organism genuinely has "
            f"no NCBI-designated official assembly. Either way, find an accession yourself at "
            f"that same URL and pass it to fetch_genome_package directly instead."
        )

    species_keys = {_species_key(r['organism']['organism_name']) for r in reports}
    if len(species_keys) > 1:
        distinct_species = sorted({r['organism']['organism_name'] for r in reports})
        raise GenomeFetchError(
            f"'{species_name_or_taxid}' matches {len(distinct_species)} different species with "
            f"an official assembly, not one: {', '.join(distinct_species[:8])}"
            f"{', ...' if len(distinct_species) > 8 else ''}. This is usually because the name "
            f"given is a genus (or other higher taxonomic rank), not a specific species. Pass a "
            f"full species or strain name (e.g. one of the ones listed above) instead, or an "
            f"explicit accession."
        )

    reports_sorted = sorted(reports, key=lambda r: r['accession'])
    chosen = reports_sorted[0]
    if len(reports_sorted) > 1:
        other_accessions = [r['accession'] for r in reports_sorted[1:]]
        warnings.warn(
            f"'{species_name_or_taxid}' has {len(reports_sorted)} official assemblies "
            f"({', '.join(r['accession'] for r in reports_sorted)}) - using {chosen['accession']!r} "
            f"(the lowest accession). Pass a specific one of {other_accessions} to "
            f"fetch_genome_package directly if you want a different one.",
            stacklevel=2,
        )
    return chosen['accession']


def fetch_genome_package_for(organism, dest_dir=None):
    """Recommended one-stop entry point: fetch an organism's genome package
    given either an explicit assembly accession (used directly, exactly
    like `fetch_genome_package`) or a bare species name/TaxID (resolved to
    its official accession first, via `resolve_assembly_accession`).

    Returns
    -------
    GenomePackage

    Raises
    ------
    GenomeFetchError
        See `resolve_assembly_accession` (for a name/TaxID) or
        `fetch_genome_package` (for the actual download) for the specific
        failure cases.
    """
    accession = organism if _ACCESSION_RE.match(str(organism)) else resolve_assembly_accession(organism)
    return fetch_genome_package(accession, dest_dir=dest_dir)
