"""tRNA Adaptation Index (tAI) - an alternative to eso.codon_usage's CAI-style
tables, scoring codon choice by genomic tRNA gene copy number instead of by
reference-gene codon frequency.

`derive_tai_weights_from_gff` is a near-verbatim port of the real, R-verified
dos Reis, Wernisch & Savva (2004) method already implemented and checked
against the method authors' own R package
(github.com/mariodosreis/tai/R/tAI.R -- get.ws()) in a sibling project's
`compute_ecoli_tai_weights.py`. Ported here rather than reimplemented from
the paper, to keep that verification meaningful.

`derive_species_optimized_tai_weights`/`derive_tai_weights` add a second,
optional way to get tAI weights: optimizing the same `get.ws()` formula's 5
free wobble parameters per-organism (matching what stAIcalc/gtAI do), rather
than using dos Reis' own fixed generic defaults. This reimplements the real
gtAI (Anwar et al. 2023) algorithm directly - NOT a wrapper around the actual
`gtAI` PyPI package, which was tried and abandoned after it crashed with a
real KeyError on real archaeon data (its own wobble-key computation produces
an invalid key for non-4-fold-degenerate codon families and looks it up with
no default). Reimplemented here using pieces already verified elsewhere in
this codebase: `eso.codon_usage._select_reference_codons` for the
highly-expressed-gene reference set (identical recipe gtAI's own `ENc_calc`
uses), and `_get_ws` itself (this module, already correct for every
degeneracy class, unlike gtAI's `abs_Wi`) for scoring any candidate wobble
parameter vector.
"""

import math
import re
import warnings
from collections import Counter

import numpy as np
from Bio import SeqIO
from Bio.Data.CodonTable import unambiguous_dna_by_id
from scipy.optimize import differential_evolution
from scipy.stats import spearmanr

_BASES = ('T', 'C', 'A', 'G')
#: Standard genetic-code-table codon order: first base slowest, third base
#: fastest, in T-C-A-G order at every position - the positional order
#: get.ws() implicitly assumes for its 64-length input vector (verified
#: against the real R source in compute_ecoli_tai_weights.py: this order
#: puts TAA at index 10, TAG at 11, TGA at 14, ATG at 35, exactly matching
#: get.ws()'s own special-cased indices for stop codons/methionine).
CODON_ORDER = tuple(b1 + b2 + b3 for b1 in _BASES for b2 in _BASES for b3 in _BASES)
assert CODON_ORDER[10] == 'TAA' and CODON_ORDER[11] == 'TAG' and CODON_ORDER[14] == 'TGA' and CODON_ORDER[35] == 'ATG'

#: Optimised s-values from dos Reis et al. (2004), exactly as hardcoded in
#: get.ws()'s own default.
_S_VALUES = (0.0, 0.0, 0.0, 0.0, 0.41, 0.28, 0.9999, 0.68, 0.89)

_COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'U': 'A', 'T': 'A'}

_KINGDOM_TO_SKING = {'prokaryote': 1, 'eukaryote': 0}

_NOTE_ANTICODON_RE = re.compile(r'Note=tRNA-[\w ]+\(([ACGU]{3})\)')
_PRODUCT_ANTICODON_RE = re.compile(r'product=tRNA-\w+\(([ACGU]{3})\)')
#: Fallback for annotations that give the anticodon's genomic position
#: instead of embedding the triplet directly (tRNAscan-SE-sourced GFFs -
#: confirmed this session on a real archaeon assembly, GCF_000091665.1,
#: where every tRNA feature looks like `anticodon=(pos:97665..97667)` or
#: `anticodon=(pos:complement(97501..97503))`, with no anticodon text
#: anywhere in the line - unlike RefSeq's own curated E. coli annotation,
#: which always embeds it directly. Needs the genome sequence itself to
#: resolve (see `derive_tai_weights_from_gff`'s `genome_fasta_path`).
_ANTICODON_POS_COMPLEMENT_RE = re.compile(r'anticodon=\(pos:complement\((\d+)\.\.(\d+)\)\)')
_ANTICODON_POS_RE = re.compile(r'anticodon=\(pos:(\d+)\.\.(\d+)\)')


def _revcomp(seq):
    return ''.join(_COMPLEMENT[base] for base in reversed(seq))


def _get_ws(trna, sking, s=_S_VALUES):
    """Faithful port of get.ws() from github.com/mariodosreis/tai's R/tAI.R.
    `trna`: a 64-length vector indexed by CODON_ORDER. `sking`: 0=Eukaryota,
    1=Prokaryota. Returns a dict of {codon: weight} for the 60 codons that
    remain after dropping the 3 stop codons and methionine (get.ws()'s own
    scope - Trp/Met have only one codon each, nothing to weight)."""
    if len(trna) != 64:
        raise ValueError(f'trna must have length 64, got {len(trna)}')

    p = [1 - value for value in s]
    W = [0.0] * 64

    for i in range(0, 61, 4):
        W[i] = p[0] * trna[i] + p[4] * trna[i + 1]
        W[i + 1] = p[1] * trna[i + 1] + p[5] * trna[i]
        W[i + 2] = p[2] * trna[i + 2] + p[6] * trna[i]
        W[i + 3] = p[3] * trna[i + 3] + p[7] * trna[i + 2]

    W[35] = p[3] * trna[35]  # Methionine (0-based index 35) - special case, not the block formula

    if sking == 1:
        W[34] = p[8]  # Bacteria-specific: ATA (Ile) wobble, a fixed constant

    remove_idx = {10, 11, 14, 35}
    kept_codons = [c for idx, c in enumerate(CODON_ORDER) if idx not in remove_idx]
    kept_weights = [w for idx, w in enumerate(W) if idx not in remove_idx]

    max_w = max(kept_weights)
    w_norm = [w / max_w for w in kept_weights]

    zero_free = [w for w in w_norm if w != 0]
    if len(zero_free) < len(w_norm):
        geometric_mean = math.exp(sum(math.log(w) for w in zero_free) / len(zero_free))
        w_norm = [w if w != 0 else geometric_mean for w in w_norm]

    return dict(zip(kept_codons, w_norm))


def _load_genome_sequences(genome_fasta_path):
    return {record.id: str(record.seq).upper() for record in SeqIO.parse(genome_fasta_path, 'fasta')}


def _anticodon_from_position(seqid, attrs, genome_seqs):
    """Resolve an anticodon given only as a genomic position
    (`anticodon=(pos:...)`) by extracting the real sequence there - the
    fallback path for tRNAscan-SE-sourced annotations (see the module-level
    comment above `_ANTICODON_POS_RE`). Returns None if `genome_seqs` wasn't
    given, the position can't be parsed, or `seqid` isn't in it."""
    if genome_seqs is None or seqid not in genome_seqs:
        return None
    match = _ANTICODON_POS_COMPLEMENT_RE.search(attrs)
    is_complement = match is not None
    if match is None:
        match = _ANTICODON_POS_RE.search(attrs)
    if match is None:
        return None
    start, end = int(match.group(1)), int(match.group(2))
    anticodon = genome_seqs[seqid][start - 1:end]  # GFF coords: 1-based, inclusive
    if len(anticodon) != 3:
        return None
    return _revcomp(anticodon) if is_complement else anticodon


def _parse_trna_gene_counts(gff_path, genome_seqs=None):
    """Return a 64-length tRNA gene copy number vector, indexed by
    CODON_ORDER, counted from a real GFF3's `tRNA` features.

    Tries the anticodon-embedded-in-text styles first (`Note=tRNA-<AA>(<anticodon>)`,
    `product=tRNA-<AA>(<anticodon>)` - RefSeq's own curated style, e.g. real
    E. coli annotations), falling back to resolving a position-only
    anticodon (`anticodon=(pos:...)`, tRNAscan-SE's own raw output style -
    confirmed this session on a real archaeon assembly) from `genome_seqs`
    (see `_load_genome_sequences`) if given.
    """
    counts = {codon: 0 for codon in CODON_ORDER}
    unresolved = 0
    with open(gff_path, encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9 or parts[2] != 'tRNA':
                continue
            seqid, attrs = parts[0], parts[8]
            match = _NOTE_ANTICODON_RE.search(attrs) or _PRODUCT_ANTICODON_RE.search(attrs)
            if match:
                anticodon = match.group(1).replace('U', 'T')
            else:
                anticodon = _anticodon_from_position(seqid, attrs, genome_seqs)
                if anticodon is None:
                    unresolved += 1
                    continue
            codon = _revcomp(anticodon)
            if codon in counts:
                counts[codon] += 1

    if unresolved and genome_seqs is None:
        raise ValueError(
            f"{unresolved} tRNA feature(s) in '{gff_path}' give the anticodon only as a genomic "
            "position (a tRNAscan-SE-sourced annotation style), not embedded directly in the "
            "feature's own attributes - pass genome_fasta_path so it can be resolved from the "
            "real sequence there, or these tRNA genes will be silently undercounted."
        )
    return [counts[codon] for codon in CODON_ORDER]


def derive_tai_weights_from_gff(gff_path, kingdom, genetic_code_num=None, genome_fasta_path=None):
    """Derive real tRNA Adaptation Index weights for an organism, from the
    tRNA gene features in its own GFF3 annotation (see
    eso.ncbi_genome.fetch_genome_package) - no reference-gene selection
    needed at all (unlike CAI - see eso.codon_usage.derive_table_from_genome),
    since tAI is grounded directly in genomic tRNA gene copy number rather
    than in codon usage of any gene set.

    `kingdom`: `'prokaryote'` or `'eukaryote'` - dos Reis et al.'s own
    `sking` parameter (bacteria get one additional wobble rule, for the ATA/
    Ile codon, that eukaryotes don't). Not reliably inferable from the
    genetic code table alone, so required explicitly.

    `genetic_code_num`: optional; if given, Trp/Met (single-codon amino
    acids, excluded from tAI's own scope) and any stop codons for that
    table are cross-checked against dos Reis' own hardcoded stop-codon/
    Met positions, and a ValueError is raised on a mismatch (e.g. a
    genetic code with reassigned stop codons) rather than silently
    returning a weight table for the wrong codons.

    `genome_fasta_path`: optional; needed only if this organism's GFF gives
    some tRNA anticodons as a genomic position rather than embedding the
    triplet directly (confirmed this session: real, common for
    tRNAscan-SE-sourced annotations, e.g. GCF_000091665.1 - RefSeq's own
    curated annotations, e.g. GCF_000005845.2/E. coli, don't need this).
    Omit it and this raises a clear error instead of silently undercounting,
    if it turns out to be needed for this particular GFF.

    Returns
    -------
    dict of {codon: weight}, one entry per codon except the 3 standard stop
    codons and ATG/Met - pass to `build_tai_score_fn` to get a scoring
    function usable as `optimization_engine`'s `custom_score_fn` argument.
    """
    if kingdom not in _KINGDOM_TO_SKING:
        raise ValueError(f"kingdom must be 'prokaryote' or 'eukaryote', got {kingdom!r}.")

    if genetic_code_num is not None:
        table = unambiguous_dna_by_id[genetic_code_num]
        real_stops = set(table.stop_codons)
        assumed_stops = {CODON_ORDER[10], CODON_ORDER[11], CODON_ORDER[14]}
        if real_stops != assumed_stops:
            raise ValueError(
                f"Genetic code table {genetic_code_num} has stop codons {sorted(real_stops)}, "
                f"not the standard {sorted(assumed_stops)} this tAI implementation assumes - "
                "this genetic code isn't supported by derive_tai_weights_from_gff."
            )

    genome_seqs = _load_genome_sequences(genome_fasta_path) if genome_fasta_path is not None else None
    trna_counts = _parse_trna_gene_counts(gff_path, genome_seqs=genome_seqs)
    if sum(trna_counts) == 0:
        raise ValueError(
            f"No tRNA gene features found in '{gff_path}' - tAI weights can't be derived "
            "without real tRNA gene copy numbers."
        )
    return _get_ws(trna_counts, sking=_KINGDOM_TO_SKING[kingdom])


def build_tai_score_fn(tai_weights):
    """Build a `(seq) -> float` scoring function from tAI weights (see
    `derive_tai_weights_from_gff`), directly usable as
    `optimization_engine`'s `custom_score_fn` argument.

    Computes the standard tAI sequence-level score: the geometric mean of
    each codon's own weight across the sequence (codons this table has no
    weight for - stop codons, Met - are skipped, matching how a coding
    sequence's stop codon and single Met codons carry no adaptive-choice
    information to score).
    """

    def score_fn(seq):
        seq = seq.upper()
        codons = [seq[i:i + 3] for i in range(0, len(seq) - 2, 3)]
        log_weights = [math.log(tai_weights[c]) for c in codons if c in tai_weights]
        if not log_weights:
            return 0.0
        return math.exp(sum(log_weights) / len(log_weights))

    return score_fn


def _rscu(reference_codons, aa_to_codons):
    """Relative Synonymous Codon Usage: for each codon, its observed count
    divided by the count expected if its amino acid's synonymous codons were
    used uniformly (`count(aa) / degeneracy(aa)`) - the standard, textbook
    formula (matches gtAI's own `CA_RSCU.RSCU`, verified directly against
    that source rather than assumed from memory). `reference_codons`: a
    flat list of codons (see `eso.codon_usage._select_reference_codons`) -
    only counts matter, not original sequence order.
    """
    codon_counts = Counter(reference_codons)
    rscu = {}
    for degenerate_codons in aa_to_codons.values():
        aa_count = sum(codon_counts.get(c, 0) for c in degenerate_codons)
        if aa_count == 0:
            continue
        expected = aa_count / len(degenerate_codons)
        for c in degenerate_codons:
            rscu[c] = codon_counts.get(c, 0) / expected
    return rscu


def derive_species_optimized_tai_weights(cds_fasta_path, gff_path, kingdom, genetic_code_num=None,
                                          genome_fasta_path=None, min_len_codons=100, top_perc=0.05,
                                          min_gene_count=50, seed=0):
    """Derive species-optimized tRNA Adaptation Index weights - reimplements
    the real gtAI (Anwar et al. 2023) algorithm: instead of dos Reis' fixed
    generic wobble parameters (`derive_tai_weights_from_gff`), optimize the
    same formula's 5 free wobble parameters (`_get_ws`'s `s[4:9]`) to
    maximize the Spearman correlation between the resulting per-codon tAI
    weights and RSCU computed over this organism's own highly-expressed
    reference genes (see module docstring for why this reimplements gtAI's
    published method rather than wrapping the `gtAI` package itself, which
    crashes on real data).

    Reuses `eso.codon_usage._select_reference_codons` for the reference set
    (the identical recipe gtAI's own `ENc_calc` uses, including this
    codebase's own length-floor fix for the short-gene contamination gtAI's
    unfixed version is equally exposed to) and this module's own `_get_ws`
    for scoring each candidate wobble-parameter vector (correct for every
    codon-degeneracy class, unlike gtAI's own `abs_Wi`).

    `seed`: `differential_evolution`'s own RNG seed, for reproducible output
    given the same inputs - this is a stochastic global optimizer, not a
    deterministic calculation.

    Returns
    -------
    dict of {codon: weight}, same shape as `derive_tai_weights_from_gff`.

    Raises
    ------
    ValueError
        If there aren't enough reference genes or tRNA genes to optimize
        against meaningfully (too few genes/RSCU-Wi overlap for a Spearman
        correlation to mean anything), or if the optimizer's result is
        somehow non-finite - callers wanting an automatic fallback to the
        generic weights on any such failure should use `derive_tai_weights`
        instead of calling this directly.
    """
    from eso.codon_usage import CustomCodonTableFileError, _select_reference_codons, detect_genetic_code_num_from_gff

    if kingdom not in _KINGDOM_TO_SKING:
        raise ValueError(f"kingdom must be 'prokaryote' or 'eukaryote', got {kingdom!r}.")
    if genetic_code_num is None:
        genetic_code_num = detect_genetic_code_num_from_gff(gff_path)

    try:
        aa_to_codons, reference_codons = _select_reference_codons(
            cds_fasta_path, genetic_code_num, min_len_codons, top_perc, min_gene_count)
    except CustomCodonTableFileError as e:
        raise ValueError(str(e)) from e

    rscu = _rscu(reference_codons, aa_to_codons)
    if len(rscu) < 10:
        raise ValueError(
            f"Only {len(rscu)} codons have a defined RSCU value from the reference gene set - "
            "too few to meaningfully optimize wobble parameters against. Try a larger genome, "
            "or use derive_tai_weights_from_gff's fixed generic weights instead."
        )

    genome_seqs = _load_genome_sequences(genome_fasta_path) if genome_fasta_path is not None else None
    trna_counts = _parse_trna_gene_counts(gff_path, genome_seqs=genome_seqs)
    if sum(trna_counts) == 0:
        raise ValueError(
            f"No tRNA gene features found in '{gff_path}' - tAI weights can't be derived "
            "without real tRNA gene copy numbers."
        )
    sking = _KINGDOM_TO_SKING[kingdom]

    def neg_fitness(x):
        candidate_weights = _get_ws(trna_counts, sking, s=(0.0, 0.0, 0.0, 0.0, *x))
        shared_codons = [c for c in candidate_weights if c in rscu]
        rscu_values = [rscu[c] for c in shared_codons]
        wi_values = [candidate_weights[c] for c in shared_codons]
        corr, _ = spearmanr(rscu_values, wi_values)
        return -corr if corr == corr else 0.0  # corr can be NaN if either side is constant

    result = differential_evolution(neg_fitness, bounds=[(0.0, 1.0)] * 5, seed=seed, polish=True)
    if not np.isfinite(result.fun):
        raise ValueError(
            "The wobble-parameter optimizer didn't converge to a usable result (non-finite "
            "fitness) - this genome/tRNA gene set may be too small or unusual for species-"
            "optimized tAI. Use derive_tai_weights_from_gff's fixed generic weights instead."
        )

    return _get_ws(trna_counts, sking, s=(0.0, 0.0, 0.0, 0.0, *(float(v) for v in result.x)))


def derive_tai_weights(gff_path, kingdom, cds_fasta_path=None, genetic_code_num=None,
                       genome_fasta_path=None, prefer_species_optimized=True, **selection_kwargs):
    """Recommended entry point for tAI weights: tries species-optimized
    weights first (`derive_species_optimized_tai_weights`, reimplementing
    gtAI's real algorithm) when `cds_fasta_path` is given and
    `prefer_species_optimized=True`, falling back to dos Reis' fixed
    generic weights (`derive_tai_weights_from_gff`) on any failure there -
    e.g. too few reference/tRNA genes, or the optimizer failing to converge
    on an unusual/tiny genome. A warning explains the fallback when it
    happens, matching `eso.optimize`'s own style for an unrecognized
    organism.

    `**selection_kwargs`: forwarded to `derive_species_optimized_tai_weights`
    (`min_len_codons`, `top_perc`, `min_gene_count`, `seed`) - ignored if the
    species-optimized path isn't attempted.

    Returns
    -------
    dict of {codon: weight}, same shape as both underlying functions.
    """
    if cds_fasta_path is not None and prefer_species_optimized:
        try:
            return derive_species_optimized_tai_weights(
                cds_fasta_path, gff_path, kingdom, genetic_code_num=genetic_code_num,
                genome_fasta_path=genome_fasta_path, **selection_kwargs)
        except Exception as e:
            warnings.warn(
                f"Species-optimized tAI weights couldn't be derived ({e!r}) - falling back to "
                "dos Reis' fixed generic weights instead.",
                stacklevel=2,
            )

    return derive_tai_weights_from_gff(
        gff_path, kingdom, genetic_code_num=genetic_code_num, genome_fasta_path=genome_fasta_path)
