"""Codon usage bias (CUB) tables for hosts not covered by python-codon-tables."""

import re
from collections import Counter, defaultdict
from importlib import resources
from os import path

import numpy as np
import pandas as pd
from Bio import SeqIO, SeqUtils
from Bio.Data.CodonTable import unambiguous_dna_by_id
from Bio.Seq import Seq


def cub_c1():
    codon_usage_table = {
        'A': {'GCA': 0.13718, 'GCC': 0.432097, 'GCG': 0.264336, 'GCT': 0.166387},
        '*': {'TAA': 0.265377, 'TAG': 0.181818, 'TGA': 0.552805}, 'W': {'TGG': 1},
        'R': {'AGA': 0.0973128, 'AGG': 0.153007, 'CGA': 0.128389, 'CGC': 0.304978, 'CGG': 0.209737,
              'CGT': 0.106575}, 'N': {'AAC': 0.792265, 'AAT': 0.207735},
        'D': {'GAC': 0.718222, 'GAT': 0.281778}, 'C': {'TGC': 0.710465, 'TGT': 0.289535},
        'Q': {'CAA': 0.361039, 'CAG': 0.638961}, 'E': {'GAA': 0.250212, 'GAG': 0.749788},
        'G': {'GGA': 0.142656, 'GGC': 0.506165, 'GGG': 0.184803, 'GGT': 0.166376},
        'H': {'CAC': 0.650846, 'CAT': 0.349154},
        'I': {'ATA': 0.0923118, 'ATC': 0.673554, 'ATT': 0.234134},
        'L': {'CTA': 0.0641845, 'CTC': 0.345609, 'CTG': 0.309073, 'CTT': 0.132086, 'TTA': 0.030553,
              'TTG': 0.118494}, 'K': {'AAA': 0.180592, 'AAG': 0.819408}, 'M': {'ATG': 1},
        'F': {'TTC': 0.672764, 'TTT': 0.327236},
        'P': {'CCA': 0.167124, 'CCC': 0.336275, 'CCG': 0.326667, 'CCT': 0.169935},
        'Y': {'TAC': 0.757778, 'TAT': 0.242222},
        'S': {'AGC': 0.234258, 'AGT': 0.0671675, 'TCA': 0.103458, 'TCC': 0.220449, 'TCG': 0.247293,
              'TCT': 0.127375},
        'T': {'ACA': 0.154802, 'ACC': 0.406355, 'ACG': 0.307613, 'ACT': 0.131231},
        'V': {'GTA': 0.0871417, 'GTC': 0.460197, 'GTG': 0.289844, 'GTT': 0.162817},
    }
    return codon_usage_table


def cub_kompas():
    """Codon usage table for Komagataella phaffii (Pichia pastoris)."""
    cub_full = {
        'Ala': {'GCA': 0.275098, 'GCC': 0.244931, 'GCG': 0.0786548, 'GCT': 0.401316},
        'Arg': {'AGA': 0.455639, 'AGG': 0.181255, 'CGA': 0.119597, 'CGC': 0.0506672, 'CGG': 0.0522151,
                'CGT': 0.140626},
        'Asn': {'AAC': 0.465812, 'AAT': 0.534188},
        'Asp': {'GAC': 0.382762, 'GAT': 0.617238},
        'Cys': {'TGC': 0.377246, 'TGT': 0.622754},
        'Gln': {'CAA': 0.603655, 'CAG': 0.396345},
        'Glu': {'GAA': 0.594953, 'GAG': 0.405047},
        'Gly': {'GGA': 0.364371, 'GGC': 0.155831, 'GGG': 0.12228, 'GGT': 0.357518},
        'His': {'CAC': 0.382548, 'CAT': 0.617452},
        'Ile': {'ATA': 0.236221, 'ATC': 0.297661, 'ATT': 0.466118},
        'Leu': {'CTA': 0.122977, 'CTC': 0.0842546, 'CTG': 0.152862, 'CTT': 0.170539, 'TTA': 0.178334,
                'TTG': 0.291033},
        'Lys': {'AAA': 0.519264, 'AAG': 0.480736},
        'Met': {'ATG': 1},
        'Phe': {'TTC': 0.42216, 'TTT': 0.57784},
        'Pro': {'CCA': 0.378486, 'CCC': 0.180179, 'CCG': 0.102394, 'CCT': 0.338942},
        'Ser': {'AGC': 0.104434, 'AGT': 0.158704, 'TCA': 0.208354, 'TCC': 0.174084, 'TCG': 0.0945539,
                'TCT': 0.259871},
        'Thr': {'ACA': 0.277392, 'ACC': 0.237645, 'ACG': 0.122946, 'ACT': 0.362017},
        'Trp': {'TGG': 1},
        'Tyr': {'TAC': 0.48627, 'TAT': 0.51373},
        'Val': {'GTA': 0.178012, 'GTC': 0.217067, 'GTG': 0.216627, 'GTT': 0.388294},
        'END': {'TAA': 0.399841, 'TAG': 0.339428, 'TGA': 0.260731},
    }
    # SeqUtils.seq1('END') returns 'X' (undefined amino acid), not '*' (stop) -
    # without this special case, the stop-codon frequencies silently ended up
    # filed under the wrong key and were never used for stop-codon scoring.
    return {('*' if x == 'END' else SeqUtils.seq1(x)): cub_full[x] for x in cub_full}


def _load_bundled_csv_cub(data_filename):
    with resources.files("eso.data").joinpath(data_filename).open("r", encoding="utf-8") as handle:
        df = pd.read_csv(handle)

    codon_usage_table = df.groupby('aa').apply(
        lambda x: x.set_index('codon')['freq_within_aa'].to_dict(), include_groups=False
    ).to_dict()

    if '*' not in codon_usage_table:
        codon_usage_table['*'] = {'TAA': 0.33, 'TAG': 0.33, 'TGA': 0.34}

    return codon_usage_table


def cub_human_antibody_heavy_chain():
    """Codon usage table for human antibody heavy chain, from the iGEM 2025 dataset."""
    return _load_bundled_csv_cub("human-antibody-heavy-chain-codon-frequencies.csv")


def cub_human_antibody_light_chain():
    """Codon usage table for human antibody light chain, from the iGEM 2025 dataset."""
    return _load_bundled_csv_cub("human-antibody-light-chain-codon-frequencies.csv")


# Selecting one of these (or any other organism_name) only changes what
# codon-usage table optimization.py's CodonOptimize objective scores
# against - it does NOT change which mutation-rate model hotspot detection
# uses to score recombination/slippage risk. That model (see
# eso.detection.recombination.calc_recombination_score's docstring) is
# calibrated specifically for E. coli's RecA-mediated recombination
# machinery, and no equivalent published model exists for these other hosts
# at this resolution - detected sites' relative ranking still transfers
# reasonably, but their absolute risk scores don't reflect this organism's
# actual recombination/repair biology.
CODON_USAGE_TABLES = {
    'C1': cub_c1,
    'kompas': cub_kompas,
    'human_antibody_heavy_chain': cub_human_antibody_heavy_chain,
    'human_antibody_light_chain': cub_human_antibody_light_chain,
}


class CustomCodonTableFileError(Exception):
    """A custom codon-usage-table file failed to load or didn't behave as
    expected.

    Raised with a plain-English message aimed at someone bringing in a
    project-specific codon usage table (e.g. exported from a codon-analysis
    tool, or one of this library's own bundled tables reused standalone) -
    not a bioinformatics/pandas expert - the goal is that this exception's
    message alone is enough to fix the problem.
    """


def load_custom_codon_table_from_file(file_path):
    """Load a user-supplied codon-usage table from a CSV file.

    Expected columns (any order, extra columns ignored): `codon`, `aa`,
    `freq_within_aa` - the same shape this module's own bundled tables use
    (see eso/data/*.csv). `aa` is the one-letter amino acid code ('*' for
    stop), `codon` a DNA triplet, `freq_within_aa` that codon's relative
    usage among synonymous codons for that amino acid (need not sum to
    exactly 1 - DNAChisel normalizes internally).

    Mirrors eso.custom_score.load_custom_score_from_file's validation
    philosophy: fails eagerly, with a message meant to be read and acted on
    directly by whoever exported this table, not a Python/pandas expert.

    This is the mechanism behind `eso-optimize --codon-usage-table-file`,
    and is also usable directly from Python - pass the result as
    `optimization_engine`'s `codon_usage_table` argument.

    Returns
    -------
    dict of the form {'*': {'TAA': 0.33, ...}, 'K': {'AAA': ..., 'AAG': ...}, ...}
    """
    if not path.isfile(file_path):
        raise CustomCodonTableFileError(
            f"Can't find the codon-usage table file '{file_path}'. Check the path is correct.")

    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        raise CustomCodonTableFileError(
            f"'{file_path}' could not be read as a CSV file. The error was: {e!r}. Make sure it's "
            "a plain comma-separated file, not e.g. an Excel file saved with a .csv extension."
        ) from e

    missing_columns = {'codon', 'aa', 'freq_within_aa'} - set(df.columns)
    if missing_columns:
        raise CustomCodonTableFileError(
            f"'{file_path}' is missing the column(s) {sorted(missing_columns)}. Expected a CSV with "
            "(at least) the columns 'codon', 'aa', and 'freq_within_aa' - one row per codon, 'aa' "
            "the one-letter amino acid code ('*' for stop), 'freq_within_aa' that codon's relative "
            "usage among synonymous codons for that amino acid."
        )

    if df.empty:
        raise CustomCodonTableFileError(f"'{file_path}' has no data rows - it's an empty table.")

    invalid_codons = sorted({
        c for c in df['codon']
        if not (isinstance(c, str) and len(c) == 3 and set(c.upper()) <= set('ACGT'))
    })
    if invalid_codons:
        raise CustomCodonTableFileError(
            f"'{file_path}' has invalid entries in its 'codon' column: {invalid_codons[:5]!r}"
            f"{' (and more)' if len(invalid_codons) > 5 else ''} - every codon must be a 3-letter "
            "DNA triplet (A/C/G/T only)."
        )

    non_numeric = df[~df['freq_within_aa'].apply(lambda v: isinstance(v, (int, float)) and not isinstance(v, bool))]
    if not non_numeric.empty:
        bad_values = sorted({repr(v) for v in non_numeric['freq_within_aa']})
        raise CustomCodonTableFileError(
            f"'{file_path}' has non-numeric value(s) in its 'freq_within_aa' column: {bad_values[:5]}"
            f"{' (and more)' if len(bad_values) > 5 else ''} - every value there must be a plain number."
        )

    codon_usage_table = df.groupby('aa').apply(
        lambda x: x.set_index('codon')['freq_within_aa'].to_dict(), include_groups=False
    ).to_dict()

    # Matches this module's own bundled tables' behavior (see
    # _load_bundled_csv_cub) - a table with no stop-codon row would otherwise
    # leave every ORF's stop codon unscored, not fail loudly, so this is
    # filled in rather than left to surface as a confusing downstream effect.
    if '*' not in codon_usage_table:
        codon_usage_table['*'] = {'TAA': 0.33, 'TAG': 0.33, 'TGA': 0.34}

    # Checked directly (rather than only discovered deep inside DNAChisel
    # during a real optimization run) - the most common real mistake here is
    # a copy-paste error or a table exported with columns in a different
    # order than expected, silently filing a codon under the wrong amino
    # acid. Mirrors this module's own bundled-table test suite
    # (tests/test_codon_usage.py's test_every_codon_translates_to_the_amino_acid_it_is_filed_under).
    mismatches = [
        (codon, aa, str(Seq(codon).translate()))
        for aa, codons in codon_usage_table.items()
        for codon in codons
        if str(Seq(codon).translate()) != aa
    ]
    if mismatches:
        examples = ', '.join(
            f"{codon!r} filed under '{aa}' (actually translates to '{actual}')"
            for codon, aa, actual in mismatches[:5]
        )
        raise CustomCodonTableFileError(
            f"'{file_path}' has codon(s) filed under the wrong amino acid: {examples}"
            f"{' (and more)' if len(mismatches) > 5 else ''}. Check the 'aa' column matches each "
            "codon's actual translation."
        )

    return codon_usage_table


_TRANSL_TABLE_RE = re.compile(r'transl_table=(\d+)')


def detect_genetic_code_num_from_gff(gff_path):
    """Look up the NCBI genetic code table number to pass to
    `derive_table_from_genome`/`eso.tai.derive_tai_weights_from_gff`, via a
    majority vote over every CDS feature's own `transl_table=` qualifier in
    the GFF that came with the same genome package (see
    eso.ncbi_genome.fetch_genome_package) - more robust than asking the
    caller to know/guess the right table number (confirmed directly this
    session: real NCBI annotations for a bacterium and an archaeon both
    state `transl_table=11`, not always the same as the domain-of-life
    default one might assume)."""
    votes = Counter()
    with open(gff_path, encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9 or parts[2] != 'CDS':
                continue
            match = _TRANSL_TABLE_RE.search(parts[8])
            if match:
                votes[int(match.group(1))] += 1
    if not votes:
        raise CustomCodonTableFileError(
            f"Couldn't find any `transl_table=` qualifier on a CDS feature in '{gff_path}' - "
            "pass genetic_code_num explicitly instead of relying on auto-detection."
        )
    return votes.most_common(1)[0][0]


#: The Wright 1990 ENc formula's own fixed coefficients (the 2, 9, 1, 5, 3 in
#: _calc_enc's final expression) aren't free parameters - they ARE the
#: standard genetic code's own count of amino-acid families at each
#: degeneracy level (1-fold: Met+Trp=2; 2-fold: 9 families; 3-fold: Ile
#: alone=1; 4-fold: 5 families; 6-fold: Leu+Arg+Ser=3). Confirmed directly
#: (checked every NCBI genetic code table Biopython bundles, not assumed)
#: that this exact family-count structure holds only for table 1 (standard)
#: and table 11 (bacterial/archaeal/plant-plastid) - the two genuinely
#: supported by this module - and genuinely differs for every one of the
#: other 25 tables (e.g. table 2, vertebrate mitochondrial, has no 3-fold
#: class at all: AGA/AGG are stop codons there, not Arg). Using this formula
#: with a mismatched genetic code wouldn't crash - _calc_enc's own fallback
#: (1/degree) for an empty degeneracy class silently absorbs the gap - it
#: would silently compute a scientifically meaningless ENc value, biasing
#: the entire reference-gene selection this table/weight derivation depends
#: on, with no error or warning at all.
_ENC_FORMULA_DEGENERACY_COUNTS = {1: 2, 2: 9, 3: 1, 4: 5, 6: 3}


def _validate_enc_formula_applies(aa_to_codons, genetic_code_num):
    degeneracy_counts = dict(Counter(len(codons) for codons in aa_to_codons.values()))
    if degeneracy_counts != _ENC_FORMULA_DEGENERACY_COUNTS:
        raise CustomCodonTableFileError(
            f"genetic_code_num={genetic_code_num} has a different codon-degeneracy structure than "
            "the standard genetic code (tables 1 and 11) - the ENc (Effective Number of Codons) "
            "formula this reference-gene selection relies on is calibrated specifically for that "
            "structure and would silently compute a meaningless result here, not raise an error on "
            "its own. derive_table_from_genome and derive_species_optimized_tai_weights currently "
            "only support genetic_code_num 1 or 11."
        )


def _calc_enc(codons, aa_to_codons):
    """Effective Number of Codons (Wright 1990) - a purely sequence-intrinsic
    measure of codon bias needing no external annotation at all. Lower means
    more strongly biased; genuinely highly-expressed genes are, empirically,
    disproportionately the most strongly biased ones. Ported from STABLES'
    own `create_he.py` (identical across all 8 of its host organisms), which
    uses exactly this formula to select its own highly-expressed-gene proxy
    set for CAI/RCA weighting."""
    degeneracy_groups = defaultdict(list)
    for aa, degenerate_codons in aa_to_codons.items():
        counts = Counter(c for c in codons if c in degenerate_codons)
        if counts:
            freqs = sum((n / sum(counts.values())) ** 2 for n in counts.values())
            degeneracy_groups[len(degenerate_codons)].append(freqs)

    d_vec = {}
    for degree in (1, 2, 3, 4, 6):
        group = degeneracy_groups.get(degree, [])
        d_vec[degree] = (1 / degree) if not group else float(np.average(group))

    return 2 + 9 / d_vec[2] + 1 / d_vec[3] + 5 / d_vec[4] + 3 / d_vec[6]


def _select_reference_codons(cds_fasta_path, genetic_code_num, min_len_codons, top_perc, min_gene_count):
    """Shared reference-gene selection behind both `derive_table_from_genome`
    (Sharp & Li CAI weights) and `eso.tai.derive_species_optimized_tai_weights`
    (RSCU for its fitness function) - both need "this genome's own
    highly-expressed-like gene set", and must use the *identical* selection
    to be meaningfully comparable/combinable, not two subtly-different
    reimplementations of the same idea.

    See `derive_table_from_genome`'s docstring for what `min_len_codons`/
    `top_perc`/`min_gene_count` mean and why these particular defaults were
    chosen - unchanged by this refactor.

    Returns
    -------
    (aa_to_codons, reference_codons) - `aa_to_codons` maps each one-letter
    amino acid code to its list of synonymous codons for `genetic_code_num`;
    `reference_codons` is every codon (with repeats) across the selected
    reference genes, flattened - sufficient for any codon-count-based
    statistic (Sharp & Li weights, RSCU, ...), no need for original
    per-gene sequence order.
    """
    aa_to_codons = defaultdict(list)
    for codon, aa in unambiguous_dna_by_id[genetic_code_num].forward_table.items():
        aa_to_codons[aa].append(codon)
    _validate_enc_formula_applies(aa_to_codons, genetic_code_num)

    genes = []
    for record in SeqIO.parse(cds_fasta_path, 'fasta'):
        seq = str(record.seq).upper()
        codons = [seq[i:i + 3] for i in range(0, len(seq) - 2, 3)]
        codons = [c for c in codons if len(c) == 3 and set(c) <= set('ACGT')]
        if len(codons) < min_len_codons:
            continue
        genes.append((_calc_enc(codons, aa_to_codons), codons))

    if not genes:
        raise CustomCodonTableFileError(
            f"No gene in '{cds_fasta_path}' has at least min_len_codons={min_len_codons} "
            "codons - try a lower min_len_codons, or check this is really a CDS FASTA."
        )

    genes.sort(key=lambda item: item[0])
    n_reference = max(round(top_perc * len(genes)), min(min_gene_count, len(genes)))
    reference_codons = [codon for _, codons in genes[:n_reference] for codon in codons]

    return aa_to_codons, reference_codons


def derive_table_from_genome(cds_fasta_path, genetic_code_num=None, min_len_codons=100,
                              top_perc=0.05, min_gene_count=50):
    """Derive a real, organism-specific codon-usage table directly from a
    genome's own CDS sequences - an alternative to the 4 bundled tables,
    python-codon-tables' Kazusa-sourced species, or a hand-supplied CSV, for
    any organism you have (or can fetch via eso.ncbi_genome.fetch_genome_package)
    a CDS FASTA for.

    Unlike a plain genome-wide-average codon-frequency table (what every
    other bundled/third-party table this library touches actually is), this
    computes the real Sharp & Li (1987) Codon Adaptation Index construction:
    codon frequencies from a reference set of highly-expressed genes only,
    not the whole genome. The highly-expressed set is found via ENc
    (Effective Number of Codons, Wright 1990; see `_calc_enc`) - a purely
    sequence-intrinsic statistic needing no functional annotation at all,
    the same method STABLES' own `create_he.py` already uses in production
    across all 8 of its host organisms.

    `min_len_codons` (default 100) and `top_perc` (default 0.05, with a
    `min_gene_count` floor of 50 genes) were both empirically tuned and
    cross-validated this session across three organisms spanning different
    domains of life (a bacterium, a eukaryote, and an archaeon): with no
    length floor, short leader peptides/toxin-antitoxin genes/transposase
    fragments contaminate the reference set (spurious low ENc from small-
    sample noise, not real translational selection) - this affected even
    STABLES' own flagship E. coli host, not just unusual organisms. A
    parameter sweep (25-400 codons; 1%-30% cutoff) found a genuine, cross-
    organism-consistent stability minimum right around these defaults -
    below it, sample-size noise dominates; above it, real ribosomal protein
    genes (inherently short) get systematically excluded and the reference
    set drifts again. `min_gene_count` guards against a small/fragmentary
    genome where `top_perc` alone would pick too few genes to reliably
    estimate rarer codons' frequencies.

    `genetic_code_num`: an NCBI genetic code table number - required, and
    currently must be 11 (bacteria/archaea/plant plastids) or 1 (the
    standard/eukaryotic-nuclear code); the ENc formula this relies on
    (`_calc_enc`) is calibrated specifically for the codon-degeneracy
    structure those two tables share, and raises a clear error for any other
    table rather than silently returning a meaningless result (confirmed
    directly: every other NCBI genetic code table has a genuinely different
    degeneracy structure - e.g. table 2, vertebrate mitochondrial, has no
    3-fold-degenerate amino acid at all). This function only reads the CDS
    FASTA, not the GFF, so it can't auto-detect this itself; look it up
    first via `detect_genetic_code_num_from_gff(gff_path)`, using the GFF
    from the same `eso.ncbi_genome.fetch_genome_package` call.

    Returns
    -------
    dict of the form {'*': {'TAA': 0.33, ...}, 'K': {'AAA': ..., 'AAG': ...}, ...} -
    the same shape as this module's bundled tables and
    `load_custom_codon_table_from_file`'s return value, directly usable as
    `optimization_engine`'s `codon_usage_table` argument.
    """
    if genetic_code_num is None:
        raise CustomCodonTableFileError(
            "genetic_code_num wasn't given - pass it explicitly (11 for bacteria/archaea, "
            "1 for the standard/eukaryotic-nuclear code, etc.), or look it up first via "
            "detect_genetic_code_num_from_gff(gff_path) using the GFF that came with this "
            "same genome package."
        )
    aa_to_codons, reference_codons = _select_reference_codons(
        cds_fasta_path, genetic_code_num, min_len_codons, top_perc, min_gene_count)

    codon_counts = Counter(reference_codons)
    codon_usage_table = {}
    for aa, degenerate_codons in aa_to_codons.items():
        max_count = max((codon_counts.get(c, 0) for c in degenerate_codons), default=0)
        if max_count == 0:
            continue
        codon_usage_table[aa] = {c: codon_counts.get(c, 0) / max_count for c in degenerate_codons}

    if '*' not in codon_usage_table:
        codon_usage_table['*'] = {'TAA': 0.33, 'TAG': 0.33, 'TGA': 0.34}

    return codon_usage_table


