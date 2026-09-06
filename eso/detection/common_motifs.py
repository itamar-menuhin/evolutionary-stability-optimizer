"""A small, ready-to-use library of commonly-referenced DNA motifs worth
checking for/avoiding when designing an engineered sequence - not just
methylation, but a few other properties that turn up across the synthetic
biology literature.

**Methylation** (E. coli's two systems, since E. coli is already a
first-class host here - see eso.codon_usage's bundled "e_coli" support):
- `dam`: GATC, N6-methyladenine (essentially universal across E. coli
  strains and many other Gammaproteobacteria).
- `dcm`: CCWGG (W = A or T), C5-methylcytosine on the internal C.
Sources: NEB, "Dam and Dcm Methylases of E. coli"; EcoSal Plus, "DNA
Methylation" (doi.org/10.1128/ecosalplus.esp-0003-2013).

**Cryptic ribosome binding**:
- `shine_dalgarno`: AGGAGG, the canonical bacterial ribosome-binding-site
  consensus. A copy of this sequence occurring *inside* a coding region
  (not at an intended start codon) is a known source of unintended internal
  translation initiation.
Source: Shine-Dalgarno sequences are measurably depleted from within
bacterial coding sequences, consistent with selection against this exact
risk (see e.g. Mol. Biol. Evol. 35(10):2487, and PMC6107199).

**Cryptic bacterial (sigma70) promoter elements**:
- `sigma70_minus35`: TTGACA, the -35 hexamer consensus.
- `sigma70_minus10`: TATAAT, the -10 ("Pribnow box") hexamer consensus.
These are the two core hexamers recognized by E. coli's housekeeping sigma
factor; an accidental occurrence of either inside a coding sequence is a
textbook source of unwanted "cryptic" transcription. A real sigma70
promoter needs BOTH hexamers at roughly the right spacing (~17±1bp apart) -
`find_sigma70_promoter_pairs` below checks specifically for that (strand-
aware, correctly-spaced pairs only), which is much stronger evidence of a
genuine cryptic promoter than either hexamer scored in isolation via the
ordinary common_motifs/find_motif_sites path (which has no concept of "two
motifs at a specific spacing" - ordinary common_motifs hits for these two
names should still be treated as a coarse, conservative screen on their
own).

**Not included, and why**: transcription terminators (rho-independent
terminators are a hairpin + poly-U structure - a secondary-structure
property, not a fixed linear sequence motif this PSSM-based approach can
represent) and the eukaryotic Kozak sequence (a *desired* translation-
initiation context to match near a real start codon, not something to
avoid - a different problem from what this detector/module is for).
Restriction enzyme sites are also not duplicated here - DNAChisel (already
a dependency of this project) already bundles a comprehensive registry via
`dnachisel.list_common_enzymes()` / `dnachisel.EnzymeSitePattern`, usable
directly with `AvoidPattern` during optimization.

For anything beyond what's bundled here, REBASE (rebase.neb.com) is the
standard reference database for restriction/methylation motifs across
organisms - use `eso.detection.motif_utils.motif_from_consensus` to turn
any REBASE-style (or other literature) consensus sequence into a usable
motif, the same way everything in this module is built.
"""

import pandas as pd

from eso.detection.motif_utils import motif_from_consensus

#: name -> IUPAC consensus sequence
COMMON_MOTIFS = {
    "dam": "GATC",
    "dcm": "CCWGG",
    "shine_dalgarno": "AGGAGG",
    "sigma70_minus35": "TTGACA",
    "sigma70_minus10": "TATAAT",
}

#: columns for find_sigma70_promoter_pairs' output.
SIGMA70_PAIR_COLUMNS = [
    'minus35_start', 'minus35_end', 'minus10_start', 'minus10_end', 'strand', 'spacing', 'combined_score',
]


def load_common_motifs(names=None, pseudocount=1):
    """Return a list of Bio.motifs Motif objects for the requested common
    motifs (default: all of them - see COMMON_MOTIFS for the full list).

    `names` is case-insensitive; see COMMON_MOTIFS for the available keys.
    `pseudocount`: see motif_from_consensus - higher values loosen how
    close a match needs to be to count, for all requested motifs at once.
    Previously there was no way to adjust this for the bundled common
    motifs specifically (only for a fully custom motif built directly via
    motif_from_consensus) without editing this module's source.
    """
    if names is None:
        names = list(COMMON_MOTIFS)

    motifs_list = []
    for name in names:
        key = name.strip().lower()
        if key not in COMMON_MOTIFS:
            raise ValueError(f"Unknown common motif {name!r}; choose from {sorted(COMMON_MOTIFS)}")
        motifs_list.append(motif_from_consensus(key, COMMON_MOTIFS[key], pseudocount=pseudocount))
    return motifs_list


def _single_hexamer_hits(motif, seq):
    """Every position (both strands) where `motif`'s PSSM scores above
    random-chance probability (> 0) - the same threshold
    eso.detection.methylation.find_motif_sites uses. Returns a list of
    (start, end, strand, score) tuples; strand kept explicit (unlike
    find_motif_sites' collapsed output) since pairing below needs it.
    """
    length = len(motif)
    valid = len(seq) - length + 1
    if valid <= 0:
        return []

    fwd_scores = motif.pssm.calculate(seq)
    rev_scores = motif.pssm.reverse_complement().calculate(seq)

    hits = []
    for i in range(valid):
        if fwd_scores[i] > 0:
            hits.append((i, i + length - 1, '+', float(fwd_scores[i])))
        if rev_scores[i] > 0:
            hits.append((i, i + length - 1, '-', float(rev_scores[i])))
    return hits


def find_sigma70_promoter_pairs(seq, spacing_min=16, spacing_max=18):
    """Detect PAIRS of sigma70_minus35 and sigma70_minus10 hexamer hits at
    the biologically correct spacing on the same strand - much stronger
    evidence of a real cryptic promoter than either hexamer scored in
    isolation (see this module's own docstring: a real sigma70 promoter
    needs both hexamers, correctly spaced ~17±1bp apart, which
    eso.detection.methylation.find_motif_sites' per-position independent
    scoring has no way to express - it can only ever report isolated hits).

    spacing_min/spacing_max: inclusive bounds, in nucleotides, on the gap
    between the end of the -35 box and the start of the -10 box (measuring
    "downstream" in the promoter's own reading direction, so on the reverse
    strand this is a gap in the direction of DECREASING sequence
    coordinate) - default 16-18, i.e. the documented ~17±1bp consensus
    spacing.

    Returns a dataframe (see SIGMA70_PAIR_COLUMNS), highest combined-score
    pair first. Empty (not an error) if no correctly-spaced pair is found -
    the individual hexamers may still be worth checking separately via the
    ordinary common_motifs path if a screen for isolated hits is wanted too.
    """
    minus35_motif, minus10_motif = load_common_motifs(['sigma70_minus35', 'sigma70_minus10'])
    minus35_hits = _single_hexamer_hits(minus35_motif, seq)
    minus10_hits = _single_hexamer_hits(minus10_motif, seq)

    pairs = []
    for m35_start, m35_end, m35_strand, m35_score in minus35_hits:
        for m10_start, m10_end, m10_strand, m10_score in minus10_hits:
            if m35_strand != m10_strand:
                continue
            if m35_strand == '+':
                spacing = m10_start - m35_end - 1
            else:
                spacing = m35_start - m10_end - 1
            if spacing_min <= spacing <= spacing_max:
                pairs.append({
                    'minus35_start': m35_start, 'minus35_end': m35_end,
                    'minus10_start': m10_start, 'minus10_end': m10_end,
                    'strand': m35_strand, 'spacing': spacing,
                    'combined_score': m35_score + m10_score,
                })

    df = pd.DataFrame(pairs, columns=SIGMA70_PAIR_COLUMNS)
    if not df.empty:
        df = df.sort_values('combined_score', ascending=False).reset_index(drop=True)
    return df
