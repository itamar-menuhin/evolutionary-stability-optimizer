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
textbook source of unwanted "cryptic" transcription. **Caveat**: a real
sigma70 promoter needs BOTH hexamers at roughly the right spacing
(~17±1bp apart) - this module flags each hexamer independently (this
detector has no concept of "two motifs at a specific spacing"), so an
isolated hit is much weaker evidence than a "both boxes, correctly spaced"
finding would be. Treat isolated hits as a coarse, conservative screen, not
a confirmed cryptic promoter.

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

from eso.detection.motif_utils import motif_from_consensus

#: name -> IUPAC consensus sequence
COMMON_MOTIFS = {
    "dam": "GATC",
    "dcm": "CCWGG",
    "shine_dalgarno": "AGGAGG",
    "sigma70_minus35": "TTGACA",
    "sigma70_minus10": "TATAAT",
}


def load_common_motifs(names=None):
    """Return a list of Bio.motifs Motif objects for the requested common
    motifs (default: all of them - see COMMON_MOTIFS for the full list).

    `names` is case-insensitive; see COMMON_MOTIFS for the available keys.
    """
    if names is None:
        names = list(COMMON_MOTIFS)

    motifs_list = []
    for name in names:
        key = name.strip().lower()
        if key not in COMMON_MOTIFS:
            raise ValueError(f"Unknown common motif {name!r}; choose from {sorted(COMMON_MOTIFS)}")
        motifs_list.append(motif_from_consensus(key, COMMON_MOTIFS[key]))
    return motifs_list
