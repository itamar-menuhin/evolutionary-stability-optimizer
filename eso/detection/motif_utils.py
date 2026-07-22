"""Build Bio.motifs Motif objects from a plain IUPAC consensus string (e.g.
"GATC", "CCWGG"), instead of requiring a full MEME-format PSSM file.

This is the easy on-ramp for custom motifs: most methylation/restriction
motifs are naturally described this way (REBASE, NEB, and the primary
literature all give them as a consensus sequence with ambiguity codes, not as
a position-probability matrix) - so there's no need to hand-author a MEME
file just to check for one.
"""

from Bio import motifs

#: IUPAC nucleotide ambiguity codes -> the bases each one represents.
IUPAC_NUCLEOTIDE_CODES = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC',
    'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT',
}


def motif_from_consensus(name, consensus, pseudocount=1):
    """Build a single Bio.motifs Motif from an IUPAC consensus string.

    Each position scores an exact match to one of its allowed bases highest;
    ambiguity codes (e.g. `W` = A or T) split probability evenly between the
    bases they represent. A small `pseudocount` (added to every base at every
    position, allowed or not) keeps every probability nonzero, so scoring
    never divides by zero - it does not meaningfully weaken the requirement
    that a real match stick to the allowed bases.

    Parameters
    ----------
    name: str
        Name for the resulting motif (appears as `matching_motif` in
        eso.detection.methylation.find_motif_sites' output).
    consensus: str
        An IUPAC nucleotide string, e.g. "GATC" (exact) or "CCWGG"
        (W = A or T).
    pseudocount: int or float
        Added to every base's raw count at every position before
        normalizing. Higher values make the motif more tolerant of
        near-matches; the default (1, against an allowed-base count of 100)
        is close to an exact-match requirement.

    Returns
    -------
    A Bio.motifs Motif, usable anywhere a MEME-file-loaded motif is (e.g.
    eso.detection.methylation.find_motif_sites' `relevant_motifs`).
    """
    consensus = consensus.strip().upper()
    if not consensus:
        raise ValueError("consensus can't be empty.")

    counts = {letter: [] for letter in 'ACGT'}
    for position, char in enumerate(consensus):
        allowed = IUPAC_NUCLEOTIDE_CODES.get(char)
        if allowed is None:
            raise ValueError(
                f"'{char}' at position {position + 1} of '{consensus}' isn't a recognized "
                f"IUPAC nucleotide code. Valid codes: A, C, G, T, or an ambiguity code "
                f"({', '.join(sorted(c for c in IUPAC_NUCLEOTIDE_CODES if len(IUPAC_NUCLEOTIDE_CODES[c]) > 1))})."
            )
        for letter in 'ACGT':
            counts[letter].append((100 if letter in allowed else 0) + pseudocount)

    motif = motifs.Motif(counts=counts)
    motif.name = name
    return motif


def motifs_from_consensus_dict(consensus_by_name, pseudocount=1):
    """Build a list of Motifs from a plain {name: consensus_string} dict -
    e.g. `motifs_from_consensus_dict({"my_site": "GANTC"})`. This, plus
    `motif_from_consensus`, is all that's needed to define custom motifs
    without a MEME file.
    """
    return [motif_from_consensus(name, consensus, pseudocount=pseudocount) for name, consensus in consensus_by_name.items()]
