"""Small sequence/region helpers shared across detection and optimization."""

import pandas as pd

COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

_VALID_BASES = frozenset("ACGT")


class InvalidSequenceError(Exception):
    """A sequence contains a character that isn't A, C, G, or T.

    Raised with a plain-English message identifying exactly where, so it can
    be fixed without knowing anything about DNAChisel or IUPAC ambiguity
    codes.
    """


def first_invalid_base(seq):
    """Returns (0-indexed position, character) of the first letter in `seq`
    that isn't A, C, G, or T (matched case-insensitively - downstream code
    uppercases the sequence, see eso.pipeline.backend's
    `str(record.seq).upper()`), or None if the whole sequence is valid.
    """
    for i, ch in enumerate(seq):
        if ch.upper() not in _VALID_BASES:
            return i, ch
    return None


def validate_dna_alphabet(seq, sequence_label="This sequence"):
    """Raises InvalidSequenceError with a plain-English message if `seq`
    contains anything other than A/C/G/T (case-insensitive) - e.g. an IUPAC
    ambiguity code (N, R, Y, ...) or a stray non-DNA character.

    Nothing downstream validates this itself, and the failure modes without
    it are bad: confirmed directly that a codon of a single ambiguity code
    (e.g. "NNN") crashes with an unhandled `dnachisel.biotools.TranslationError`
    once EnforceTranslation is in play, and a GC-content-only run (no
    translation constraint reached yet) crashes with a bare `KeyError`
    instead - both deep inside DNAChisel, with no eso-level message at all.
    Separately, eso.io_utils.exclusion_gc_tester's GC-content check doesn't
    crash on an ambiguity code, but silently miscalculates instead: it
    treats an ambiguous base (e.g. "R" = A-or-G) as neither G nor C, giving a
    confidently wrong GC fraction rather than flagging the base as unknown -
    confirmed directly (a sequence of all "R"s "measured" as ~4.8% GC despite
    GC content being genuinely undefined for it). See
    docs/detector-comparisons.md for the full investigation.
    """
    if not seq:
        # Confirmed directly: an empty sequence reaches DNAChisel's own
        # EnforceTranslation.restrict_nucleotides with a zero-length ORF
        # region (from optimize.py's default `orf_regions` computation) and
        # crashes with a bare, unhandled `IndexError: string index out of
        # range` deep inside DNAChisel - no eso-level message at all.
        raise InvalidSequenceError(f"{sequence_label} is empty - there's nothing to optimize.")

    invalid = first_invalid_base(seq)
    if invalid is not None:
        position, char = invalid
        raise InvalidSequenceError(
            f"{sequence_label} contains a letter other than A, C, G, T at position {position + 1} "
            f"({char!r}) - please check for typos or ambiguous bases."
        )


def reverse_complement_seq(seq):
    return ''.join(COMPLEMENT[x] for x in seq[::-1])


def add_backward_sites(df):
    """Duplicate each row of a {sequence, start, end} dataframe with its reverse complement."""
    df_forward = df.copy()
    df_backward = df.copy()
    df_backward.loc[:, 'sequence'] = df_backward.sequence.apply(reverse_complement_seq)
    return pd.concat([df_forward, df_backward], ignore_index=False)


def shorten_sequences(df):
    df_short = df.copy()
    df_short.loc[:, 'sequence'] = df_short.sequence.apply(lambda x: x[:-1])
    df_short.loc[:, 'end'] = df_short.end.apply(lambda x: x - 1)
    return df_short


def parse_region(region_string):
    """Parse a region string like "start_1-end_1,start_2-end_2,..." (1-indexed, inclusive)
    into a list of 0-indexed (start, end) tuples. Returns () for '' or 'None', 'error' if malformed.
    """
    if region_string in ('', 'None'):
        return ()

    region_list = region_string.split(',')
    try:
        return [
            (int(region.strip().split('-')[0]) - 1, int(region.strip().split('-')[1]))
            for region in region_list
        ]
    except (ValueError, IndexError):
        return 'error'
