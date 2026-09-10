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
    """Reverse-complement a DNA sequence (A/C/G/T, case-insensitive - always
    returned uppercase, matching this codebase's convention of uppercasing
    sequences once at the point they're first read - see e.g.
    eso.pipeline.backend's `str(record.seq).upper()`).

    Every current caller only ever passes a substring of an already-validated
    (validate_dna_alphabet) sequence, but nothing enforced that here - a
    stray ambiguity code or non-DNA character would otherwise hit a raw,
    unhelpful `KeyError` deep in this function with no eso-level message at
    all, the exact failure mode validate_dna_alphabet's own docstring
    describes fixing elsewhere. Raises the same InvalidSequenceError here
    instead, so this function is safe to call on its own, not just implicitly
    trusted to only ever see pre-validated input.
    """
    try:
        return ''.join(COMPLEMENT[x.upper()] for x in seq[::-1])
    except KeyError as e:
        raise InvalidSequenceError(
            f"Cannot reverse-complement {seq!r}: contains a letter other than A, C, G, T "
            f"({e.args[0]!r})."
        ) from None


def add_backward_sites(df):
    """Duplicate each row of a {sequence, start, end} dataframe with its
    reverse complement.

    `ignore_index=True`: without it, the concatenated result has each
    original row index repeated twice (once per half) rather than a clean,
    unique index - harmless for every current caller (they only ever access
    columns by position/zip(), never by index label), but a real trap for
    any future or direct caller that does a label-based `.loc[idx]` lookup
    expecting one row back, which would silently get two instead. Hardened
    regardless of current reachability, not just left as a latent risk.
    """
    df_forward = df.copy()
    df_backward = df.copy()
    df_backward.loc[:, 'sequence'] = df_backward.sequence.apply(reverse_complement_seq)
    return pd.concat([df_forward, df_backward], ignore_index=True)


def shorten_sequences(df):
    """Drop the last character of each row's `sequence` and decrement `end`
    by one - used to build the "one deletion away" candidate set.

    Raises ValueError for any row whose sequence has length < 2: shortening
    it further would silently produce a degenerate empty-string pattern
    passed on to detection/DNAChisel code downstream instead of any
    genuinely useful "one deletion away" candidate. Every current caller
    only ever passes 16-17nt windows, so this isn't reachable today -
    hardened anyway rather than relying on that staying true forever.
    """
    if not df.empty and (df.sequence.str.len() < 2).any():
        raise ValueError("shorten_sequences requires every sequence to have length >= 2.")
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
        regions = []
        for region in region_list:
            parts = region.strip().split('-')
            # Confirmed a real, previously-silent bug: with `[0]`/`[1]` pulled
            # directly from an unchecked split(), an extra dash (e.g. a typo
            # like "10-20-30", or "10--5") silently parsed as (9, 20) -
            # dropping "-30"/producing a bogus negative end entirely
            # unflagged - rather than being caught here as malformed, exactly
            # the kind of input this function's own docstring promises to
            # reject. Confirmed directly before this fix.
            if len(parts) != 2:
                return 'error'
            start, end = parts
            regions.append((int(start) - 1, int(end)))
        return regions
    except (ValueError, IndexError):
        return 'error'
