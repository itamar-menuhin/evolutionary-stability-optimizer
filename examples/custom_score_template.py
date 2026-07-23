"""Template for scoring your own DNA sequences during optimization, instead
of the built-in codon-usage (CAI/tAI) scoring.

HOW TO USE THIS FILE
--------------------
1. Copy this file and rename it (e.g. `my_score.py`).
2. Edit the `score` function below to compute your own number for a piece
   of DNA. Bigger number = better/more desirable sequence. That's it -
   you don't need to know anything about DNAChisel or how optimization
   works internally.
3. Run:
       eso-optimize --input-folder my_sequences --custom-score-file my_score.py
   instead of passing --organism-name.

WINDOW: "per-codon" vs "whole sequence"
----------------------------------------
Leave WINDOW = 3 below if your score only needs to look at 3 letters
(one codon) at a time to decide how good that piece is - this is how the
built-in CAI/tAI scoring works.

If your score genuinely needs to look at the WHOLE gene at once to make
sense (e.g. it depends on some property of the entire sequence, not just
one codon), set WINDOW = None instead.

This is a CORRECTNESS choice, not a dependable speed choice - don't assume
WINDOW = 3 is "the fast one." Confirmed directly: for a cheap score
function, both settings perform about the same. But if your score function
does real work per call (e.g. it calls out to a model), windowed mode can
end up calling it far more times overall than whole-sequence mode, making
it meaningfully SLOWER in practice, not faster. If you're choosing WINDOW to
make things faster, benchmark both settings on your actual sequence first -
don't assume.

The only thing that matters for choosing WINDOW is whether your score is
CORRECT when computed this way: if you can imagine writing your score by
looking at 3 letters at a time and adding up the results, and that sum truly
equals what you'd compute by looking at the whole gene at once, use
WINDOW = 3. If you can't, use WINDOW = None - it's always correct.

Your score only ever sees the coding region(s) (ORF) being optimized, never
any flanking sequence outside them (e.g. a UTR) - same as the built-in
CAI/tAI scoring. If WINDOW doesn't evenly divide the length of what's being
scored, the leftover few letters at the end are silently skipped (never
passed to this function at all) - ESO will warn you if this happens, though
it can't happen for the default WINDOW = 3 case, since a coding region's
length is always a multiple of 3 already.

There's no automatic way for ESO to check whether your score genuinely
decomposes this way - it's on you to know. If you're not sure, WINDOW = None
is always correct (just slower); only reach for a WINDOW value once you're
confident your score really is a sum of independent per-piece contributions.
"""

WINDOW = 3  # 3 = score is computed per-codon and added up.
            # None = score is computed on the whole sequence at once.
            # (This is about correctness, not speed - see above.)


def score(seq):
    """Return a number for this piece of DNA - higher is better.

    `seq` is a plain text string of the letters A, C, G, T (e.g. "ATG").
    If WINDOW = 3 above, `seq` will always be exactly one codon (3 letters).
    If WINDOW = None, `seq` will be the entire gene, every time.

    EXAMPLE BELOW: this one just rewards sequences with more G and C
    letters in them. Replace it with your own logic.
    """
    return sum(1 for letter in seq if letter in "GC")
