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
built-in CAI/tAI scoring works, and it's FAST, because the optimizer only
has to re-check the codons it actually changed.

If your score genuinely needs to look at the WHOLE gene at once to make
sense (e.g. it depends on some property of the entire sequence, not just
one codon), set WINDOW = None instead. This always works, but is much
SLOWER, because the optimizer has to recompute your score for the entire
sequence every single time it tries a change, however small.

If you're not sure which applies to you: if you can imagine writing your
score by looking at 3 letters at a time and adding up the results, use
WINDOW = 3. If you can't, use WINDOW = None.
"""

WINDOW = 3  # 3 = score is computed per-codon and added up (fast).
            # None = score is computed on the whole sequence at once (slower).


def score(seq):
    """Return a number for this piece of DNA - higher is better.

    `seq` is a plain text string of the letters A, C, G, T (e.g. "ATG").
    If WINDOW = 3 above, `seq` will always be exactly one codon (3 letters).
    If WINDOW = None, `seq` will be the entire gene, every time.

    EXAMPLE BELOW: this one just rewards sequences with more G and C
    letters in them. Replace it with your own logic.
    """
    return sum(1 for letter in seq if letter in "GC")
