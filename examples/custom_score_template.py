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

Your score only ever sees the coding region(s) (ORF) being optimized, never
any flanking sequence outside them (e.g. a UTR) - same as the built-in
CAI/tAI scoring. `score` is called on the whole ORF at once, on every trial
change tried during optimization - this can be slow for a long sequence or
an expensive `score` function (a warning is raised when you use this).
"""


def score(seq):
    """Return a number for this piece of DNA - higher is better.

    `seq` is a plain text string of the letters A, C, G, T (e.g. "ATG..."),
    the whole coding region being optimized.

    EXAMPLE BELOW: this one just rewards sequences with more G and C
    letters in them. Replace it with your own logic.
    """
    return sum(1 for letter in seq if letter in "GC")
