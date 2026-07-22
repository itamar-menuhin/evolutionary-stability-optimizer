# How this repo's behavior compares to the paper

Compiled by an AI assistant working on this repo during its consolidation and bug-fixing
pass - **not a formal line-by-line diff against the paper's actual text**, since the
assistant doesn't have access to the paper itself, only its citation. Treat this as a
starting point to verify, not a finished audit - you (as an author) are the one who can
actually check these against what the paper states.

The purpose of this doc: if someone reads the paper and then uses this tool, a few things
could look like discrepancies. This sorts them into what's actually a change vs. what
only looks like one.

## Confirmed NOT a divergence (despite looking like one)

Two things changed during this repo's development that could raise a flag if compared
naively against the paper, but turned out - on closer investigation - to match it, or to
provably not change any result:

- **The recombination-risk constant.** A cross-check against the EFM Calculator tool's
  own code initially suggested this codebase's `a=5.8` should be `a=8.8` instead. Reading
  the actual primary source (Oliveira et al. 2008, Plasmid 60:159-165) directly showed the
  EFM Calculator's own code conflates two different rows of that paper's Table 3 - `5.8`
  is genuinely the correct value from the recA+ row this codebase's other constants
  (`b=1465.6`, `alpha=29`) already match. The value was briefly changed to `8.8`, then
  reverted after this was caught. Current behavior: unchanged from before this whole
  detour, and matching what Oliveira et al. 2008 actually specifies.
- **The homopolymer (single-nucleotide repeat) detection threshold**, raised from 4 to 12
  during a performance pass. This sounds like it could change which sites get reported,
  but was proven output-neutral: every homopolymer run of 6-11 nt is *already* detected
  and reported via a different code path (as a 2-nucleotide-repeat site), and the -9
  log10-probability risk filter already rejects every run shorter than 6 outright
  regardless of the seed. The threshold change only skips generating candidates that were
  always going to be discarded anyway - confirmed via a 300-trial randomized fuzz test
  (0 mismatches) plus a targeted sweep of runs length 4-20 (0 mismatches). Current
  behavior: identical final detected hotspots to before this change.

## Capabilities added beyond what a 2022 paper would describe

These are real additions to the tool, not something the paper's Methods section would
cover, since they didn't exist yet:

- **Custom scoring** (`--custom-score-file` / `custom_score_fn=`) - scoring sequences by
  an arbitrary user-supplied function instead of CAI/tAI codon-usage tables.
- **A second, independently-developed implementation** of the recombination and slippage
  detectors (`eso.detection.staubility_variant`, selectable via `--recombination-mode
  fast` / `--slippage-mode fast`), alongside the original.
- **Additional motif types for `--compute-motifs`** beyond methylation: a cryptic
  ribosome-binding-site scanner (Shine-Dalgarno) and cryptic bacterial promoter element
  scanner (sigma70 -35/-10 hexamers), plus a way to define your own motif from a plain
  IUPAC consensus string without a MEME/PSSM file.
- **`--indexes-file`** - a CLI-level way to specify per-sequence ORF/exclusion regions
  from a file, rather than only programmatically.

If asked "is this in the paper," the honest answer for any of the above is no.

## One thing that needs your direct judgment, not a GPT's

**Methylation-motif avoidance was found to be a complete no-op** in this codebase - an
off-by-one in how a detected motif's location was converted into a DNAChisel constraint
meant the constraint could structurally never detect the motif was still present, so
`optimize()` never actually removed any methylation motif it was asked to avoid, silently
(confirmed directly: before the fix, a real motif survived optimization completely
untouched - 0 edits at its exact position; after, it's genuinely removed). This is now
fixed.

**What's genuinely unclear** (and only you can resolve, since it depends on the history
of this exact code vs. whatever produced the paper's results): was this bug present in
whatever version of the tool actually generated the paper's reported results, or was it
introduced later - e.g. during this repo's consolidation/port from an earlier codebase?
If any published figure, table, or claim in the paper depended on methylation-motif
avoidance actually working, and the underlying code at the time had this same bug, that
specific result may not be reproducible with the code as it existed then (only with the
current, fixed version). If the bug was introduced during this repo's later
consolidation instead, the paper's original results are unaffected.

This is exactly the kind of question a GPT (or anyone without access to the original
pre-consolidation codebase and exact result-generation scripts) can't answer - if your
colleague or anyone else asks about this, the honest response is "that needs checking
against the original results directly," not a guess.
