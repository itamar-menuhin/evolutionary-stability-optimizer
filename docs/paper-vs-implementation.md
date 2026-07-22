# How this repo's behavior compares to the paper

Checked directly against the paper's text (Menuhin-Gruman et al., 2022, *ACS Synthetic
Biology* 11(3):1142-1151) - quotes below are taken from it directly, not reconstructed
from memory. This replaces an earlier version of this doc that was written before the
paper's text was available and got two things imprecise as a result (noted below).

## Confirmed NOT a divergence

Two things changed during this repo's development that could raise a flag if compared
naively against the paper, but are confirmed - directly against the paper's own stated
values - to match it, or to provably not change any result:

- **The recombination-risk constant.** The paper's Methods section states explicitly:
  "*A = 5.8 ± 0.4, B = 1465.6 ± 50.0, and α = 29.0 ± 0.1 were found empirically*" (citing
  Oliveira et al. 2008). This codebase's `calc_recombination_score` uses exactly
  `a, b, c, alpha = 5.8, 1465.6, 0, 29` - a direct match. (A brief attempt to "fix" this to
  `8.8`, based on a different tool's code rather than the paper/primary source itself, was
  caught and reverted before ever being part of a release.)
- **The homopolymer (single-nucleotide repeat) detection threshold.** The paper states
  SSR sites are considered at "*(N ≥ 3, L ≥ 2)... or (N ≥ 4, L = 1)*" - i.e. a run of 4+
  identical nucleotides is the threshold for reporting. This codebase's *internal
  candidate-generation* seed was raised from 4 to 12 purely as a performance
  optimization - but this doesn't change what gets reported at N≥4: any run of 4-11
  identical nucleotides that would be reported is generated and scored via a different,
  equivalent code path (as a 2-nucleotide-repeat site) either way. Verified output-neutral
  via a 300-trial randomized fuzz test (0 mismatches) and a targeted sweep of runs length
  4-20 (0 mismatches) - the paper's stated N≥4 threshold is still exactly what the tool
  reports today.

## Confirmed genuine divergences

Two real differences between what the paper describes and what the current defaults do:

- **Default number of sites considered.** The paper states: "*For computational
  considerations, we offer avoidance of the 10 most probable sites from each type (SSR,
  RMD, and methylation or custom motif).*" The current codebase's default is
  `num_sites=inf` (unbounded) in both `eso.pipeline.main()` and the `--num-sites` CLI
  flag - every detected site above the risk-score cutoff is considered, not just the top
  10 per category. **If you're trying to reproduce behavior described in the paper, pass
  `num_sites=10` (or `--num-sites 10`) explicitly** - the current default does not match
  what the paper describes as the tool's behavior.
- **The current *default* recombination-detection algorithm may not be the one the paper
  describes.** The paper's pseudocode (step 3, Methods) describes finding recombination
  sites by: "*Divide the sequence into subsequences of length 16, find those appearing
  more than once, and merge together if they are subsequent*" - exact-match merging, with
  no mention anywhere in the Methods section of tolerance for mismatches. The current
  **default** (`--recombination-mode thorough` / `eso.detection.recombination`) instead
  requires only Levenshtein distance ≤1 between two sites, not an exact match - a more
  permissive, different algorithm. The **non-default** `fast` mode
  (`eso.detection.staubility_variant` - note the name) does plain exact-match merging,
  matching the paper's pseudocode far more closely. This strongly suggests `thorough`'s
  Levenshtein tolerance is a **post-publication improvement**, not what was used to
  generate the paper's reported results. **If reproducing the paper's case study or
  results specifically, `--recombination-mode fast` is the closer match to what's
  described**, not the current default.

Both of these mean: **running this tool with default settings today will not necessarily
reproduce what the paper reports** - it will generally find and avoid *more* sites (no
cap, plus a more permissive recombination-match criterion) than the tool as originally
described. This isn't a bug in either direction - the current defaults are a reasonable,
arguably better-motivated choice for actual use (why silently ignore real hotspots past
an arbitrary top-10 cutoff, why miss near-duplicates that differ by one substitution) -
but it does mean paper reproduction and "best practice for a new sequence today" are two
different configurations, and worth being explicit about which one you want.

## Capabilities added beyond what the paper describes

Checked directly against every mention of "custom" in the paper: they all refer to
custom **motifs/PSSM matrices to avoid** ("users may provide their own PSSM matrices for
sites to be avoided", "This standard format allows users to import custom motifs") - the
same concept as today's `--motifs-path`/`--common-motifs`/`motif_from_consensus`, which
already existed and was just made easier to use this session (bundled common motifs,
IUPAC-consensus shorthand). The paper does **not** describe anything like a custom
scoring *function* replacing codon-usage-based optimization - `eso.custom_score` /
`--custom-score-file` (score sequences by arbitrary Python logic instead of CAI/tAI) has
no precedent in the paper at all. If asked "is custom scoring in the paper," the honest
answer is no - the paper's "custom sites" means custom motifs to avoid, a different
feature from custom scoring.

Other things not in the paper, since they didn't exist yet:
- **A second, independently-developed implementation** of slippage detection
  (`eso.detection.staubility_variant`, also selectable via `--slippage-mode fast`) -
  unlike recombination, this one is confirmed a pure speed choice with no sensitivity
  difference from the default.
- **Additional motif types for `--compute-motifs`** beyond methylation: a cryptic
  ribosome-binding-site scanner (Shine-Dalgarno) and cryptic bacterial promoter element
  scanner (sigma70 -35/-10 hexamers).
- **`--indexes-file`** - a CLI-level way to specify per-sequence ORF/exclusion regions
  from a file, rather than only through the GUI the paper describes (this repo has no
  GUI at all - see below).
- **No GUI.** The paper describes and screenshots a downloadable graphical application
  (Figure 5). This repo is CLI/Python-library only - there is no GUI here.

## Confirmed NOT a divergence (additional checks)

- **GC-content window size**: paper states "*GC content will be enforced on each
  subsequence with length 50*"; current default `window_size_gc=50` matches.
- **Codon optimization methods**: paper lists "use best codon", "match codon usage", and
  "harmonize RCA"; current `--method` choices (`use_best_codon`, `match_codon_usage`,
  `harmonize_rca`) match exactly.
- **Recombination site minimum length**: paper states RMDs are "*long (L ≥ 16)*"; current
  detectors use the same 16nt minimum.

## One thing that needs your direct judgment, not a GPT's

**Methylation-motif avoidance was found to be a complete no-op** in this codebase before
this session's fix - an off-by-one in how a detected motif's location was converted into
a DNAChisel constraint meant the constraint could structurally never detect the motif was
still present, so `optimize()` never actually removed any methylation motif it was asked
to avoid, silently (confirmed directly: before the fix, a real motif survived
optimization completely untouched - 0 edits at its exact position; after, it's genuinely
removed).

**What's genuinely unclear** (and can't be checked against the paper's text - only
against the actual code and results that generated the paper's figures): was this bug
present in whatever version of the tool actually produced the case study in Figure 4/the
conservation-score analysis in Figures 6/8, or was it introduced later, e.g. during this
repo's consolidation from an earlier codebase? If the paper's reported case study
(BBa_I13604: "186 nucleotides were changed... removing five potential recombination
sites and 26 potential slippage sites") or the conservation-score validation depended on
methylation-motif avoidance actually functioning, and the code at the time had this same
bug, that specific result may not be reproducible with the code as it existed then (only
with the current, fixed version). If the bug was introduced during later consolidation
instead, the paper's original results are unaffected. This is exactly the kind of
question a GPT (or anyone without the original pre-consolidation codebase and exact
result-generation scripts) can't answer - if your colleague or anyone else asks about
this, the honest response is "that needs checking against the original results
directly," not a guess.
