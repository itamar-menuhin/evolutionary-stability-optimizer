# Detector comparisons and mode decisions

ESO carries independently-developed implementations of several hotspot detectors
(see the main README and the repo's provenance annotations). This doc records,
one detector type at a time, what was actually tested, what bugs were found and
fixed along the way, and the resulting recommendation for which mode to use
when.

## Cross-cutting bug: `ranges_overlap` treated touching ranges as overlapping

Found while building and stress-testing a (still-experimental, not adopted -
see below) regex prefilter for slippage detection, not during the original
overlap-collapse work itself. `eso.detection._overlap.ranges_overlap` used
`a[0] <= b[1] and b[0] <= a[1]`, which is correct for *inclusive*-end ranges
but wrong here: `(start, end)` throughout this codebase uses Python-slicing
*exclusive*-end semantics (`seq[start:end]`), so a range ending at 950 and one
starting at 950 share no actual character - they merely touch. The `<=`
version counted that as overlapping, so `collapse_overlapping_intervals`
(used by both recombination and slippage) would silently discard a
genuinely distinct, adjacent, lower-scoring hotspot whenever it happened to
sit immediately next to a higher-scoring one - e.g. `"A"*950 + "T"*300 +
"A"*950` reported only the two flanking A-runs and dropped the T-run
entirely. Fixed to strict inequality (`a[0] < b[1] and b[0] < a[1]`).

This was a real correctness regression in already-shipped code (introduced
when `collapse_overlapping_intervals` was built for the redundancy fixes
above), not something the regex prefilter work introduced - it was already
live before this was found. Verified via `tests/test_detection_overlap.py`
and a regression test in `tests/test_detection_slippage.py`, plus the full
300-trial fuzz sweep and test suite. One narrow, pre-existing, now-*exposed*
difference surfaced as a side effect: `eso.detection.slippage` and
`eso.detection.staubility_variant` can independently pick a different
1-position phase for a length>1 periodic region that sits immediately next
to another repeat (e.g. reporting `"TCTCTC"` vs `"CTCTCT"` for the same
physical alternating run) - both readings are individually valid, and the
buggy `<=` check had been masking this by discarding both in the same
direction. Not a correctness bug, just a documented quirk of two
independently-developed algorithms; observed in 1/300 fuzz trials.

## Recombination detection

### Cross-validated against the authoritative reference implementation - initially "fixed" a constant, then had to revert that fix after reading the primary source

Prompted by a direct ask to look for libraries or reference sources that
would *increase trust*, not just speed - the first trust signal checked was
the EFM Calculator tool itself
([github.com/barricklab/efm-calculator](https://github.com/barricklab/efm-calculator)),
whose `get_recombo_rate()` implements the same formula shape
`calc_recombination_score` does. Its hardcoded constant is `a=8.8`, vs. our
`a=5.8` - which was initially changed to match (see git history), based on
that tool's own docstring attribution to "Oliveira et al." But that
attribution was itself second-hand - not yet checked against the actual
paper - and a direct challenge to go verify against the primary source
before trusting the "fix" turned out to be warranted.

**Reading Oliveira et al. 2008** (Plasmid 60:159-165,
doi:10.1016/j.plasmid.2008.06.004) directly, Table 3 gives two full
parameter sets for Eq. (4), `FR(LR,LS) = (A+LS)^(-a/LR) * LR/(1+B*LR+C*LS)`:

| | A | B | C | a (exponent) |
|---|---|---|---|---|
| recA⁻ strains | 200.4 | 2163.0 | 14438.6 | 8.8 |
| recA⁺ strains | 5.8 | 1465.6 | - | 29.0 |

Our `b=1465.6` and `alpha=29` match the recA⁺ row exactly - but that row's
`A` constant is 5.8, not 8.8. The value 8.8 only appears in this table as
the recA⁻ row's *exponent*, a different parameter of a different model.
The EFM Calculator's own `a=8.8` therefore combines the recA⁺ row's B and
exponent with the recA⁻ row's exponent mislabeled as `A` - a transcription
mix-up in the reference tool between two rows of the same source table, not
a correction to this codebase's original value. **`a=5.8` (the original
ESO_curr/STABLES value) is what the primary source actually supports; the
change to 8.8 was a mistake and has been reverted.**

This is the standing lesson from this whole detour: cross-checking against a
tool's *code* found a real-looking discrepancy, but only reading the actual
*paper* the code claims to implement caught that the code itself doesn't
match its own cited source. "Verified against the reference implementation"
and "verified against the primary literature" are not the same claim, and
shouldn't be reported interchangeably.

Reverted in both `eso.detection.recombination.calc_recombination_score` and
`eso.detection.staubility_variant.find_recombination_sites` (the formula is
duplicated between them), with docstrings now citing Oliveira et al. 2008
directly instead of the EFM Calculator tool. `calc_recombination_score`'s
reference test now pins to `a=5.8`
(`test_calc_recombination_score_matches_efm_calculator_reference`), and a
100-trial fuzz sweep confirms this only affects the risk score, not the
underlying near-duplicate *detection* ability (still 0 misses).

**Two other library candidates were investigated for this same "increase
trust" goal and ruled out** (not adopted, no code changes):
- **`rapidfuzz`**'s bounded/early-exit batch distance search, as a possible
  correctness cross-check or replacement for the neighbor-generation
  matching in `thorough` - algorithmically worse at scale (still O(n^2) all-pairs
  under the hood) despite a faster per-comparison constant; see the
  "Implementation-level optimization" section below for the benchmark.
- **A regex backreference-based rewrite of slippage detection** - a real
  ~20-44x raw speedup on the candidate-scan sub-step, but a genuine
  correctness gap (nested short-period-inside-a-longer-repeat cases missed
  in ~9.9% of randomized trials) that would need real design work to close;
  a follow-up chunked-prefilter design was fully validated for correctness
  (0 mismatches across ~1,100 fuzz + boundary-stress trials) but provided
  no real-world speedup, since realistic/random DNA almost always contains
  *some* incidental short repeat in any multi-kb window, so the "skip empty
  regions" strategy essentially never fires. Not adopted either way.

Two implementations, both scoring with the identical empirical EFM Calculator
formula, differing only in how they decide two sites are "the same":

- **`thorough`** (`eso.detection.recombination`, from `ESO_curr`) - Levenshtein-tolerant:
  generates every single-edit (insertion/deletion/substitution) neighbor of
  each candidate 16-17mer and checks whether any neighbor occurs elsewhere in
  the sequence. Catches near-duplicates, not just exact ones.
- **`fast`** (`eso.detection.staubility_variant`, from STABLES) - exact-match only:
  counts every 16-mer via a vectorized n-gram count and flags any that occur
  more than once. Cannot tolerate even a single mismatch between two otherwise
  identical sites.

Both are exposed through `eso.detection.dispatch.find_recombination_sites(seq,
num_sites, mode="thorough" | "fast")`, and threaded through
`eso.pipeline.main(..., recombination_mode=...)` / `eso-optimize
--recombination-mode {thorough,fast}`.

### What was actually tested

Ad hoc scripts (not fixtures kept in the repo - the resulting behavior is
captured in `tests/test_detection_dispatch.py` and
`tests/test_detection_recombination.py` instead) covered:

1. **Exact duplicate** - both modes catch it.
2. **Near-duplicate, mutation near an edge** - both modes catch it. This was a
   near-miss finding: `fast` isn't actually blind to single-nucleotide changes
   in general, it just needs the mismatch to leave a 16+nt exact stretch
   somewhere in the site. A mutation 17 characters into a 20nt site still
   leaves a 16nt exact prefix, which `fast` happily matches.
3. **Near-duplicate, mutation dead-center** (prefix and suffix both <16nt) -
   `thorough` still catches it (1 distinct pair); `fast` returns nothing. This
   is the actual, narrower condition under which the two modes disagree.
   Codified as `test_thorough_catches_centered_near_duplicate_that_fast_misses`.
4. **No repeats** - both correctly return empty.
5. **Runtime vs. sequence length**, 300nt-6,500nt, one embedded repeat per
   length (see benchmark below).

### Bugs found and fixed during testing

- **`staubility_variant.find_recombination_sites` crashed on any repeat longer
  than 16nt.** Root cause: `df.loc[:, 'end'] = values` does not change an
  existing column's dtype in pandas - after a `None`-then-backfill step left
  `'end'` as `float64`, the subsequent `.astype(int)` was silently discarded on
  reassignment, and slicing the sequence with a float index raised
  `TypeError`. Traced to the same pattern in the original STABLES source, so
  this predates the port; it just happens to only trigger when a repeat
  exceeds the 16nt window (the common, meaningful case). Fixed by reassigning
  the whole column (`df['end'] = ...`) instead of `.loc`-assigning into it.
- **`recombination.find_recombination_sites` returned many duplicate/near-duplicate
  rows per real hotspot** - up to 1,023 rows for a single real ~150nt repeat in
  one benchmark case. Two causes: (a) many different candidate seed 16-mers
  converge, via `elongate_sites`, on the exact same final site pair, with no
  deduplication before returning; (b) even after removing exact duplicates,
  different seeds can converge on *slightly different* boundary extents for
  the same real hotspot. Fixed with `_collapse_overlapping_pairs`: sort by
  score descending, keep a pair only if neither its site_1 nor site_2 range
  overlaps an already-kept pair's corresponding range (classic non-max
  suppression). This also removed a related bug where the old `num_sites`
  capping logic (`drop_duplicates(...).merge(...)` back onto the
  undeduplicated frame) silently reintroduced the very duplicates it was
  trying to cap.

Neither bug affected the *optimizer's* correctness (downstream constraint
conversion already deduplicated on exact coordinates before building DNAChisel
constraints) - but both made the raw `recombination_sites.csv` a human would
open unreadable, and the `staubility_variant` crash blocked using it at all
for realistic repeat sizes.

### Benchmark (superseded - see "Implementation-level optimizations" below)

One embedded repeat per sequence length, wall-clock time to detect it, **as
originally ported, before the pandas-overhead fix documented below**:

| Length (nt) | `thorough` | `fast` | Speedup |
|---|---|---|---|
| 300 | 1.481s | 0.028s | 54x |
| 500 | 2.383s | 0.025s | 95x |
| 900 | 3.844s | 0.021s | 180x |
| 1,700 | 16.277s | 0.051s | 323x |
| 3,300 | 19.134s | 0.024s | 791x |
| 6,500 | 76.354s | 0.040s | 1,902x |

`thorough` grew roughly quadratically with sequence length; `fast` stayed
flat. **These numbers are now stale** - `thorough` is 5-78x faster than this
table shows. See below for current numbers and updated recommendation (and
further below, for numbers up to 1,000,000nt showing `thorough` never
actually becomes impractical in the first place).

### Recommendation (updated after the optimization below, and again after the 1Mb re-test below)

- **`thorough` (default)** for essentially any realistic sequence length -
  gene-length through whole-plasmid scale and beyond, up to at least
  1,000,000nt where it still finishes in ~2 minutes (see the 1Mb benchmark
  further below). `fast`'s sensitivity gap (missing a centrally-mutated
  near-duplicate) is a real cost, and no length was found in testing where
  `thorough`'s own runtime forces the tradeoff - so there's no length-driven
  reason to switch by default.
- **`fast`** when the ~19-34x speed gap actually matters for a given
  workload - many-sequence batch processing where every second compounds, or
  sequences well beyond the 1Mb tested here - not because `thorough` becomes
  impractical at any tested length, but because `fast` is simply faster.
- Kept as two explicit modes rather than merged into one canonical
  implementation, per project decision (2026-07) - revisit if/when the
  slippage and methylation detector pairs are worked through and a pattern for
  unification (or not) emerges.

### Implementation-level optimization: stop using pandas inside the hot loop

Prompted by a direct question - *this code was hand-written; are there easy
wins, standard libraries doing the same work?* - profiled `thorough` with
`cProfile` on an ~830nt sequence instead of guessing. Result: `_generate_relevant_pairs`
was 5.1s of 7.0s total, almost entirely pandas call overhead (`.isin()`,
`_ixs`, `_get_value`, `__getattr__`) from filtering `df_substitutions` /
`df_deletions` / `df_insertions` with a fresh `df.sequence.isin(neighbor_set)`
call **once per candidate site** (up to O(n) sites) - the classic
"pandas-in-a-Python-loop" anti-pattern, not an algorithmic problem.

Fixed by building a plain `dict` mapping `sequence -> [(start, end), ...]` for
each of the three candidate frames *once*, and iterating the input frame via
`zip()` instead of `.iloc[ii, col]` - same algorithm, same pairs produced,
just without pandas overhead reintroduced on every iteration. Verified
equivalent to the original (kept temporarily, not committed) implementation
via a 100-trial fuzz sweep comparing exact `(start_1,end_1,start_2,end_2)`
sets - 0/100 mismatches - since dict iteration order differs from the
original DataFrame row order, but every downstream step (collapse, sort) is
independent of input order.

No standard library swap was actually needed here - the fix is entirely
"stop paying pandas' per-call overhead in a loop," not a smarter algorithm.
(`rapidfuzz`, already an indirect dependency via `python-Levenshtein`, was
considered for the underlying edit-distance work, but the profiling showed
the cost was in the *pandas plumbing*, not the Levenshtein distance
computation itself, so replacing that wouldn't have addressed the actual
bottleneck.)

**Updated benchmark**, same methodology as the stale table above:

| Length (nt) | old | new | speedup |
|---|---|---|---|
| 300 | 1.481s | 0.291s | 5.1x |
| 500 | 2.383s | 0.409s | 5.8x |
| 900 | 3.844s | 0.422s | 9.1x |
| 1,700 | 16.277s | 0.480s | 33.9x |
| 3,300 | 19.134s | 0.645s | 29.7x |
| 6,500 | 76.354s | 0.973s | 78.5x |

And pushed further, since it's now fast enough to test at scale:

| Length (nt) | `thorough` | `fast` |
|---|---|---|
| 12,900 | 2.191s | 0.066s |
| 25,700 | 3.433s | 0.161s |
| 51,300 | 6.527s | 0.229s |

Growth is now roughly linear-to-mildly-superlinear rather than quadratic -
51,300nt in 6.5s, where the old implementation's quadratic trend would have
put it at tens of minutes. This is why the recommendation above changed:
`thorough` stopped being the kind of implementation where reaching for
`fast` at moderate lengths was ever necessary. (At the time this was
written, the plan was to reach for `fast` "in the tens of kb" once
`thorough` presumably got inconvenient beyond that - but that was an
extrapolation, not a measurement; see the 1,000,000nt re-test further below,
which found no such point up to 1Mb.)

**Re-verified after the `ranges_overlap` fix and the revert back to `a=5.8`**,
since neither is a no-op on timing even though neither changes algorithmic
complexity (random ACGT sequences, `random.seed(0)`, single run per length,
same methodology as the tables above):

| Length (nt) | `thorough` (a=5.8, current) | `fast` (a=5.8, current) |
|---|---|---|
| 400 | 0.034s | 0.001s |
| 3,400 | 0.275s | 0.005s |
| 6,600 | 0.447s | 0.012s |
| 13,000 | 1.027s | 0.023s |
| 25,800 | 2.025s | 0.046s |
| 51,400 | 4.256s | 0.124s |

These are, if anything, a little faster than the transient `a=8.8` table
above rather than slower - despite `a=5.8` giving systematically
*higher* (less-negative) scores and thus letting more marginal candidates
past the `-9` filter, as expected. With only one random sequence sampled per
length (no injected repeat, so which incidental candidates show up is down
to chance), run-to-run variance from the actual sequence content dominates
over the constant's effect at this sample size - a second run at these same
lengths landed within the same range (e.g. 51,400nt: 4.256s vs. 4.324s), so
the numbers above are stable, but they shouldn't be read as isolating the
`a=5.8` vs. `a=8.8` effect on its own.

**Extended to 1,000,000nt (2026-07), because the "crossover to `fast` sits
somewhere in the tens of kb" line above was never actually tested past
51,400nt** - it was a plausible-sounding extrapolation, not a measurement. It
turned out not to hold up: no crossover to "impractical" was found anywhere
in the tested range. Same methodology as the table above (random ACGT,
`random.seed(0)`, single run per length, no injected repeat):

| Length (nt) | `thorough` | `fast` | thorough/fast |
|---|---|---|---|
| 100,000 | 10.633s | 0.457s | 23x |
| 200,000 | 20.407s | 0.889s | 23x |
| 400,000 | 45.597s | 2.079s | 22x |
| 700,000 | 86.175s | 4.172s | 21x |
| 1,000,000 | 126.651s | 6.815s | 19x |

Growth from 51,400nt through 1,000,000nt stays close to linear throughout -
local exponent (`log(time ratio)/log(size ratio)` between each successive
pair of points, including the 51,400nt point above as the starting anchor)
comes out to 1.38, 0.94, 1.16, 1.14, 1.08 across the five doubling-ish steps
above. That's noisy (single run per point, no injected repeat, so which
incidental near-duplicate candidates show up - and how many pairwise
comparisons they trigger - varies with the specific random sequence, same
caveat as the table above), but there's no trend of the exponent climbing
as length grows, which is what true superlinear/quadratic blowup would look
like. **`thorough` never became impractical at any length tested up to
1,000,000nt** - the worst case measured was 126.7s (~2.1 minutes) at 1Mb,
nowhere close to a 5+ minute threshold, let alone the tens-of-minutes the
original pre-optimization quadratic implementation would have produced at
far shorter lengths.

**This means the "practical crossover to `fast` in the tens of kb" framing
from the sections above was an overstatement and should be walked back.**
No actual crossover point exists in the range tested (300nt to 1Mb) -
`thorough` stays comfortably fast the entire way, growing roughly linearly
rather than hitting a wall. `fast` is still meaningfully faster in absolute
terms (~19-23x at these larger lengths, similar to the ~34x gap already
noted at 51,400nt), so it remains worth reaching for when every second
compounds - many-sequence batch workloads, or truly whole-genome-scale
inputs well beyond 1Mb that weren't tested here. But `fast` is not required
for tractability at any realistic construct size (whole plasmids included):
it's a speed optimization on top of an already-practical `thorough`, not a
rescue from an impractical one. The recommendation above ("`thorough`
default, `fast` once its runtime actually matters") still holds, but the
justification changes - not "`thorough` becomes inconvenient out in the
tens of kb," but "`thorough` stays inconvenient-free through at least
1,000,000nt, so `fast` is an optional speedup, not a requirement, unless a
workload's realistic scale or volume makes even ~2 minutes per sequence add
up."

## Slippage detection

Two implementations of the same SSR (short tandem repeat) concept:

- **`eso.detection.slippage`** (`ESO_curr`) - for each base-unit length 1-15,
  pre-filters to candidate repeat units that already occur 3+ times back-to-back
  somewhere in the sequence (`seq.find(subunit*3) > -1`), then finds every
  non-overlapping occurrence of each candidate.
- **`eso.detection.staubility_variant`** (STABLES) - for each length and every
  frameshift offset within that length, splits the sequence into fixed-size
  chunks and scans for 3 (or 6, for length 1) identical chunks in a row.

Both exposed through `eso.detection.dispatch.find_slippage_sites(seq,
num_sites, mode="default" | "fast")`, and threaded through
`eso.pipeline.main(..., slippage_mode=...)` / `eso-optimize --slippage-mode
{default,fast}`.

### What was actually tested

1. The three cases already covered by `tests/test_detection_slippage.py`
   (dinucleotide repeat, single-nucleotide run, no repeats), run through both
   implementations side by side.
2. A 300-trial randomized fuzz sweep: random prefix + a repeated unit (length
   1-8, count 3-7) + random suffix, checking whether both implementations
   agree on (a) whether the true repeat region is detected at all, and (b) the
   resulting row count. **0/300 mismatches** on both counts after the fixes
   below - unlike recombination, these two implementations are fully
   equivalent in what they detect.
3. Runtime vs. sequence length, 315nt-9,615nt (see benchmark below).

### Bugs found and fixed during testing

Four, all in `eso.detection.staubility_variant` (its `eso.detection.slippage`
sibling only needed the redundancy fix):

- **Crashed on any repeat long enough to trigger the back-to-back merge
  logic** - the identical `df.loc[:, 'end'] = ...` dtype bug found in the
  recombination variant (see above), present in the same shape here.
- **Crashed with `IndexError` whenever a length>1 repeat sat at the boundary
  of its scan range** - `is_followed1`'s `curr_seq_split[ii + 3]` access ran
  unconditionally whenever the first two equality checks passed, *regardless
  of `length`*, because the `and length == 1` guard sat last in the boolean
  chain rather than first - so it didn't short-circuit before the
  out-of-bounds access for `length > 1`. Confirmed present in the original
  STABLES source (`Staubility_Code_shimshi.py` lines 1081-1084) via a
  side-by-side read, so this predates the port; found only by randomized fuzz
  testing (`CACGCATTTCCCCCCTACATCACCAGAGAG`), since it requires a genuine
  repeat to sit exactly at a length/frameshift boundary. Fixed by reordering
  the boolean chain so `length == 1` / `length > 1` is checked *first*.
- **Missing the `-9` risk-score cutoff entirely** - a 4-5nt homopolymer run
  (log10 prob ~-9.3 to -9.98, below the cutoff both `eso.detection.slippage`
  and this module's own `find_recombination_sites` apply) was reported
  unfiltered. Traced to the original STABLES codebase: the filter genuinely
  isn't inside `find_slippage_sites` there either - it's applied by the
  *caller*, elsewhere in `Staubility_Code_shimshi.py`
  (`df_slip.loc[df_slip.log10_prob_slippage_ecoli > -9]  # only 'severe'
  constraints`). Extracting `find_slippage_sites` in isolation lost that
  downstream step. Restored inside the function so it's self-contained and
  consistent with its sibling.
- **No deduplication at all** - unlike `eso.detection.slippage`'s
  `drop_duplicates(subset="start")`, this had no collapsing step, so
  same-start competing base-unit-length representations were both reported
  (see next bug for why "same start" dedup isn't sufficient either).

And one in `eso.detection.slippage` itself:

- **`drop_duplicates(subset="start")` missed phase-shifted overlapping
  detections.** A run like `"GCGCGCGC"` is a valid length-2 site starting at
  position N (`"GC"` x4), *and* its 1-shifted reading `"CGCGCG"` is a separate
  valid length-2 site starting at N+1 (`"CG"` x3) - different `start`, same
  physical repeat, so exact-start dedup didn't catch it.

Fixed both with a shared `eso.detection._overlap.collapse_overlapping_intervals`
helper (non-max suppression by score, generalizing the pair-based collapse
already built for recombination): sort by score descending, keep a candidate
only if its range doesn't overlap an already-kept one's.

### Benchmark (superseded - see "Implementation-level optimization" below)

One embedded repeat per sequence length, wall-clock time to detect it, **as
originally ported**, before the O(n^2)->O(n) fix documented below (row
counts matched exactly at every length: 3/3, 3/3, 5/5, 8/8, 18/18, 36/36):

| Length (nt) | `default` | `fast` | Faster mode |
|---|---|---|---|
| 315 | 0.026s | 0.055s | default (2.1x) |
| 615 | 0.042s | 0.076s | default (1.8x) |
| 1,215 | 0.092s | 0.112s | default (1.2x) |
| 2,415 | 0.266s | 0.188s | fast (1.4x) |
| 4,815 | 0.916s | 0.438s | fast (2.1x) |
| 9,615 | 3.495s | 0.663s | fast (5.3x) |

**This crossover no longer exists** - `default`'s per-candidate substring
search had an O(n^2) algorithmic problem, not just a slower constant factor;
fixing it (below) makes `default` faster than `fast` at every length tested,
including well past where this table shows `fast` "winning".

### Recommendation (updated after the optimization below)

- **`default`** essentially always, per the updated benchmark below - there
  is no longer a length range where `fast` is actually faster.
- Unlike recombination, this remains a **pure speed choice with no
  sensitivity tradeoff** (verified via the 300-trial fuzz sweep) - so
  switching modes never changes which hotspots get reported, only how long
  it takes (and now, `default` is simply the better choice on both counts).

### Skipping candidates that can never survive the risk filter or the collapse

The `-9` log10-prob cutoff is a function of `num_base_units` alone, so it's
possible to work out ahead of time exactly which candidates it will always
reject - and skip detecting them at all, rather than finding, scoring, and
then discarding them. This went through two rounds, the second one prompted
by a sharper follow-up question: *since any homopolymer run of n>=6 is also
always a valid length-2 site, doesn't length-2 detection already catch it -
so when is length-1 detection actually needed?*

**Round 1 - the filter itself.** `log10_prob = -12.9 + 0.729n` for
length_base_unit==1 crosses `-9` at **n=6** - runs of exactly 4 or 5 always
fail outright, no matter what else exists in the sequence. (`length_base_unit
> 1`'s formula, `-4.749 + 0.063n`, is already `-4.56 > -9` at the smallest
detectable `n=3` - every length>1 candidate always survives; there's nothing
to skip there.) Raised the length-1 seed from 4-in-a-row to 6-in-a-row.

**Round 2 - the collapse.** Every homopolymer run of n>=6 is *also* always
detected as a length-2 site (6+ identical characters trivially contain "XX"
repeated 3+ times), scored as `-4.749 + 0.063*floor(n/2)`. Since
`collapse_overlapping_intervals` only keeps the single highest-scoring
representation of an overlapping region, and the length-1 score grows
~0.729/nt vs. length-2's ~0.0315/nt, **length-2 always wins for n=6..11**
(e.g. n=8: length-1 -7.068 vs. length-2 -4.497) - so length-1's detection at
those lengths is real but its result is *always* discarded. The crossover is
**n=12 exactly** (length-1 -4.152 vs. length-2 -4.371) - only from there does
length-1's reading start being the one that survives, correctly reflecting a
long homopolymer run's much higher risk than length-2's linear-in-repeats
score would suggest. Raised the seed again, from 6 to 12.

Verified both rounds output-neutral: re-ran the 300-trial fuzz sweep (0/300
mismatches) plus a targeted sweep of homopolymer runs n=4..20, 10 trials each
(0 mismatches) - and confirmed directly against the actual detector (not just
the formula) that n=6..11 are still found via length-2, and n=12 is the exact
point where length-1 first wins (see
`test_homopolymer_runs_6_to_11_are_found_via_length_2_not_length_1` and
`test_homopolymer_run_of_12_is_the_first_to_win_as_length_1`).

Implemented with a safer slice-based check in `staubility_variant`
(`len(curr_seq_split[ii:ii+12]) == 12 and len(set(...)) == 1`) rather than
extending the previous chain of explicit index comparisons to 12 terms -
slicing never raises `IndexError` past the list end, so it needs no
equivalent of the bounds-safe ordering the length>1 branch still requires.

**Isolated benchmark** (length-1 detection only, on random sequences with no
injected repeats - i.e. only the naturally-occurring homopolymer runs random
sequence content produces by chance):

| Length (nt) | seed=6 | seed=12 | candidates (seed=6 \| seed=12) | speedup |
|---|---|---|---|---|
| 20,000 | 0.0003s | 0.0018s | 17 \| 0 | n/a (both trivial) |
| 100,000 | 0.0029s | 0.0015s | 72 \| 0 | 1.9x |
| 500,000 | 0.0458s | 0.0030s | 336 \| 0 | 15.3x |
| 2,000,000 | 3.820s | 0.0104s | 1,521 \| 1 | 366x |

At 2,000,000nt, seed=6 found 1,521 "candidates" - all but one of them runs of
6-11 that were always going to lose to length-2 and never appear in the
output. seed=12 correctly finds close to none (a run of 12+ identical
nucleotides by pure chance is astronomically rare - real engineered
sequences with e.g. poly-A tails are a different story) while still catching
the one that's actually there.

**Caveat as of the two rounds above (superseded immediately below)**:
length-1 detection was a tiny fraction of total `find_slippage_sites`
runtime - the dominant cost was candidate generation for length_base_unit
2-15 (`_find_relevant_subunits_len_l`'s substring generation + `.find()`
call per candidate), which had no equivalent filter-based shortcut, since
every length>1 candidate always survives the `-9` filter on its own merits.

### Implementation-level optimization: the real bottleneck was an O(n^2) algorithm

Profiling `default` with `cProfile` on a 9.6kb random sequence (prompted by
the same *"are there easy wins?"* question that led to the recombination fix
above) found `_find_relevant_subunits_len_l` responsible for 93% of total
runtime (2.591s of 2.798s), essentially all of it in 88,151 calls to
`str.find()`. The function builds every unique length-`L` substring, then
calls `seq.find(subunit*3)` once per unique candidate to check whether it
repeats 3+ times back-to-back - each call rescans the whole sequence, so with
O(n) unique candidates (typical for non-repetitive DNA) this is O(n^2)
overall, not the O(n) the surrounding code's structure suggests.

Fixed by building a single position dict (`{substring: [positions]}`) in one
O(n) pass, then checking for 3 consecutive occurrences via O(1) set
membership (`p`, `p+L`, `p+2L` all present means `subunit*3` occurs at `p`,
exactly what the old `.find()` call was testing indirectly). Verified
output-neutral via the same 300-trial fuzz sweep (0/300 mismatches).

**This changed the mode comparison completely.** Re-running the `default` vs.
`fast` benchmark from earlier:

| Length (nt) | `default` | `fast` | Faster mode |
|---|---|---|---|
| 330 | 0.025s | 0.042s | default (1.7x) |
| 1,230 | 0.039s | 0.086s | default (2.2x) |
| 9,630 | 0.155s | 0.572s | default (3.7x) |
| 20,015 | 0.358s | 1.192s | default (3.3x) |
| 100,015 | 2.120s | 9.764s | default (4.6x) |
| 300,015 | 5.721s | 58.262s | default (10.2x) |

`default` is now faster than `fast` at **every** tested length, with the gap
*widening* rather than crossing over - the opposite of the original
benchmark's shape. `fast`'s name no longer describes its speed relative to
`default`; it remains a distinct, independently-implemented algorithm (still
useful as a cross-check / second opinion, and identical in what it detects -
row counts matched exactly at every length above), but there is no longer a
performance reason to reach for it. Left the mode API and default
(`mode="default"`) as-is rather than renaming or deprecating `"fast"` -
worth revisiting if `staubility_variant`'s own `.find()`-per-candidate
pattern (used for its recombination detector too) gets the same fix.

**Re-verified after the `ranges_overlap` and recombination-formula fixes**
(neither touches slippage's code paths, so this is mainly a check that
nothing regressed, not an expected improvement):

| Length (nt) | `default` (current) | `fast` (current) |
|---|---|---|
| 345 | 0.014s | 0.032s |
| 1,245 | 0.022s | 0.075s |
| 9,645 | 0.126s | 0.434s |
| 20,030 | 0.248s | 0.847s |
| 100,030 | 1.178s | 6.023s |
| 300,030 | 3.661s | 39.066s |

Within normal run-to-run variance of the table above (as expected - neither
fix touches this code), and the conclusion is identical: `default` wins at
every length, gap widening with size (~10.7x at 300k here vs. ~10.2x
previously - same story).

---

## Methylation motif detection

### Not the same shape of split as recombination/slippage - and, unlike those two, this one ended in deletion rather than a mode choice

Recombination and slippage's two implementations are independently-arrived-at
algorithms that (after fixing bugs) detect the *same* hotspots, differing only
in speed or, for recombination, exact-vs-near-duplicate sensitivity. Motif
detection's split was different in kind: `eso.detection.methylation` (Biopython
PSSM) and `eso.detection.staubility_variant`'s `site_motif_grader`/
`calc_max_site` scored sites by two genuinely different formulas, and neither
had a top-level `find_motif_sites`-shaped wrapper or any test coverage before
this comparison - `staubility_variant`'s motif code had never been exercised
against real motif data or wired to anything.

**Update**: `staubility_variant`'s motif scorer (`site_motif_grader`,
`calc_max_site`, and the `find_motif_sites` wrapper built around them) has
since been **deleted**, along with the `mode="raw_probability"` dispatch
option, `--motif-mode` CLI flag, and `motif_mode` pipeline parameter that were
initially built around it. The findings below led to the recommendation
"use `pssm`, it's both more accurate and faster" - which meant there was no
real reason to keep the inferior mode as a selectable option at all, rather
than leaving it in as a cross-check. The rest of this section is kept as a
historical record of that investigation, not as documentation of live code -
`eso.detection.methylation` is now the only methylation motif detector.

**`eso.detection.methylation`** scores each candidate site by its Biopython
PSSM log-odds against the motif file's own stated background letter
frequencies (`log2(observed_probability / background_probability)`, summed
over the motif's positions) - the standard formulation, and background-aware.

**`eso.detection.staubility_variant`** scores raw sequence probability under
each motif's own per-position distribution
(`sum(log10(P(observed letter at position i)))`), with **no background
correction at all**.

### A confirmed off-by-one bug, found before any correctness comparison was possible

`calc_max_site` computed `end_index = start_index + num_nucleotides + 1`,
reporting a site one nucleotide longer than the actual match. Verified
directly: a clean 4nt exact match at position 4 in a hand-built sequence
returned `actual_site = "ACGTG"` (5 characters) instead of `"ACGT"`. The
score itself (computed inside `site_motif_grader`, using a correctly-sized
slice) was unaffected - this only corrupted the reported coordinates and site
string, which would have corrupted any downstream `AvoidPattern` constraint
built from it, had this code ever been wired up. Fixed to
`end_index = start_index + num_nucleotides` (exclusive, matching this
codebase's convention elsewhere). Regression test:
`test_calc_max_site_no_longer_off_by_one`.

### The real finding: these two scores are only equivalent under a uniform background

Reasoned first, then verified directly rather than trusted on reasoning
alone (per the standing lesson from the recombination-constant detour):
under a *uniform* background (0.25/letter), Biopython's log-odds score is a
purely monotonic transform of the raw log10-probability score
(`log-odds = (raw_log10_score - length*log10(0.25)) / log10(2)`), so the two
implementations rank and threshold candidates identically. This was used to
build a uniform-background-equivalent default threshold
(`min_score = length * log10(0.25)` per motif) for `staubility_variant`'s new
`find_motif_sites` wrapper, matching `eso.detection.methylation`'s implicit
`PSSM_score > 0` filter - and confirmed empirically that both flag the same
site under a uniform background
(`test_both_implementations_find_the_same_exact_match_under_uniform_background`).

Under a **realistic non-uniform background** (e.g. a GC-rich genome,
A=0.2/C=0.3/G=0.3/T=0.2), this equivalence breaks down - not just in scale,
in actual ranking. With one such background and motif, `"AATT"` and `"CCGG"`
score *identically* under raw log10 probability (both -2.743, an artifact of
this motif's PWM symmetry combined with these sequences' letter composition),
while Biopython's log-odds score clearly separates them (0.175 vs. -2.165,
over 2 apart - AATT is rare under this background, so a decent match to it is
more surprising/noteworthy than an equally "good" raw match built from
background-common C/G letters that the raw score can't distinguish from it).
Verified in `test_implementations_diverge_under_non_uniform_background`.
**`eso.detection.methylation` (`mode="pssm"`) is the principled default
whenever a motif file's background isn't uniform - not a stylistic
preference, a real accuracy difference.**

### Benchmark: no speed advantage to offset the accuracy gap, either

5 motifs (lengths 6-15), random ACGT sequences, `random.seed(0)`, single run
per length:

| Length (nt) | `pssm` (Biopython, current) | `raw_probability` (pure-Python, current) | ratio |
|---|---|---|---|
| 400 | 0.023s | 0.053s | 2.3x |
| 3,400 | 0.148s | 0.194s | 1.3x |
| 13,000 | 0.346s | 0.646s | 1.9x |
| 51,400 | 1.324s | 2.077s | 1.6x |

`raw_probability` is consistently slower (Biopython's `.calculate()` is
C-optimized; `staubility_variant`'s is a pure-Python per-position loop),
though the gap is modest (1.3-2.3x, not an order of magnitude) - this was
`raw_probability`'s own bottleneck (or lack of one) relative to `pssm`'s
existing implementation, not the same thing as `pssm` being internally
optimal - see below, where a real, separate implementation-level bottleneck
*was* found and fixed in `pssm` itself.

### Implementation-level optimization: pssm's own bottleneck, found later

The benchmark above compares `pssm` against `raw_probability` using only a
handful of motifs - too few to expose `pssm`'s own scaling problem. With a
more realistic multi-motif set (20 motifs, lengths 4-15, matching a curated
methylation-motif database's typical size and length range), profiling
(`cProfile`) found `eso.detection.methylation.find_motif_sites` spent about
half its time in a `df.apply(axis=1)` call extracting each site's sequence -
the same pandas-per-row-overhead anti-pattern already found and fixed in
`eso.detection.recombination`/`slippage` - and most of the rest building and
converting an intermediate long-format dataframe with
`2 * num_motifs * len(seq)` rows (one row per motif, per strand, per
position), just to reduce it down to "the best-scoring match per position."

Rewrote to skip that intermediate dataframe entirely: build one
`2*num_motifs x len(seq)` numpy score matrix directly (`-inf` where a motif
doesn't reach that far), reduce it with a vectorized `np.argmax`, and extract
site strings via a plain `zip()` loop instead of `df.apply`. Verified via 500
randomized trials against the original logic - 0 genuine mismatches; the only
divergences found were cases where two *different* motifs scored *exactly*
equally at the same position, an inherent ambiguity the original code
resolved only by accident (via `pandas.sort_values`'s default,
non-deterministic quicksort) rather than by any documented rule - this
rewrite resolves those ties deterministically (lowest `motif_number` first)
instead, and both behaviors are pinned in `tests/test_detection_methylation.py`.

Benchmarked with the same 20-motif set, `random.seed(0)`, single run per length:

| Length (nt) | before | after | speedup |
|---|---|---|---|
| 3,400 | 0.254s | 0.018s | 14.1x |
| 13,000 | 0.691s | 0.050s | 13.8x |
| 51,400 | 2.850s | 0.197s | 14.5x |
| 200,000 | 11.367s | 0.788s | 14.4x |
| 1,000,000 | (not measured - old scaling makes this impractical) | 3.723s | - |

A consistent ~14x speedup, roughly flat across scales (both implementations
are linear in sequence length; this only removed constant-factor Python/pandas
overhead, not a complexity-class change) - the same qualitative story as the
recombination/slippage per-implementation fixes earlier in this doc.

### Recommendation

`eso.detection.methylation.find_motif_sites` (background-corrected PSSM) is
the only methylation motif detector in this codebase now - it was both more
accurate (background-aware) and faster than the alternative, so the
alternative was removed rather than kept as a selectable mode. If a future
need arises for a raw-probability/uniform-background variant, this section
documents the exact formula difference and the off-by-one bug to avoid
reintroducing.

## eso.constraints - the detection-to-optimization boundary

Not a detector, and not a comparison between implementations (there's only
one) - included here because this is where the "what bugs were found and
fixed" log for this codebase already lives, and this module sat completely
untested despite being load-bearing: it converts detected hotspot dataframes
into the actual DNAChisel `AvoidPattern` constraints that get optimized
against, and trims sites around user-specified exclusion (locked) regions.
A bug here silently changes what gets optimized, with no visible symptom -
exactly the kind of thing worth reviewing the same way as the detectors.

Two real bugs found, both at range boundaries under this codebase's
exclusive-end convention (matching Python slicing,
`eso.sequence_utils.parse_region`'s output, and
`eso.detection.recombination`'s output after `_elongate_sites`):

1. **`has_overlap_exclusion` treated touching ranges as overlapping** - the
   same bug class already found and fixed in
   `eso.detection._overlap.ranges_overlap`, but found independently here, in
   a separate module that fix never touched. A site ending exactly where an
   exclusion region begins (or starting exactly where one ends) shares no
   actual nucleotide with it, but the `<`/`>` comparison (instead of
   `<=`/`>=`) treated it as a conflict anyway. Fixed to `<=`/`>=`.

2. **`exclusion_site_correcter` trimmed one nucleotide more than necessary**
   on both sides of an exclusion region - `region[0] - 1` / `region[1] + 1`
   instead of `region[0]` / `region[1]`. Not a data-corruption bug (the
   trimmed `sequence` field always matched the true substring at its
   resulting `start`/`end` - verified directly, not assumed, via 2,000
   randomized trials checking `sequence == full_seq[start:end]`), but a real,
   needless loss of one still-modifiable, non-excluded nucleotide per
   exclusion boundary encountered - meaning slightly less of the sequence
   than necessary was actually available to fix a detected hotspot whenever
   it sat near a locked region.

Both confirmed as real (not hypothetical) via direct reproduction before
fixing, then verified via the same 2,000-trial randomized sweep: for every
trial, the trimmed site's `sequence` matches the true substring, no longer
overlaps any exclusion, and is *maximal* (extending it by one nucleotide on
either side, when the original untrimmed site had room to, immediately hits
either an exclusion or the original site's own boundary). See
`tests/test_constraints.py` for the full test suite (18 tests) added
alongside this - the module had none before.

The rest of the module - `_indel_recombinations`/`_substitution_recombinations`
(deciding which side of a detected recombination pair to mutate, and at what
candidate sequences, to break its Levenshtein-distance-1 relationship) and
`convert_df_to_constraints` (the final dataframe -> `AvoidPattern` list
conversion) - was reviewed and tested but no further bugs were found.

## eso.optimize - methylation-motif avoidance has never actually worked

Found while reviewing `eso.optimize.optimization_engine` (also previously
untested beyond one pipeline smoke test and the custom-score-specific tests):
**the `df_motifs` -> `AvoidPattern` constraint path was completely inert.**
Every methylation motif ESO ever detected and was asked to avoid during
optimization was silently ignored - the constraint always reported itself as
already satisfied, regardless of whether the motif was actually still present
in the sequence, so nothing ever forced it to change.

**Root cause**: `eso.detection.methylation.find_motif_sites` returns
`end_index` as *inclusive* (the index of the motif's last nucleotide -
`end_index = start_index + motif_length - 1`). Everywhere else in this
codebase - `eso.sequence_utils.parse_region`'s output, `eso.constraints`'
`has_overlap_exclusion`/`exclusion_site_correcter`, `convert_df_to_constraints`,
and DNAChisel's own `Location` - uses *exclusive* end (matching Python
slicing: `Location(3, 10).extract_sequence(seq) == seq[3:10]`, confirmed
directly). `optimize.py`'s `df_mot` construction renamed `end_index` straight
to `end` with no adjustment, so every motif's `AvoidPattern` location was
built exactly one nucleotide short of the motif's actual length.

DNAChisel's `AvoidPattern.evaluate()` searches for a pattern *only within its
given location* (`SequencePattern.find_matches`'s docstring: "Only patterns
entirely included in the segment will be returned" - confirmed by reading its
source directly). A location one nucleotide shorter than the pattern it's
supposed to contain can *structurally never* contain a match, so the
constraint always scored 0 ("Passed. Pattern not found!") no matter what was
actually in the sequence - not a subtle statistical effect, an unconditional
no-op.

**Confirmed directly, not assumed**: built a sequence with a real `GATC`
(Dam) motif site inside GC-balanced flanks (so `EnforceGCContent` wouldn't
force unrelated edits nearby that could mask the result either way), passed
it through `optimization_engine` via `df_motifs`, and checked the exact
4-nucleotide window at the motif's original position before and after.
Before the fix: 0 edits, `GATC` completely untouched. After: 1 edit, the
motif genuinely gone from that position, translation still preserved.

**Fix**: `df_mot.loc[:, 'end'] = df_mot['end'].astype(int) + 1` (was
`.astype(int)` with no adjustment) - converts the inclusive `end_index` to
the exclusive `end` everything downstream expects. This is a narrow,
targeted fix at the specific consumption site rather than changing
`eso.detection.methylation`'s own `end_index` contract, since that module is
self-consistent on its own (its `actual_site` slicing already correctly uses
`end_index + 1`, verified by its own test suite) and changing its column
semantics would be a wider, backward-incompatible change (e.g. for anyone
already reading `end_index` out of the CSV output) for a bug that's really
about how two already-self-consistent conventions get merged.

**Impact**: any past run of this pipeline with `--compute-motifs` set would
have detected methylation motifs correctly (that detector was fine) and
reported them in `motif_sites.csv` correctly, but the actual DNA sequence
optimization step would never have removed or altered any of them - the
final "optimized" sequence could still contain every methylation hotspot the
tool itself found and reported. Recombination and slippage avoidance were
not affected (`eso.detection.recombination`/`slippage`'s own `end` columns
are already exclusive - see their own docs above).

Verified via `tests/test_optimize.py` (8 tests, this module had none before),
including the regression test for this bug, one confirming exclusion regions
are still respected for motif sites too, and basic sanity checks (translation
preservation, GC-content enforcement, unknown-organism handling,
recombination/slippage avoidance mechanics).

## eso.io_utils - exclusion-region validation checked the wrong variable

Also previously untested. `test_input`'s exclusion-region validation looped
over `orf_indexes` (the ORF regions, already validated in the block just
above it) instead of the just-parsed `region_indexes` - so a malformed
exclusion region (an unparseable string, or a valid-looking one with
`start >= end`) was never actually checked. It silently passed validation
and only surfaced later as a confusing, unrelated GC-limit error message
(`"the maximal GC limit must be at least 1.0 and the minimal GC limit must
be at most 0.0"`, from `exclusion_gc_tester` being fed a nonsensical region)
instead of the clear, specific validation message that already existed and
was meant to catch exactly this. Confirmed directly:
`test_input(0.3, 0.7, {('f', '0'): ('1-9', '10-5')}, [])` returned the
GC-limit error before the fix; after, the actual
`'Start index must be smaller than end index, also for exclusion sites!'`
message. Fixed by changing the loop variable. Verified via
`tests/test_io_utils.py` (17 tests, this module had none before), covering
GC-bounds validation, ORF/exclusion region validation (including this
regression), file discovery (`relevant_file_paths`), and file
opening/gzip handling (`file_opener`).

(Noted but not chased further, given it's cosmetic and pre-existing: a
`RuntimeWarning: invalid value encountered in scalar divide` surfaces from
DNAChisel's `biotools.gc_content` when `exclusion_gc_tester` computes the GC
content of an empty window overlap (0/0) - doesn't affect the actual
validation result in the cases tested, but could be worth a guard clause if
`exclusion_gc_tester` gets revisited.)

## eso.codon_usage - the kompas table's stop-codon data was silently unusable

Also previously untested. `cub_kompas()` builds its table from 3-letter
amino acid codes (`'Ala'`, `'Arg'`, ..., and `'END'` for the stop codons),
converted to the single-letter codes DNAChisel's `CodonOptimize` expects via
`Bio.SeqUtils.seq1`. `seq1('END')` returns `'X'` (undefined amino acid), not
`'*'` (stop) - confirmed directly, not assumed. So the kompas table's
TAA/TAG/TGA frequencies were silently filed under the wrong key, and
`CodonOptimize` never had real stop-codon usage data to score against for
this host - not a crash, just quietly-absent data (`'*' not in table`,
confirmed directly). Fixed with an explicit `'END' -> '*'` special case.

While reviewing, also verified (not assumed) that all four bundled tables -
`C1`, `kompas`, `human_antibody_heavy_chain`, `human_antibody_light_chain` -
have exactly the 20 standard amino acids plus stop as keys, that every codon
actually translates (via `Bio.Seq.translate()`) to the amino acid it's filed
under, and that each amino acid's codon frequencies sum to ~1.0. No further
issues found. The two CSV-backed antibody tables don't include stop-codon
data at all (confirmed by inspecting the CSVs directly) and correctly fall
back to the documented default `{'TAA': 0.33, 'TAG': 0.33, 'TGA': 0.34}`.
Verified via `tests/test_codon_usage.py` (14 tests, this module had none
before).

## eso.report - a repeated sequence name could swallow its own page break

Also previously untested (and python-docx, needed to test it, wasn't even
installed in the dev environment - added to `pyproject.toml`'s dev
dependency group). `create_word_document_with_highlighted_differences`
decided whether to add a page break after each entry via
`if seq_name != sequences_data[-1][0]` - comparing by *name*, not position.
If an entry's name happened to match the true last entry's name (e.g. two
files, or two records, that coincidentally share a stem), it would wrongly
skip its own page break too. Confirmed directly: three entries where the
first and third share a name (`'gene'`) produced only 1 page break instead
of the expected 2. Fixed by comparing position (`index != len(sequences_data) - 1`)
instead of name. Verified via `tests/test_report.py` (6 tests, this module
had none before - gracefully skipped if the optional `python-docx` dependency
isn't installed, matching `eso.report`'s own optional-dependency design),
including direct inspection of the generated `.docx`'s actual page breaks
and per-character highlighting (not just "did it not crash").

This closes out the module-by-module review pass across the rest of `eso/`
(`constraints.py`, `optimize.py`, `io_utils.py`, `codon_usage.py`,
`report.py`) that followed the same detector comparisons above - five
previously-untested modules, five real bugs found and fixed, all documented
here alongside the detector work.

## eso.pipeline - a crash on unrecognized organisms, and a filename-stem bug

Following on directly from the module review above: `pipeline.py` and
`cli.py` hadn't been directly exercised yet. Reviewing `pipeline.py` found
two real bugs.

`_extract_cai` parsed DNAChisel's `objectives_text_summary()` output as
`text.split('\n')[0].split(':')[1].strip()` unconditionally. When
`organism_name` isn't recognized (and no `custom_score_fn` is given),
`_codon_optimization_objectives` returns `[]` and DNAChisel summarizes the
(empty) objective list as `'===> No specifications\n\n'` - no colon to split
on, so `.split(':')[1]` raised an unhandled `IndexError`. Confirmed directly
via a full `pipeline.main()` call with `organism_name='not_a_real_organism_xyz'`.
Fixed by returning `None` when there's no colon to parse, and updating
`backend()`'s report writer to say plainly that no codon-usage objective was
applied, instead of crashing.

Separately, `backend()` (and `io_utils.exclusion_gc_tester`) derived a
file's stem via `path.basename(file[0]).split('.')[0]` to look it up in the
user-supplied `indexes` dict - `.split('.')[0]` truncates at the *first*
dot, so a file named e.g. `sample.2024.fasta` silently became the stem
`sample` instead of `sample.2024`. Any `indexes` entry keyed on the intended
full stem would then simply never match, silently falling back to no
ORF/exclusion regions rather than erroring. Fixed by adding
`eso.io_utils.file_stem` (strips a trailing `.gz` if present, then takes one
`path.splitext` level) and using it in both call sites. Verified directly
against `gene.fasta`, `sample.2024.fasta`, `sample.fasta.gz`, and
`sample.2024.fasta.gz`.

## eso.optimize - a genuine DNAChisel-internal crash, investigated to its root cause

While testing `pipeline.py` against real detector output (rather than
hand-crafted dataframes), a specific but realistic scenario crashed:
resolving an `AvoidPattern` constraint whose location doesn't align to the
codon grid under `EnforceTranslation` (e.g. `modify_df_slippage`'s
per-occurrence windows for a 2nt "GC" repeat, or a 1nt homopolymer run,
starting mid-codon) could make DNAChisel's
`DnaOptimizationProblem.resolve_constraint` compute a "jointly mutable"
mutation-space region (accounting for codon interdependencies) that ends up
not overlapping the constraint's own declared location at all.
`AvoidPattern.localized()` correctly returns `None` for that (they genuinely
don't overlap) - but DNAChisel's own caller doesn't guard against that
`None` before calling `.evaluate()` on it, crashing with
`AttributeError: 'NoneType' object has no attribute 'evaluate'`. Confirmed
via monkeypatch tracing of the exact call and its `None` return, and
reproduced directly: a codon-unaligned 2nt "GC" repeat and, separately, a
1nt homopolymer run both crash; a codon-aligned 3nt repeat unit never does
(see `test_slippage_avoidance_disrupts_the_repeat`, pre-existing and never
affected).

The first fix attempted was to widen every avoid-window outward to the
nearest enclosing codon boundaries before building the `AvoidPattern`
(confirmed separately that `AvoidPattern` doesn't require its search pattern
to fill the whole `location` - only the *location* needs widening, not the
literal pattern). This did stop the crash for the "GC" repeat case
(9 edits, repeat disrupted) - but testing it against a pure `T`-homopolymer
landing on all-Phe codons (`TTT`/`TTC`, both of which always contain a T)
exposed a subtler problem: widening turns "avoid this pattern at this exact
position" into "avoid this pattern *anywhere* in the whole codon", which for
a single-nucleotide pattern in a homopolymer run is unsatisfiable for every
codon in the run (there's no Phe codon without a T), so every widened
constraint got dropped by the retry loop and the homopolymer survived
completely untouched - no crash, but also silently no effect
(`num_edits=0`), which is arguably as bad as the original crash for the
tool's actual purpose.

Given that tension, the codon-widening approach was abandoned in favor of a
narrower fix: catch this exact crash (matched by its message, so no other
`AttributeError` gets silently swallowed) in `optimization_engine`'s retry
loop, and handle it exactly like DNAChisel's own "final consistency check
failed without naming a culprit constraint" case already did - re-evaluate
every constraint and drop whichever ones are still actually failing. This
never asks DNAChisel for anything harder than what was originally detected;
it just means a site whose *only* narrow, well-targeted disruption
constraint is what triggers the crash gets dropped (silently, for now - see
below) rather than the tool crashing or falsely reporting success while
quietly doing nothing. Verified directly: the "GC" repeat case still
resolves with edits made (`num_edits=9`, deterministic across 5 repeated
runs); the `T`-homopolymer case no longer crashes and produces a valid,
translation-preserving sequence (`num_edits=0` - genuinely unsatisfiable
under this codon usage table, not a silent regression).

Investigating this crash also surfaced two independent, pre-existing bugs in
the same retry loop, both now fixed:

- `NoSolutionError.constraint` defaults to `None`, and DNAChisel's own
  `perform_final_constraints_check()` (run at the end of
  `resolve_constraints()`) raises it this way - without setting
  `.constraint` - when constraints that each individually resolved during
  the per-constraint pass end up regressing each other by the time solving
  reaches the last one. This is a real, deterministic scenario (confirmed
  via 5 repeated attempts all failing identically, not flaky), not a
  hypothetical - the retry loop's `cnst.remove(e.constraint)` crashed with
  `ValueError: list.remove(x): x not in list` whenever this happened.
- The fix for that (re-evaluating every constraint directly) crashed
  differently: `AttributeError: 'NoneType' object has no attribute 'start'`,
  because specs like `EnforceGCContent` have `location=None` until
  `.initialized_on_problem(problem)` is called (which returns a *copy*, not
  an in-place mutation) - evaluating the raw, never-initialized constraint
  crashes the same way DNAChisel's own machinery would. Fixed by evaluating
  `c.initialized_on_problem(problem, role='constraint')` instead of `c`
  directly.

Residual limitation, left as future work rather than chased further here: a
homopolymer (or other short-pattern repeat) landing entirely on
highly-constrained codons (worst case: Met/`ATG` or Trp/`TGG`, which have
only one codon each, no synonymous alternative at all) will have every one
of its per-position disruption constraints dropped as unsatisfiable, and
currently there's no visible signal to the user that this happened beyond
`num_edits` coming back lower than the number of detected slippage sites
would suggest. If this turns out to matter in practice, the fix would be to
surface which specific sites got dropped (e.g. in `final_sequence.txt` or a
returned list) rather than staying silent about it.

## eso.cli - `--num-sites` crashed whenever it actually limited anything

Last of the previously-untested modules. `--num-sites` is parsed via
argparse as `type=float` (its default, `np.inf`, isn't representable as an
`int`), and reaches `eso.detection.recombination.find_recombination_sites`
/ `eso.detection.slippage.find_slippage_sites` unchanged. Both called
`df.head(num_sites)` directly whenever a finite limit was given - but
`pandas.DataFrame.head()` rejects a float positional count outright
(`TypeError: cannot do positional indexing ... of type float`). The
`staubility_variant` "fast"-mode siblings of both had the identical bug.
`eso.detection.methylation.find_motif_sites` already guarded this correctly
with `order[:int(num_sites)]` - the inconsistency across otherwise-parallel
functions is what made this worth checking directly rather than assuming
it worked. Confirmed the crash directly (a real `--num-sites 5`-style call,
not a mock) before fixing; fixed by adding the same `int(num_sites)` cast at
all four remaining call sites. The reason this was never caught before: the
one existing end-to-end test (`tests/test_pipeline_integration.py`) calls
`eso.pipeline.main()` directly with an already-int `num_sites=50`, bypassing
argparse's float conversion entirely.

Also added `tests/test_cli.py` (10 tests - `--num-sites` end-to-end via the
real detection pipeline, plus `--no-optimize`, `--mini-gc`/`--maxi-gc`,
`--recombination-mode`/`--slippage-mode`, `--organism-name`/`--method`, and
argparse's own choice validation), on top of the pre-existing
`tests/test_cli_custom_score.py`. No further issues found in `cli.py`
itself - it's a thin, mostly-mechanical argparse-to-`pipeline.main()`
wrapper.

## Follow-up: closing out the four remaining known gaps

Four items were left open after the review above - one residual limitation
and three known-but-deferred items. All four are now addressed.

**A dropped constraint now warns the user, with a reason.** Previously,
when `optimization_engine`'s retry loop dropped a constraint (see the
DNAChisel-crash section above), it did so silently - the only visible
symptom was `num_edits` coming back lower than expected, with nothing
telling the user *which* site survived or *why*. Added
`_warn_dropped_constraint`, called from both places the retry loop drops a
constraint, using Python's `warnings.warn`: it names the specific
constraint and gives a one-paragraph, minimal explanation (disrupting the
site would require a change that conflicts with another hard constraint,
most commonly translation preservation - e.g. a homopolymer landing on a
Met/Trp codon that has no synonymous alternative at all). Verified directly
that the T-homopolymer scenario from the earlier section now emits one
`UserWarning` per dropped position, and added a regression test asserting
both the warning fires and its explanation mentions translation
preservation.

**Two cosmetic warnings, fixed.** `eso.codon_usage._load_bundled_csv_cub`'s
`df.groupby('aa').apply(...)` triggered a pandas `FutureWarning` about
grouping columns being included in the group operation - fixed with
`include_groups=False` (the lambda only ever used the `codon`/
`freq_within_aa` columns, never `aa` itself, so this changes nothing about
the result). `eso.io_utils.exclusion_gc_tester` called
`biotools.gc_content('')` whenever a sliding window genuinely didn't
overlap any exclusion region (confirmed this is a real, reachable path, not
hypothetical - e.g. a window entirely before the region under test) - `0/0`
raised a `RuntimeWarning` and returned `NaN`. The `NaN` turned out to be
harmless in practice (it's immediately multiplied by `len(overlap) == 0`
inside a `max()`/`min()` comparison, and both of Python's `max`/`min` treat
a `NaN` second argument as "not greater/less", so it's silently dropped
without corrupting `curr_max_gc`/`curr_min_gc`) - but relying on that
comparison-with-`NaN` quirk to be harmless isn't something to depend on, so
fixed by skipping the `gc_content` call entirely when `overlap` is empty
(`curr_gc = 0.0` instead). Verified both fixes directly by re-running the
affected tests with `-W error::FutureWarning -W error::RuntimeWarning` -
they now pass clean.

**`--indexes-file`, a CLI flag for ORF/exclusion regions.** `indexes` (which
regions of each sequence are the translation-preserving ORF, and which are
locked from editing entirely) was only reachable by calling
`eso.pipeline.main` directly from Python - there was no CLI equivalent,
since the parameter's native shape (a dict keyed by `(file_stem,
seq_index)` tuples) doesn't serialize to JSON as-is. Added
`eso.io_utils.load_indexes_from_file` (JSON list of `{"file", "seq_index",
"orf_regions", "exclusion_regions"}` objects, converted into the dict shape
`pipeline.main` expects) and `--indexes-file` in `cli.py`, following the
same eager-validation-with-a-friendly-message pattern as
`--custom-score-file`/`CustomScoreFileError`. Added
`examples/indexes_template.json` as a copyable starting point and a new
README section. Verified with unit tests for the loader itself (missing
file, invalid JSON, wrong top-level type, missing required keys) and an
end-to-end CLI test confirming a declared exclusion region genuinely
survives optimization untouched when set via `--indexes-file`.

## Packaging - `python-Levenshtein`'s version constraint silently blocked a fresh install

Found while verifying the repo is actually installable by someone new,
ahead of handing it off. `pyproject.toml` declared
`python-Levenshtein = "^0.25.0"` - but Poetry's caret operator treats
pre-1.0 versions specially: for a `0.x.y` version, `^0.25.0` means
`>=0.25.0,<0.26.0`, not "any 0.25+ release" the way it would for a `^1.x.y`
constraint. Confirmed directly: a genuinely fresh `pip install
<this-repo>` in a brand-new virtualenv resolved `python-Levenshtein` to
`0.25.1` (the newest version inside that narrow range), which requires
`Levenshtein==0.25.1` exactly - and that release predates Python 3.13
(released after it), so no prebuilt wheel exists for it on this platform.
pip fell back to building from source, which failed outright without
Visual Studio installed, with a fairly opaque `CMake`/`scikit-build` error
that gives no hint the actual problem is a version-constraint syntax issue
three layers up. Fixed by widening the constraint to
`>=0.26.0,<1.0.0` (0.26.0 is the first version with prebuilt wheels
covering current Python versions). Re-verified end-to-end from a clean
virtualenv after the fix: `pip install <repo>` succeeds with no compiler
needed, `eso-optimize --help` and a real optimization run both work, the
bundled antibody example (`examples/antibody_optimization/run_example.py`)
runs successfully, and the full test suite passes (158 passed, 1 skipped -
the skip is the optional `docx-report` test, expected since `python-docx`
isn't installed by the base install).

## Two genuine cross-platform bugs (Mac/Linux), found via targeted audit

Everything up to this point had only ever been tested on Windows. Ahead of
handing the repo to a colleague on an unknown OS, an audit specifically for
Windows-only assumptions across `eso/` found two real, reachable bugs (no
hardcoded path separators, no `os.system`/`subprocess`, no Windows-only env
vars/registry access, and no hardcoded line-ending assumptions - those were
all already fine).

**File discovery silently dropped uppercase/mixed-case files on Mac/Linux.**
`eso.io_utils.relevant_file_paths` matched files via
`glob('*.fasta')`-style patterns against lowercase-only extensions in
`FILE_ENDINGS`. `glob`'s case sensitivity is filesystem-dependent, not
Python-version-dependent: Windows' default NTFS/FAT filesystems are
case-insensitive, so `glob('*.fasta')` there also matches `GENE.FASTA` -
but the identical call on a Mac/Linux (case-sensitive) filesystem does
*not*. Confirmed directly: a file named `GENE.FASTA` was found on this
(Windows) machine, and there is no "0 files found" check anywhere in
`pipeline.main`/`test_input`, so on Mac/Linux this would have been silent
data loss - the pipeline runs, reports `'Success!'`, and simply never
touches the file, with no error or warning. Fixed by rewriting
`relevant_file_paths` to enumerate directory entries directly and match
extensions case-insensitively in Python (`_matching_filetype`), rather than
relying on the filesystem's own case-folding behavior - verified directly
that `GENE.FASTA`, `Mixed.Fa`, and a nested `NESTED.GBK` are all now found,
and added a regression test.

**Missing explicit `encoding=` on several `open()` calls.** FASTA/GenBank
files (`file_opener`), `--indexes-file` JSON, `--motifs-path` MEME files,
`final_sequence.txt` output, and the bundled codon-usage CSVs were all
opened without `encoding=`, meaning Python falls back to
`locale.getpreferredencoding()` - `cp1252` by default on Windows, `utf-8`
almost everywhere else. GenBank records in particular routinely carry
non-ASCII metadata (accented author/species names in `DEFINITION`/`SOURCE`
lines, "µ", "°", etc.), so a file saved as UTF-8 (the near-universal
default outside Windows) and then read on Windows without an explicit
encoding - or the reverse - risks a `UnicodeDecodeError` or silently
mis-decoded characters. This is the same bug class already hit and fixed
once before in this session's own diagnostic scripts (see the DNAChisel
crash-investigation section above), just never swept across the actual
package source until now. Fixed by adding `encoding="utf-8"` to every
`open()`/`gzip.open()` call in the package (`io_utils.py`, `pipeline.py`,
`codon_usage.py`, `detection/methylation.py`) - no test added for this one
specifically (it would require actually exercising a non-UTF-8-locale
environment to observe the original failure), but it's a mechanical,
low-risk change and the existing suite continues to pass unchanged.

Also added to the README: an explicit "Quickstart" section at the very top
(exact prerequisites, one canonical install command, a copy-pasteable
first run, how to read the result) and a "Troubleshooting" section for the
most likely first-run failure modes (`eso-optimize` not found → use
`python -m eso.cli` instead; a Visual-Studio/CMake build error → upgrade
pip; wrong Python for the install; empty results). Previously, install and
usage instructions were split across several sections and assumed the
reader already knew what `poetry`, a virtualenv, or a CLI flag was.
