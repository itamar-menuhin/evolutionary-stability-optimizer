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

### Cross-validated against the authoritative reference implementation - and found a real bug

Prompted by a direct ask to look for libraries or reference sources that
would *increase trust*, not just speed - the strongest trust signal
available wasn't a library at all, but the actual reference tool this whole
scoring model is based on: the EFM Calculator
([github.com/barricklab/efm-calculator](https://github.com/barricklab/efm-calculator)),
whose `get_recombo_rate()` implements the same formula shape
`calc_recombination_score` does. Comparing them directly (fetched the source
via `gh`/web fetch, not from memory) found a real discrepancy: our constant
was `a=5.8`; the reference uses `a=8.8`.

Checked the reference repo's full commit history before concluding anything
- `8.8` has been the value since the very first commit (2015-03-31),
already attributed in that commit's own docstring to "the Oliviera, et al.
formula" (a citation distinct from the Jack et al. 2015 EFM Calculator paper
this codebase's docstrings cite for the tool overall). There is no point in
the reference tool's history where `5.8` was ever correct - so this isn't a
case of two legitimately different published versions; `5.8` appears to be a
plain transcription error that entered somewhere in ESO_curr/STABLES'
history and was faithfully carried through the port.

**Quantified impact**: the `a=5.8` version systematically *overestimates*
recombination risk (less-negative log10 score) by up to ~0.29 at short
range/site length (e.g. location_delta=1, site_length=16: -4.675 vs. the
correct -4.963) - enough to flip a result right at the `-9` filter cutoff,
i.e. a real false-positive risk at the margin, not just a cosmetic score
difference.

Fixed in both `eso.detection.recombination.calc_recombination_score` and
`eso.detection.staubility_variant.find_recombination_sites` (the formula is
duplicated between them). Verified two ways: `calc_recombination_score` now
matches the reference formula to floating-point precision
(`test_calc_recombination_score_matches_efm_calculator_reference`), and a
100-trial fuzz sweep confirms the change doesn't affect the underlying
near-duplicate *detection* ability (still 0 misses) - only the risk score
itself, as intended.

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
table shows, and the crossover point for reaching for `fast` moved from
~1-2kb out to the tens of kb. See below for current numbers and updated
recommendation.

### Recommendation (updated after the optimization below)

- **`thorough` (default)** for anything from gene-length up through tens of kb
  (see updated benchmark below) - `fast`'s sensitivity gap (missing a
  centrally-mutated near-duplicate) is a real cost, so there's no reason to
  pay it until `thorough`'s runtime actually becomes inconvenient.
- **`fast`** once sequences get long enough that even the optimized
  `thorough` cost matters (very large multi-kb+ constructs, whole plasmids,
  or workloads processing many sequences where every second compounds).
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
put it at tens of minutes. This is why the recommendation above changed: the
practical crossover to `fast` moved from ~1-2kb out to the tens of kb.

**Re-verified after the `ranges_overlap` and `a=5.8→8.8` fixes below**, since
neither is a no-op on timing even though neither changes algorithmic
complexity: `a=8.8` gives systematically lower (more negative) scores, so
fewer marginal candidates survive the `-9` filter before the expensive
elongation/collapse steps, which measurably *helps*:

| Length (nt) | `thorough` (current) | `fast` (current) |
|---|---|---|
| 400 | 0.255s | 0.017s |
| 3,400 | 0.431s | 0.014s |
| 6,600 | 0.748s | 0.020s |
| 13,000 | 1.099s | 0.041s |
| 25,800 | 2.202s | 0.077s |
| 51,400 | 4.311s | 0.219s |

Consistently faster than the table above (e.g. 51,300nt: 6.527s -> 4.311s),
but the *shape* of the recommendation is unchanged - `thorough` vs. `fast`'s
ratio is still ~20-27x in this range, so the crossover guidance above still
holds exactly as stated.

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

*Methylation detector comparison: not yet done - see the open decisions list
in the repo status summary.*
