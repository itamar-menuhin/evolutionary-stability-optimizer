# Detector comparisons and mode decisions

ESO carries independently-developed implementations of several hotspot detectors
(see the main README and the repo's provenance annotations). This doc records,
one detector type at a time, what was actually tested, what bugs were found and
fixed along the way, and the resulting recommendation for which mode to use
when.

## Recombination detection

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

### Benchmark

One embedded repeat per sequence length, wall-clock time to detect it:

| Length (nt) | `thorough` | `fast` | Speedup |
|---|---|---|---|
| 300 | 1.481s | 0.028s | 54x |
| 500 | 2.383s | 0.025s | 95x |
| 900 | 3.844s | 0.021s | 180x |
| 1,700 | 16.277s | 0.051s | 323x |
| 3,300 | 19.134s | 0.024s | 791x |
| 6,500 | 76.354s | 0.040s | 1,902x |

`thorough` grows roughly quadratically with sequence length (consistent with
its per-site candidate-generation-and-matching approach scanning the whole
sequence); `fast` stays flat regardless of length in this range (consistent
with a single vectorized counting pass).

### Recommendation

- **`thorough` (default)** for gene-length sequences (hundreds to low-thousands
  of nt - ESO's typical case), where a near-duplicate hotspot matters
  biologically and the ~1-20s cost is a non-issue.
- **`fast`** for much longer sequences (multi-kb constructs, whole plasmids),
  where `thorough`'s quadratic-ish cost becomes minutes-to-hours - accepting
  that it can miss a centrally-mutated near-duplicate.
- Kept as two explicit modes rather than merged into one canonical
  implementation, per project decision (2026-07) - revisit if/when the
  slippage and methylation detector pairs are worked through and a pattern for
  unification (or not) emerges.

---

*Slippage and methylation detector comparisons: not yet done - see the open
decisions list in the repo status summary.*
