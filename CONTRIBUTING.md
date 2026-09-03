# Contributing to ESO

Thanks for considering a contribution. This document covers how to get set
up, how this codebase likes to work, and what a good pull request looks
like here.

## Getting set up

```bash
git clone https://github.com/itamar-menuhin/evolutionary-stability-optimizer.git
cd evolutionary-stability-optimizer
poetry install --with dev
pytest
```

No Poetry? `pip install -e ".[docx-report]" && pip install pytest hypothesis`
works too. See the README's Quickstart/Troubleshooting sections if the
install itself doesn't work - that's a documentation bug, please open an
issue.

## Before you open a PR

- **Run the full test suite** (`pytest`) and make sure it's green on your
  change. CI runs it again across Ubuntu/Windows/macOS x Python
  3.11/3.12/3.13 (`.github/workflows/tests.yml`) - if something only fails
  on one OS or Python version, that's worth understanding before merging,
  not working around.
- **Every real bug fix needs a regression test.** Not just "does the code
  run" - a test that fails on the old code and passes on the new code, with
  a comment/docstring explaining the concrete scenario that triggers it (a
  specific input, not just "edge case"). Look at any `test_*.py` file in
  `tests/` for the pattern this repo already uses everywhere.
- **Property-based tests for anything with a real invariant.** Several
  modules here (`eso.detection._overlap`, `eso.detection.slippage`) have
  `tests/test_*_properties.py` files using
  [Hypothesis](https://hypothesis.readthedocs.io/) to fuzz random inputs
  against invariants ("no two output ranges overlap," "every candidate's
  `sequence` field matches the real substring at its own coordinates") -
  this repo's actual bug history is mostly boundary/off-by-one errors that
  hand-picked examples missed the first time, so a new piece of range or
  index arithmetic is a good candidate for one of these, not just an
  example-based test.
- **Explain the *why*, not the *what*, in comments and docstrings.** A
  reader can already see what the code does; what they can't see is why an
  edge case is handled a particular way, what alternative was tried and
  rejected, or what real bug a piece of logic guards against. Most
  docstrings in this repo follow that pattern - match it.
- **Don't silently paper over a correctness risk with a safer default.**
  If you find a case where behavior is wrong or undefined, fix the actual
  logic (or explicitly reject the input with a clear error/warning) rather
  than picking a default that happens to avoid the symptom. If you're not
  sure a fix is complete, say so in the PR description rather than
  presenting it as settled.
- **Record real investigations in `docs/detector-comparisons.md`**, not just
  the final fix. This file is a running log of bugs found, benchmarks run,
  and even wrong turns later corrected - it exists so a future contributor
  (or AI assistant) can see *why* a decision was made instead of
  re-litigating it. Don't sanitize it down to only the right answer if the
  path there mattered.

## Code style

- No enforced linter/formatter currently - match the surrounding file's
  style (this repo generally avoids one-letter variable names, keeps
  functions focused, and prefers explicit code over clever one-liners).
- Cite the actual source for any formula or constant that comes from a
  paper or external tool (see `eso.detection.recombination.calc_recombination_score`
  for the pattern: exact equation, exact table/section reference, and a
  note about where a *different* value appears elsewhere and why it's not
  used here).
- Don't add a feature flag, config option, or abstraction for a
  hypothetical future need - this repo's own history includes at least one
  investigated, benchmarked, and then fully *removed* feature (windowed
  custom scoring) after nothing showed it was actually needed. Prefer
  adding it back later, with evidence, over speculative generality now.

## Two independently-developed detector implementations

`eso.detection.recombination`/`eso.detection.slippage` and
`eso.detection.staubility_variant` are two separate implementations of the
same underlying detection problem (see `eso.detection.dispatch` and
`docs/detector-comparisons.md` for why both exist and how they're
compared). If you're fixing a bug in one, check whether the same bug class
exists in the other - several fixes in this repo's history were found in
one implementation and then confirmed (or ruled out) in the other via
direct comparison, not assumed to be isolated.

## Adding a new organism / codon usage table

See `eso.codon_usage.CODON_USAGE_TABLES` - add a new entry there rather
than special-casing the organism elsewhere; `organism_name` already
supports arbitrary NCBI TaxIDs/species names via `python-codon-tables` for
anything in that database, so a custom table is only needed for organisms
outside it.

## Reporting a bug

Open a GitHub issue with: the exact input sequence (or the smallest one
that reproduces it), the exact command/function call, and the full error
message or unexpected output. If you can pin down which specific
detector/mode is involved (`--recombination-mode`/`--slippage-mode`), that
narrows it down significantly - see `eso.detection.dispatch`'s docstrings
for what each mode does differently.

## License

By contributing, you agree your contribution is licensed under this
repo's MIT license (see `LICENSE`).
