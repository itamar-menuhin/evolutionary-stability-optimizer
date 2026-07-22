# ESO - Evolutionary Stability Optimizer

ESO detects hypermutable sites in engineered DNA sequences and optimizes them away with
[DNAChisel](https://github.com/Edinburgh-Genome-Foundry/DNAChisel), while preserving the
amino-acid translation and optimizing for host codon usage.

Genes built from repetitive or duplicated DNA elements (a side effect of standard codon
optimization, which tends to reuse the same "best" codon repeatedly) are prone to mutate
away during propagation in a host organism, through mechanisms like replication slippage
and recombination-mediated deletion. ESO detects these hotspots and asks DNAChisel to
route around them while it optimizes, using the empirical mutation-rate model from the
[EFM Calculator](https://doi.org/10.1021/acssynbio.5b00068) (Jack et al., 2015, ACS
Synthetic Biology).

## What it detects

- **Replication slippage** - short tandem repeats (a base unit of length 1-15 repeated
  3+ times, or a single nucleotide repeated 4+ times) that polymerase can skip or
  duplicate during replication.
- **Recombination-mediated deletion (RMD)** - pairs of near-identical sites (16+ nt,
  within Levenshtein distance 1 of each other) that homologous recombination can delete
  between.
- **Methylation motifs** - sequence motifs recognized by host methylation machinery,
  which can trigger repair-associated mutation. Given as a PSSM/MEME file, a bundled
  common motif (E. coli's Dam/Dcm systems), or your own IUPAC consensus string
  (e.g. `"GATC"`) - no file needed for either of the latter two.

Two independently-developed implementations of the recombination/slippage detectors are
included (`eso.detection.recombination`/`eso.detection.slippage`, and
`eso.detection.staubility_variant`) - they haven't yet been reconciled into a single
canonical algorithm. Both are routed through `eso.detection.dispatch`
(`find_recombination_sites(seq, num_sites, mode="thorough" | "fast")` and
`find_slippage_sites(seq, num_sites, mode="default" | "fast")`), also exposed via
`eso.pipeline.main(..., recombination_mode=..., slippage_mode=...)` /
`eso-optimize --recombination-mode --slippage-mode`. See
[`docs/detector-comparisons.md`](docs/detector-comparisons.md) for the tradeoffs,
benchmarks, and bugs found while comparing them - recombination's two modes trade off
sensitivity for speed, but slippage's are fully equivalent in what they detect (a pure
speed choice). Methylation motif detection only has one implementation
(`eso.detection.methylation`) - a second was built and compared, but deleted after it
turned out to disagree in accuracy, not just speed (see the same doc for the writeup).

**Not included**: STABLES also carried a junction-linker hotspot checker
(`linker_suspect_utils.py`), applying the same slippage/recombination-detection idea to a
different, narrower problem - checking whether joining a target gene to a host sequence
via a linker creates a *new* hotspot at the junction. It was ported into an earlier version
of this repo (`eso.detection.junction_linker`) but removed: its detection thresholds (a
4-repeat homopolymer check, a 12-mer recombination check) don't match the calibrated -9
filter this codebase's other detectors use, and tracing them back into the original
STABLES `select_fusion_linkers.py` pipeline showed they act as a near-absolute veto in
STABLES' own linker-selection algorithm - deliberate, STABLES-specific conservatism for
screening many cheap, redesignable linker candidates, not a general-purpose hotspot
detector. That's a different problem from ESO's (evaluating one already-chosen sequence),
so it was judged out of scope for this library rather than reconciled with the other
detectors' thresholds.

## Install

```bash
poetry install
# or: pip install -e .
```

Word-document diff reports (`eso.report`) need the optional `docx-report` extra:

```bash
poetry install -E docx-report
```

## Usage

```python
from eso import main

message, results = main(
    input_folder="path/to/fasta_files",
    output_path="path/to/output",
    optimize=True,
    mini_gc=0.3,
    maxi_gc=0.7,
    method="use_best_codon",
    organism_name="e_coli",  # or a TaxID, or a bundled table (see eso.codon_usage)
)
```

Or from the command line:

```bash
eso-optimize --input-folder path/to/fasta_files --output-path path/to/output --organism-name e_coli
```

For each FASTA/GenBank file found in `input_folder`, this writes to `output_path/<file_stem>/`:

- `final_sequence.txt` - the optimized sequence, plus CAI-before/after and edit-count stats.
- `recombination_sites.csv` / `slippage_sites.csv` / `motif_sites.csv` - detected hotspots.
- `sequence_comparison.docx` - a diff view of original vs. optimized sequence (if the
  `docx-report` extra is installed).

`organism_name` accepts anything supported by
[python-codon-tables](https://github.com/Edinburgh-Genome-Foundry/python_codon_tables)
(a species name or NCBI TaxID), or one of the bundled custom tables in
`eso.codon_usage.CODON_USAGE_TABLES` (`C1`, `kompas`, `human_antibody_heavy_chain`,
`human_antibody_light_chain`) for hosts not in that database.

See [`examples/antibody_optimization`](examples/antibody_optimization) for a complete
worked example (human antibody heavy/light chain optimization).

## Scoring sequences your own way, instead of CAI/tAI

By default, optimization scores codon choices against a codon-usage table (CAI/tAI-style,
via `organism_name`). To score sequences with your own logic instead, write a Python file
defining a `score(seq)` function (`seq` is a plain DNA string like `"ATGCGT..."`; return a
number, higher = better) - copy
[`examples/custom_score_template.py`](examples/custom_score_template.py) as a starting
point, which also explains the `WINDOW` setting described below.

From the command line:

```bash
eso-optimize --input-folder path/to/fasta_files --custom-score-file my_score.py
```

(`--custom-score-file` overrides `--organism-name`/`--method`.) A `--custom-score-file`
with a mistake in it - a missing `score` function, a typo, a function that crashes or
returns the wrong type - fails immediately with a plain-English message, before any
optimization runs, rather than surfacing later as a Python traceback.

From Python, either call `optimization_engine`/`eso.pipeline.main` directly with
`custom_score_fn`/`custom_score_window`/`custom_score_minimize`, or reuse the same
file-loading + validation the CLI uses via `eso.custom_score.load_custom_score_from_file`:

```python
from eso.custom_score import load_custom_score_from_file
from eso.optimize import optimization_engine

score_fn, window = load_custom_score_from_file("my_score.py")  # same validation as the CLI
final_seq, _, _ = optimization_engine(seq, custom_score_fn=score_fn, custom_score_window=window)

# or, skipping the file entirely:
final_seq, _, _ = optimization_engine(
    seq,
    custom_score_fn=lambda codon: codon.count("G") + codon.count("C"),
    custom_score_window=3,
)
```

**The `WINDOW` setting** controls how `score_fn` gets called during optimization, and is
the main thing to get right:

- **`WINDOW = 3`** (or any positive integer N) - "windowed" mode: `score_fn` is called on
  each successive, non-overlapping N-nt chunk of the sequence, and the results are summed.
  Fast, because DNAChisel only has to re-score the chunks actually touched by a given trial
  edit, not the whole sequence - this is the same mechanism CAI/tAI-style per-codon scoring
  already uses (there, N=3). Only give a *correct* total score if your real score genuinely
  decomposes as a sum over fixed-size chunks - it doesn't have to be codon-sized (N=3);
  any fixed window size works, as long as that decomposition assumption holds.
- **`WINDOW = None`** (the default if omitted) - "global" mode: `score_fn` is called once
  on the entire sequence, from scratch, on *every* trial mutation during optimization.
  Always correct, no assumption about how the score decomposes, but can be very slow on
  long sequences or an expensive `score_fn` - a warning is raised when this mode is used,
  for exactly that reason.

Independently of `WINDOW`, `custom_score_minimize=True` (`--custom-score-minimize` on the
CLI) treats a *lower* `score_fn` value as better, instead of higher.

**Current limitation**: the underlying `eso.custom_score.CustomScore` also accepts a
`location` to restrict scoring to a sub-region of the sequence, but this isn't exposed via
the file-based loader or the CLI yet - only relevant if you need to score just part of a
sequence rather than the whole thing.

## Motif sources: methylation, cryptic ribosome binding, cryptic promoters, and your own

Motif detection (`--compute-motifs`) isn't only about methylation - the same PSSM-based
scanner (`eso.detection.methylation.find_motif_sites`) works for any short sequence motif
you want to flag. Needs at least one motif, from any combination of three sources:

- **A MEME-minimal-format PSSM file** (`--motifs-path` / `motifs_path=`) - the original
  option, for a curated or experimentally-derived motif set you already have as a file.
- **Bundled common motifs** (`--common-motifs dam,dcm` / `common_motifs=["dam", "dcm"]`)
  - no file needed. See [`eso/detection/common_motifs.py`](eso/detection/common_motifs.py)
  (`COMMON_MOTIFS`) for the full list and sources; currently:
  - **Methylation** (E. coli, already a first-class host here - see `eso.codon_usage`'s
    bundled `e_coli` table): `dam` (GATC, N6-methyladenine), `dcm` (CCWGG,
    C5-methylcytosine on the internal C).
  - **Cryptic ribosome binding**: `shine_dalgarno` (AGGAGG) - flags a copy of the
    bacterial RBS consensus occurring *inside* a coding region, a known source of
    unintended internal translation initiation.
  - **Cryptic bacterial promoter elements**: `sigma70_minus35` (TTGACA),
    `sigma70_minus10` (TATAAT) - the two sigma70 hexamers; an accidental occurrence of
    either inside a coding sequence is a classic source of unwanted transcription. Each
    hexamer is flagged independently (a real promoter needs both, correctly spaced,
    which this detector doesn't check for) - treat an isolated hit as a coarse screen,
    not a confirmed cryptic promoter.
  - **Not included, and why** (see the module docstring for the full explanation):
    transcription terminators (a secondary-structure property, not a fixed linear
    motif), the Kozak sequence (something to match, not avoid - a different problem),
    and restriction enzyme sites (already covered comprehensively by DNAChisel's own
    `dnachisel.list_common_enzymes()` / `EnzymeSitePattern`, usable directly with
    `AvoidPattern` during optimization - no need to duplicate it here).
- **Your own IUPAC consensus string** - the easy custom-motif path, no PSSM/MEME file
  needed:

  ```python
  from eso.detection.motif_utils import motif_from_consensus, motifs_from_consensus_dict

  my_motif = motif_from_consensus("my_site", "GANTC")  # N = any base, standard IUPAC codes
  # or define several at once:
  my_motifs = motifs_from_consensus_dict({"site_a": "GATC", "site_b": "CCWGG"})
  ```

  Pass the result as `relevant_motifs` to `eso.detection.methylation.find_motif_sites`
  directly, or combine with the other two sources yourself before calling it - there's no
  CLI flag for this one yet (a MEME file is still the CLI's path for anything beyond the
  bundled common motifs).

Most published methylation/restriction motifs (REBASE, NEB's technical notes, the primary
literature) are naturally given exactly this way - a short consensus sequence with IUPAC
ambiguity codes (`W` = A or T, `N` = any base, etc.) - not as a position-probability
matrix, so this avoids hand-authoring a MEME file just to check for one. For organisms
beyond E. coli, [REBASE](http://rebase.neb.com) is the standard reference database for
methylation motifs - look up the motif there and pass it straight to
`motif_from_consensus`.

## Development

```bash
poetry install --with dev
pytest
```
