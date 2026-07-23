# ESO - Evolutionary Stability Optimizer

[![Tests](https://github.com/itamar-menuhin/evolutionary-stability-optimizer/actions/workflows/tests.yml/badge.svg)](https://github.com/itamar-menuhin/evolutionary-stability-optimizer/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)](pyproject.toml)

ESO detects hypermutable sites in engineered DNA sequences and optimizes them away with
[DNAChisel](https://github.com/Edinburgh-Genome-Foundry/DNAChisel), while preserving the
amino-acid translation and optimizing for host codon usage. This implementation follows
the approach introduced in Menuhin-Gruman et al. (2022, *ACS Synthetic Biology*) - see
[Citation](#citation) below.

Genes built from repetitive or duplicated DNA elements (a side effect of standard codon
optimization, which tends to reuse the same "best" codon repeatedly) are prone to mutate
away during propagation in a host organism, through mechanisms like replication slippage
and recombination-mediated deletion. ESO detects these hotspots and asks DNAChisel to
route around them while it optimizes, using the empirical mutation-rate model from the
[EFM Calculator](https://doi.org/10.1021/acssynbio.5b00068) (Jack et al., 2015, ACS
Synthetic Biology).

## Quickstart

**1. Requirements**: Python 3.11 or newer. Check what you have installed:

```bash
python --version
```

If that prints something below `Python 3.11`, or fails with "command not found"/"not
recognized", install a current version from [python.org/downloads](https://www.python.org/downloads/)
first (the installer's default settings are fine) before continuing.

**2. Install ESO.** In a terminal, `cd` into this folder (the one this README is in), then:

```bash
pip install .
```

This can take a minute or two the first time. If it ends with a line that doesn't contain
the word `error`, it worked - skip to step 3. If you do see an error, check
[Troubleshooting](#troubleshooting) below before asking for help.

**3. Try it on a real sequence.** Put a FASTA file (a plain text file starting with a `>`
line, then the DNA sequence) in a folder by itself, e.g. `my_sequences/gene.fasta`
containing:

```
>my_gene
ATGGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTTAA
```

Then run:

```bash
eso-optimize --input-folder my_sequences --output-path my_results
```

**4. Check the result.** Look inside `my_results/gene/` - `final_sequence.txt` has the
optimized sequence, `recombination_sites.csv`/`slippage_sites.csv` list what was detected
and fixed. If the very last line printed to the terminal was `Success!`, it worked.

That's the whole loop. Everything below this point is reference material for going
further (different host organisms, custom scoring, locking regions from editing, and so
on) - not required to get a first result.

If you'd rather be walked through setup interactively instead of reading this file, try the
[ESO Onboarding Assistant](https://chatgpt.com/g/g-6a60d7b081108191a8bc208b89958267-eso-onboarding-assistant),
a ChatGPT assistant configured specifically for onboarding onto this tool (see
[`docs/custom-gpt-setup.md`](docs/custom-gpt-setup.md) for its configuration).

## Troubleshooting

- **`'eso-optimize' is not recognized` / `command not found: eso-optimize`** (common on
  Windows right after a fresh install): the command was installed, but its folder isn't on
  your terminal's PATH yet. Use this instead - it always works, no PATH needed:

  ```bash
  python -m eso.cli --input-folder my_sequences --output-path my_results
  ```

  (Every `eso-optimize ...` example in this README works identically as
  `python -m eso.cli ...`.)
- **An error mentioning `Microsoft Visual Studio`, `CMake`, or "failed building wheel"
  during `pip install .`**: this means `pip` tried to compile a dependency from source
  instead of using a prebuilt version - almost always fixed by upgrading pip first
  (`python -m pip install --upgrade pip`) and trying the install again, since an older
  `pip` can miss prebuilt wheels that a newer one finds.
- **`ModuleNotFoundError` or `ImportError` right after install**: double check you're
  running `python`/`eso-optimize` from the *same* Python installation you ran
  `pip install .` with. If you're not sure, run
  `python -m pip show eso` - if that fails, you installed into a different Python than the
  one you're now running.
- **"No such file or directory" for your input folder**: `--input-folder` is relative to
  wherever your terminal's current directory is - either `cd` there first, or use a full
  path (e.g. `C:\Users\you\my_sequences` or `/Users/you/my_sequences`).
- **Nothing in the output folder / an empty `results` list**: ESO looks for files ending
  in `.fasta`/`.fna`/`.ffn`/`.faa`/`.frn`/`.fa`/`.gb`/`.gbk`/`.genbank` (optionally
  gzipped), directly inside `--input-folder` or one level under it - check your file's
  extension matches one of those.
- Anything else: the error messages this tool prints are meant to be read directly and
  acted on (not just a Python traceback to decode) - if one isn't clear, that's a bug in
  the tool itself, worth reporting rather than working around.

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

## Using ESO as a library, instead of files

`main()` (above) is file-in, file-out - convenient for a one-off CLI run, but if you
already have sequences in memory as part of your own code (e.g. generated, fetched from a
database, or produced by an earlier step in your own pipeline), you don't need to write
them to disk first. `eso.optimize.optimization_engine` (also importable as
`eso.optimization_engine`) takes a plain DNA string and returns a plain DNA string - no
files involved:

```python
from eso import optimization_engine, suspect_site_extractor

seq = "ATG" + "GCT" * 15 + "TAA"  # any DNA string you already have

# 1. detect hotspots (skip this step and the dataframes below entirely if you only
#    want codon/GC optimization, with no hotspot avoidance)
sites = suspect_site_extractor(seq, compute_motifs=False, num_sites=50)

# 2. optimize, avoiding what was detected
final_seq, objectives_summary, num_edits = optimization_engine(
    seq,
    organism_name="e_coli",
    df_recombination=sites["df_recombination"],
    df_slippage=sites["df_slippage"],
)
```

`final_seq` is a plain `str` you can feed straight back into whatever your own code does
next (write it out yourself, pass it to another function, etc.) - nothing here touches
the filesystem. `suspect_site_extractor` (also importable as
`eso.suspect_site_extractor`) is the same detection step `main()` runs internally; call
it on its own if you only want the hotspot dataframes, with no optimization at all.

This composes directly with a custom scoring function - just pass `custom_score_fn=...`
to `optimization_engine` instead of `organism_name`. **If you're plugging in your own
model** (an ML model, or any function that scores the sequence as a whole rather than one
codon at a time - true for most real models), leave `custom_score_window` unset:

```python
final_seq, _, num_edits = optimization_engine(
    seq,
    custom_score_fn=my_model.predict,  # any function: whole seq (str) -> a number, higher = better
    df_recombination=sites["df_recombination"],
    df_slippage=sites["df_slippage"],
)
```

**Do not** pass `custom_score_window=3` (or any number) unless your score is *provably*
just a sum of independent per-codon contributions - the same way built-in CAI/tAI scoring
works. Passing a window when that assumption doesn't hold isn't just "a bit less
accurate" - it silently computes something structurally unrelated to your real score
(confirmed directly: a model rewarding the sequence for containing a long repeated run
anywhere scored a real 9-nucleotide run as `9` in the correct, unwindowed mode, but as
`15` - a meaningless sum of per-codon maxes - when a window was wrongly applied), with no
error or warning, since DNAChisel has no way to know your window choice was wrong. The
next section explains exactly when a window *is* safe to use (mostly: never, for an
external model - it's really only for scores that are inherently per-codon, like the
built-in codon-usage scoring itself) and the speed/correctness trade-off it's for.

See `optimization_engine`'s docstring (`eso/optimize.py`) for every parameter
(`mini_gc`/`maxi_gc`, `orf_regions`/`exclusion_regions`, `method`, and so on) - everything
available via `main()`/the CLI is available here too, just without the file layer.

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
    custom_score_window=3,  # only correct because THIS score really is per-codon - see below
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
  any fixed window size works, as long as that decomposition assumption holds. **If the
  scored region's length isn't itself a multiple of N**, the trailing remainder is
  silently excluded from every `score_fn` call (only whole windows are ever scored) - ESO
  warns about this once per optimization run when it happens, so it won't pass silently.
  This can't occur for the common case (N=3 against a translation-preserving ORF, which is
  always a multiple of 3 already) - it's only reachable with an N that doesn't evenly
  divide your scored region's length.
- **`WINDOW = None`** (the default if omitted) - "global" mode: `score_fn` is called once
  on the entire sequence, from scratch, on *every* trial mutation during optimization.
  Always correct, no assumption about how the score decomposes, but can be very slow on
  long sequences or an expensive `score_fn` - a warning is raised when this mode is used,
  for exactly that reason.

Independently of `WINDOW`, `custom_score_minimize=True` (`--custom-score-minimize` on the
CLI) treats a *lower* `score_fn` value as better, instead of higher.

**Scope**: custom scoring is automatically restricted to `orf_regions` (one scored region
per ORF, matching how the built-in CAI/tAI codon-usage scoring is already scoped via
DNAChisel's `CodonOptimize(location=orf, ...)`) - `score_fn` never sees any non-ORF
flanking sequence (UTRs, locked/excluded regions). If you don't pass `orf_regions`, this
is the whole sequence (trimmed to a multiple of 3), same as everywhere else in ESO.

## Restricting ORF and exclusion regions per sequence

By default, the entire sequence is treated as one in-frame, translation-preserving ORF
with nothing locked. To instead give each sequence its own ORF region(s) (e.g. skip a
UTR) and/or exclusion regions that must never be edited (e.g. a known regulatory
element), pass `indexes` (from Python) or `--indexes-file` (from the CLI).

From the command line, point at a JSON file - copy
[`examples/indexes_template.json`](examples/indexes_template.json) as a starting point:

```json
[
  {
    "file": "my_gene",
    "seq_index": "0",
    "orf_regions": "1-6, 51-68",
    "exclusion_regions": "1-6, 50-68"
  }
]
```

```bash
eso-optimize --input-folder path/to/fasta_files --indexes-file indexes.json
```

- `file` is the FASTA/GenBank file's stem (no extension, e.g. `"my_gene"` for
  `my_gene.fasta`).
- `seq_index` is which record within that file, 0-indexed in file order, as a string.
- `orf_regions`/`exclusion_regions` are 1-indexed, inclusive region strings (e.g.
  `"1-6, 51-68"` for two separate regions); omit `exclusion_regions` (or use `""`) for no
  exclusions.

A malformed `--indexes-file` (bad JSON, missing `file`/`seq_index`) fails immediately with
a plain-English message; malformed region strings themselves are validated the same way
whether `indexes` came from a file or was passed directly to `eso.pipeline.main`/`main`.

From Python, pass the equivalent dict directly:

```python
from eso import main

message, results = main(
    input_folder="path/to/fasta_files",
    indexes={("my_gene", "0"): ("1-6, 51-68", "1-6, 50-68")},
)
```

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

## Citation

If you use this tool, please cite the paper it implements:

> Menuhin-Gruman, I., Arbel, M., Amitay, N., Sionov, K., Naki, D., Katzir, I., Edgar, O.,
> Bergman, S., & Tuller, T. (2022). Evolutionary Stability Optimizer (ESO): A Novel
> Approach to Identify and Avoid Mutational Hotspots in DNA Sequences While Maintaining
> High Expression Levels. *ACS Synthetic Biology*, 11(3), 1142-1151.
> https://doi.org/10.1021/acssynbio.1c00426

```bibtex
@article{menuhingruman2022eso,
  title   = {Evolutionary Stability Optimizer (ESO): A Novel Approach to Identify and
             Avoid Mutational Hotspots in DNA Sequences While Maintaining High
             Expression Levels},
  author  = {Menuhin-Gruman, Itamar and Arbel, Matan and Amitay, Niv and Sionov, Karin
             and Naki, Doron and Katzir, Itai and Edgar, Omer and Bergman, Shaked
             and Tuller, Tamir},
  journal = {ACS Synthetic Biology},
  volume  = {11},
  number  = {3},
  pages   = {1142--1151},
  year    = {2022},
  doi     = {10.1021/acssynbio.1c00426}
}
```

## License

[MIT](LICENSE)
