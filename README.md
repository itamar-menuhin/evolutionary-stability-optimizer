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
- **Methylation motifs** - sequence motifs (given as a PSSM/MEME file) recognized by
  host methylation machinery, which can trigger repair-associated mutation.

Two independently-developed implementations of the slippage/recombination detectors are
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
speed choice).

There's also `eso.detection.junction_linker`, which applies the same hotspot-detection
idea to a different problem: checking whether joining a target gene to a host sequence
via a linker creates a *new* hotspot right at the junction.

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

## Development

```bash
poetry install --with dev
pytest
```
