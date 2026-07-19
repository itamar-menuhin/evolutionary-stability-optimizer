#!/usr/bin/env python3
"""Example: optimize a human antibody heavy/light chain against recombination,
slippage, and GC-content constraints, using the bundled human-antibody codon
usage tables.

Usage: python run_example.py
"""

from pathlib import Path

from eso.pipeline import main as run_pipeline

DATA_DIR = Path(__file__).parent / "data"


def _write_fasta(text_path, fasta_path, record_id):
    seq = text_path.read_text().strip()
    with open(fasta_path, 'w') as f:
        f.write(f">{record_id}\n")
        for i in range(0, len(seq), 70):
            f.write(seq[i:i + 70] + "\n")


def main():
    output_dir = Path(__file__).parent / "output"

    # Each chain gets its own input subdirectory - eso scans the whole
    # input_folder for FASTA files, so mixing both chains in one directory
    # would optimize each under the other's codon-usage table too.
    heavy_dir = DATA_DIR / "heavy_chain"
    light_dir = DATA_DIR / "light_chain"
    heavy_dir.mkdir(exist_ok=True)
    light_dir.mkdir(exist_ok=True)

    _write_fasta(DATA_DIR / "HeavyCORF.txt", heavy_dir / "heavy_chain.fasta", "heavy_chain_antibody")
    _write_fasta(DATA_DIR / "LightCORF.txt", light_dir / "light_chain.fasta", "light_chain_antibody")

    print("Optimizing heavy chain (human_antibody_heavy_chain codon usage)...")
    message, heavy_results = run_pipeline(
        input_folder=str(heavy_dir),
        output_path=str(output_dir / "heavy_chain"),
        compute_motifs=False,
        num_sites=50,
        optimize=True,
        mini_gc=0.3,
        maxi_gc=0.7,
        method='use_best_codon',
        organism_name='human_antibody_heavy_chain',
    )
    print(message)

    print("Optimizing light chain (human_antibody_light_chain codon usage)...")
    message, light_results = run_pipeline(
        input_folder=str(light_dir),
        output_path=str(output_dir / "light_chain"),
        compute_motifs=False,
        num_sites=50,
        optimize=True,
        mini_gc=0.3,
        maxi_gc=0.7,
        method='use_best_codon',
        organism_name='human_antibody_light_chain',
    )
    print(message)

    print(f"\nResults written to: {output_dir}")
    print("See 'final_sequence.txt' in each subdirectory for the optimized sequence, "
          "and 'recombination_sites.csv'/'slippage_sites.csv' for detected hotspots.")


if __name__ == "__main__":
    main()
