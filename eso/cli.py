"""Command-line entry point: `eso-optimize --input-folder ... --output-path ...`."""

import argparse

import numpy as np

from eso.pipeline import main as run_pipeline


def build_parser():
    parser = argparse.ArgumentParser(
        description="Detect and remove hypermutable sites (recombination, slippage, "
                    "methylation hotspots) from DNA sequences while preserving translation.")
    parser.add_argument('--input-folder', default='.', help="Directory containing FASTA/GenBank files.")
    parser.add_argument('--output-path', default=None, help="Directory to write results into (default: <input-folder>/output).")
    parser.add_argument('--compute-motifs', action='store_true', help="Also detect methylation motif sites.")
    parser.add_argument('--motifs-path', default=None, help="Path to a MEME-minimal-format PSSM file (required with --compute-motifs).")
    parser.add_argument('--num-sites', type=float, default=np.inf, help="Max hotspots to report/constrain per category (default: all).")
    parser.add_argument('--no-optimize', action='store_true', help="Only detect hotspots, skip sequence optimization.")
    parser.add_argument('--mini-gc', type=float, default=0.3)
    parser.add_argument('--maxi-gc', type=float, default=0.7)
    parser.add_argument('--method', default='use_best_codon',
                        choices=['use_best_codon', 'match_codon_usage', 'harmonize_rca'])
    parser.add_argument('--organism-name', default='not_specified',
                        help="Host organism for codon optimization (species name, TaxID, or a custom table name "
                             "e.g. kompas/human_antibody_heavy_chain/human_antibody_light_chain).")
    parser.add_argument('--recombination-mode', default='thorough', choices=['thorough', 'fast'],
                        help="'thorough' (default): Levenshtein-tolerant, catches near-duplicate hotspots, "
                             "best for gene-length sequences. 'fast': exact-match only, ~100x faster, "
                             "for much longer sequences (see eso.detection.dispatch).")
    parser.add_argument('--slippage-mode', default='default', choices=['default', 'fast'],
                        help="Both detect identical hotspots (unlike --recombination-mode); "
                             "'fast' just scales better for longer (multi-kb+) sequences.")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    message, results = run_pipeline(
        input_folder=args.input_folder,
        output_path=args.output_path,
        compute_motifs=args.compute_motifs,
        num_sites=args.num_sites,
        motifs_path=args.motifs_path,
        optimize=not args.no_optimize,
        mini_gc=args.mini_gc,
        maxi_gc=args.maxi_gc,
        method=args.method,
        organism_name=args.organism_name,
        recombination_mode=args.recombination_mode,
        slippage_mode=args.slippage_mode,
    )
    print(message)
    return 0 if message == 'Success!' else 1


if __name__ == '__main__':
    raise SystemExit(main())
