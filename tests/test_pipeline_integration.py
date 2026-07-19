"""End-to-end smoke test: detection -> constraint conversion -> DNAChisel
optimization -> file output, via eso.pipeline.main().
"""

from eso.pipeline import main


def test_optimizes_out_repetitive_codon_hotspot(tmp_path):
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    output_dir = tmp_path / "output"

    # 17x the same Ala codon in a row is a severe slippage/recombination hotspot.
    seq = "ATG" + "GCT" * 17 + "TAA"
    (input_dir / "test_gene.fasta").write_text(f">test_gene\n{seq}\n")

    message, results = main(
        input_folder=str(input_dir),
        output_path=str(output_dir),
        compute_motifs=False,
        num_sites=50,
        optimize=True,
        mini_gc=0.3,
        maxi_gc=0.7,
        method='use_best_codon',
        organism_name='kompas',
    )

    assert message == 'Success!'
    assert len(results) == 1

    _, seq_index, optimized_seq = results[0]
    assert seq_index == 0
    assert len(optimized_seq) == len(seq)
    # same amino acid translation preserved
    assert optimized_seq[:3] == "ATG"
    # codons diversified away from the all-GCT hotspot
    assert optimized_seq != seq

    out_files = {p.name for p in (output_dir / "test_gene").iterdir()}
    assert "final_sequence.txt" in out_files
    assert "recombination_sites.csv" in out_files or "slippage_sites.csv" in out_files
