"""End-to-end smoke test: detection -> constraint conversion -> DNAChisel
optimization -> file output, via eso.pipeline.main().
"""

from eso.pipeline import main, suspect_site_extractor


def test_suspect_site_extractor_exposes_raw_alongside_collapsed():
    # regression test for the overlap-collapse coverage-gap fix: optimization
    # must be fed the uncollapsed, num_sites-unlimited view (df_*_raw), not
    # just the reported/collapsed one - see eso.pipeline.backend and
    # docs/detector-comparisons.md.
    seq = "ATG" + "GCT" * 17 + "TAA"

    result = suspect_site_extractor(seq, compute_motifs=False, num_sites=1)

    assert 'df_recombination_raw' in result
    assert 'df_slippage_raw' in result
    # num_sites=1 limits the reported view...
    assert len(result['df_slippage']) <= 1
    # ...but never the raw view, which is what optimization actually uses.
    assert len(result['df_slippage_raw']) >= len(result['df_slippage'])


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
    # regression test: the "_corrected" CSVs (every candidate actually given a
    # correction constraint, not just the collapsed report view) must be
    # written alongside the report CSVs whenever optimize=True.
    assert "recombination_sites_corrected.csv" in out_files or "slippage_sites_corrected.csv" in out_files


def test_no_corrected_csvs_written_when_optimize_is_off(tmp_path):
    # the "_corrected" CSVs describe what optimization actually touched - with
    # optimize=False nothing was corrected, so writing them would be misleading.
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    output_dir = tmp_path / "output"

    seq = "ATG" + "GCT" * 17 + "TAA"
    (input_dir / "test_gene.fasta").write_text(f">test_gene\n{seq}\n")

    message, _ = main(
        input_folder=str(input_dir),
        output_path=str(output_dir),
        compute_motifs=False,
        num_sites=50,
        optimize=False,
    )

    assert message == 'Success!'
    out_files = {p.name for p in (output_dir / "test_gene").iterdir()}
    assert "recombination_sites_corrected.csv" not in out_files
    assert "slippage_sites_corrected.csv" not in out_files
