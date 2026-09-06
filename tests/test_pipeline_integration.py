"""End-to-end smoke test: detection -> constraint conversion -> DNAChisel
optimization -> file output, via eso.pipeline.main().
"""

import pandas as pd

from eso import pipeline as pipeline_module
from eso.pipeline import main, reoptimize_until_stable, suspect_site_extractor

_REOPT_KWARGS = dict(
    compute_motifs=False, num_sites=50, motifs_path=None, common_motifs=None,
    recombination_mode='thorough', slippage_mode='default', mini_gc=0.3, maxi_gc=0.7,
    method='use_best_codon', organism_name='kompas', custom_score_fn=None,
    custom_score_minimize=False, orf_regions=(), exclusion_regions=(),
)


def _fake_sites(df_slippage_raw=None):
    empty = pd.DataFrame()
    slippage = df_slippage_raw if df_slippage_raw is not None else empty
    return {
        'df_recombination': empty, 'df_recombination_raw': empty,
        'df_slippage': slippage, 'df_slippage_raw': slippage,
    }


_HOTSPOT_DF = pd.DataFrame([{
    # 'sequence' must be the FULL repeated string (all num_base_units copies
    # of the base unit), not just one copy - modify_df_slippage (called from
    # optimize.py) slices it in `length_base_unit`-sized chunks per base unit
    # (see eso/detection/slippage.py:229-248); passing just one copy left
    # every chunk past the first as an empty string, an effectively-empty
    # AvoidPattern that made an earlier version of this test hang for
    # several CPU-minutes instead of genuinely testing anything.
    "start": 3, "end": 33, "length_base_unit": 3, "sequence": "GCT" * 10,
    "num_base_units": 10, "log10_prob_slippage_ecoli": -1.0,
}])


def test_reoptimize_until_stable_stops_as_soon_as_a_round_finds_nothing(monkeypatch):
    # Regression test for the real gap this closes: a single fixed
    # optimize->detect->optimize recipe never re-checks its own final edits.
    # Round 1 "finds" a hotspot and fixes it; round 2's detection pass
    # (faked here to isolate the looping logic from DNAChisel's actual,
    # nondeterministic behavior) reports nothing left - convergence should
    # stop there, not loop needlessly.
    seq = "ATG" + "GCT" * 10 + "TAA"
    calls = {"n": 0}

    def fake_extractor(curr_seq, *args, **kwargs):
        calls["n"] += 1
        return _fake_sites(_HOTSPOT_DF if calls["n"] == 1 else None)

    monkeypatch.setattr(pipeline_module, "suspect_site_extractor", fake_extractor)

    final_seq, obj_description, num_edits, cumulative, rounds_used = reoptimize_until_stable(
        seq, **_REOPT_KWARGS)

    assert calls["n"] == 2  # one round that found something, one confirming round that didn't
    assert rounds_used == 2
    assert obj_description is not None
    assert num_edits > 0
    assert len(cumulative['df_slippage']) == 1  # only round 1's finding, not an empty round 2 row
    assert final_seq != seq


def test_reoptimize_until_stable_respects_the_round_cap(monkeypatch):
    # A pathological "fixing one hotspot always reveals another" scenario
    # must not loop forever - it stops at max_rounds, with every round's
    # finding preserved in the cumulative history. Reports the pattern that's
    # ACTUALLY at that position in the current sequence each round (not a
    # stale, no-longer-present one) - an AvoidPattern constraint for a
    # pattern that genuinely isn't there anymore is trivially satisfiable and
    # would defeat the point of "always finds something"; worse, a stale
    # constraint DNAChisel can never actually locate is exactly the kind of
    # malformed input that made an earlier version of this test hang for
    # several CPU-minutes rather than genuinely testing the round cap.
    seq = "ATG" + "GCT" * 10 + "TAA"

    def fake_extractor(curr_seq, *args, **kwargs):
        window = curr_seq[3:33]
        return _fake_sites(pd.DataFrame([{
            "start": 3, "end": 33, "length_base_unit": 3, "sequence": window,
            "num_base_units": 10, "log10_prob_slippage_ecoli": -1.0,
        }]))

    monkeypatch.setattr(pipeline_module, "suspect_site_extractor", fake_extractor)

    _, _, _, cumulative, rounds_used = reoptimize_until_stable(seq, max_rounds=3, **_REOPT_KWARGS)

    assert rounds_used == 3
    assert len(cumulative['df_slippage']) == 3


def test_reoptimize_until_stable_leaves_a_clean_sequence_untouched():
    # A real (not faked) detection pass on an already-clean sequence should
    # find nothing in round 1 and stop immediately, with no obj_description
    # (nothing was ever re-optimized) and the sequence unchanged. Confirmed
    # directly (not just assumed) that this specific sequence has zero
    # slippage/recombination candidates before relying on that here - a
    # shorter all-A run looked "clean" by eye but is actually a real,
    # detectable length-3 repeat ("AAA" x3).
    seq = "ATG" + "GCATCGATCGTAGCTAGCTA" + "TAA"

    final_seq, obj_description, num_edits, cumulative, rounds_used = reoptimize_until_stable(
        seq, **_REOPT_KWARGS)

    assert rounds_used == 1
    assert obj_description is None
    assert num_edits == 0
    assert final_seq == seq
    assert len(cumulative['df_slippage']) == 0
    assert len(cumulative['df_recombination']) == 0


def test_reoptimize_until_stable_reports_progress_via_callback(monkeypatch):
    seq = "ATG" + "GCT" * 10 + "TAA"
    calls = {"n": 0}
    seen_rounds = []

    def fake_extractor(curr_seq, *args, **kwargs):
        calls["n"] += 1
        return _fake_sites(_HOTSPOT_DF if calls["n"] == 1 else None)

    monkeypatch.setattr(pipeline_module, "suspect_site_extractor", fake_extractor)

    reoptimize_until_stable(
        seq, on_round_start=lambda round_index, sites: seen_rounds.append(round_index), **_REOPT_KWARGS)

    assert seen_rounds == [1, 2]


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
