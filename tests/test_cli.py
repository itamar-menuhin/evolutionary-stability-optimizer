"""Tests for eso.cli's flag surface beyond --custom-score-file/--common-motifs
(already covered by tests/test_cli_custom_score.py).

Includes a regression test for a real bug: --num-sites is parsed via
argparse as `type=float` (its default is `np.inf`, not representable as an
int), but eso.detection.recombination/slippage and their staubility_variant
"fast" siblings called `df.head(num_sites)` directly whenever a finite limit
was given - pandas' DataFrame.head() rejects a float positional count and
crashes with TypeError. eso.detection.methylation already guarded this
correctly with `int(num_sites)`; the other four call sites didn't. Confirmed
directly before fixing: passing a real, non-mocked --num-sites through the
CLI crashed with "cannot do positional indexing ... of type float". Existing
end-to-end coverage (tests/test_pipeline_integration.py) never caught this
because it calls eso.pipeline.main() directly with an int, bypassing the
CLI's argparse float conversion entirely.
"""

import eso.cli as cli
from eso.pipeline import main as run_pipeline_module_main


def _fake_pipeline(capture):
    def fake_run_pipeline(**kwargs):
        capture.update(kwargs)
        return 'Success!', []
    return fake_run_pipeline


def test_num_sites_flag_reaches_pipeline_as_a_float(monkeypatch):
    captured = {}
    monkeypatch.setattr(cli, 'run_pipeline', _fake_pipeline(captured))

    exit_code = cli.main(['--num-sites', '5'])

    assert exit_code == 0
    assert captured['num_sites'] == 5.0


def test_num_sites_actually_limits_detected_sites_end_to_end(tmp_path, monkeypatch):
    # regression test for the float-vs-int .head() crash described in the
    # module docstring - exercises the real detection pipeline (not a mock)
    # with --num-sites passed as a CLI string, the way a real user would.
    monkeypatch.chdir(tmp_path)
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    output_dir = tmp_path / "output"

    # several repeats of the same 3nt unit: multiple independent slippage
    # hotspots, enough to actually test that --num-sites=1 truncates them.
    seq = "ATG" + ("GCT" * 6 + "AAA" * 6) * 3 + "TAA"
    (input_dir / "gene.fasta").write_text(f">gene\n{seq}\n")

    exit_code = cli.main([
        '--input-folder', str(input_dir),
        '--output-path', str(output_dir),
        '--num-sites', '1',
        '--no-optimize',
        '--organism-name', 'kompas',
    ])

    assert exit_code == 0
    slippage_csv = output_dir / "gene" / "slippage_sites.csv"
    if slippage_csv.exists():
        import pandas as pd
        assert len(pd.read_csv(slippage_csv)) <= 1


def test_no_optimize_flag_skips_optimization(monkeypatch):
    captured = {}
    monkeypatch.setattr(cli, 'run_pipeline', _fake_pipeline(captured))

    exit_code = cli.main(['--no-optimize'])

    assert exit_code == 0
    assert captured['optimize'] is False


def test_optimize_is_the_default(monkeypatch):
    captured = {}
    monkeypatch.setattr(cli, 'run_pipeline', _fake_pipeline(captured))

    exit_code = cli.main([])

    assert exit_code == 0
    assert captured['optimize'] is True


def test_mini_maxi_gc_flags_reach_pipeline(monkeypatch):
    captured = {}
    monkeypatch.setattr(cli, 'run_pipeline', _fake_pipeline(captured))

    exit_code = cli.main(['--mini-gc', '0.2', '--maxi-gc', '0.6'])

    assert exit_code == 0
    assert captured['mini_gc'] == 0.2
    assert captured['maxi_gc'] == 0.6


def test_recombination_and_slippage_mode_flags_reach_pipeline(monkeypatch):
    captured = {}
    monkeypatch.setattr(cli, 'run_pipeline', _fake_pipeline(captured))

    exit_code = cli.main(['--recombination-mode', 'fast', '--slippage-mode', 'fast'])

    assert exit_code == 0
    assert captured['recombination_mode'] == 'fast'
    assert captured['slippage_mode'] == 'fast'


def test_unknown_mode_choice_is_rejected_by_argparse(capsys):
    try:
        cli.main(['--recombination-mode', 'not_a_real_mode'])
        assert False, "argparse should have exited on an invalid choice"
    except SystemExit as e:
        assert e.code == 2
        assert "invalid choice" in capsys.readouterr().err


def test_pipeline_failure_message_is_reported_and_exit_code_is_nonzero(monkeypatch):
    def fake_run_pipeline(**kwargs):
        return "The minimal GC content must be less than the maximum!", []
    monkeypatch.setattr(cli, 'run_pipeline', fake_run_pipeline)

    exit_code = cli.main(['--mini-gc', '0.8', '--maxi-gc', '0.2'])

    assert exit_code == 1


def test_organism_and_method_flags_reach_pipeline(monkeypatch):
    captured = {}
    monkeypatch.setattr(cli, 'run_pipeline', _fake_pipeline(captured))

    exit_code = cli.main(['--organism-name', 'kompas', '--method', 'match_codon_usage'])

    assert exit_code == 0
    assert captured['organism_name'] == 'kompas'
    assert captured['method'] == 'match_codon_usage'


def test_indexes_file_flag_reaches_pipeline_as_the_expected_dict(tmp_path, monkeypatch):
    captured = {}
    monkeypatch.setattr(cli, 'run_pipeline', _fake_pipeline(captured))

    indexes_path = tmp_path / "indexes.json"
    indexes_path.write_text(
        '[{"file": "my_gene", "seq_index": "0", "orf_regions": "1-6, 51-68", '
        '"exclusion_regions": "1-6, 50-68"}]'
    )

    exit_code = cli.main(['--indexes-file', str(indexes_path)])

    assert exit_code == 0
    assert captured['indexes'] == {("my_gene", "0"): ("1-6, 51-68", "1-6, 50-68")}


def test_no_indexes_file_flag_leaves_it_none(monkeypatch):
    captured = {}
    monkeypatch.setattr(cli, 'run_pipeline', _fake_pipeline(captured))

    exit_code = cli.main([])

    assert exit_code == 0
    assert captured['indexes'] is None


def test_bad_indexes_file_fails_fast_with_friendly_message(tmp_path, capsys):
    indexes_path = tmp_path / "indexes.json"
    indexes_path.write_text('{not valid json')

    exit_code = cli.main(['--indexes-file', str(indexes_path)])

    assert exit_code == 1
    assert "isn't valid JSON" in capsys.readouterr().err


def test_indexes_file_actually_locks_the_exclusion_region_end_to_end(tmp_path, monkeypatch):
    # exercises the real pipeline (not a mock): an exclusion region declared
    # via --indexes-file must genuinely survive optimization untouched.
    monkeypatch.chdir(tmp_path)
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    output_dir = tmp_path / "output"

    seq = "ATG" + "TTT" * 10 + "TAA"
    (input_dir / "gene.fasta").write_text(f">gene\n{seq}\n")

    indexes_path = tmp_path / "indexes.json"
    indexes_path.write_text(
        '[{"file": "gene", "seq_index": "0", "orf_regions": "1-33", "exclusion_regions": "4-33"}]'
    )

    exit_code = cli.main([
        '--input-folder', str(input_dir),
        '--output-path', str(output_dir),
        '--indexes-file', str(indexes_path),
        '--organism-name', 'kompas',
    ])

    assert exit_code == 0
    final_sequence_path = output_dir / "gene" / "final_sequence.txt"
    contents = final_sequence_path.read_text()
    final_seq = ''.join(
        line.strip() for line in contents.splitlines()[contents.splitlines().index("The final sequence is:") + 1:]
    )
    assert final_seq[3:33] == seq[3:33]


def test_cli_main_module_actually_imports_pipeline_main():
    # sanity check that eso.cli.run_pipeline is really eso.pipeline.main,
    # not a stale copy - the monkeypatch-based tests above would pass even
    # if cli.py called a dead reference, since they replace cli.run_pipeline
    # directly rather than eso.pipeline.main.
    assert cli.run_pipeline is run_pipeline_module_main
