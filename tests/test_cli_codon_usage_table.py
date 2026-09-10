"""Tests that --codon-usage-table-file is wired correctly through
eso.cli.main: a bad file must fail fast with a friendly message (no
pipeline run), and a good file's table must reach eso.pipeline.main
unchanged.
"""

import eso.cli as cli


def _write(tmp_path, name, content):
    file_path = tmp_path / name
    file_path.write_text(content)
    return str(file_path)


def test_bad_codon_usage_table_file_fails_fast_with_friendly_message(tmp_path, capsys):
    file_path = _write(tmp_path, "bad.csv", "codon,aa\nATG,M\n")

    exit_code = cli.main(['--codon-usage-table-file', file_path])

    assert exit_code == 1
    assert "missing the column" in capsys.readouterr().err


def test_good_codon_usage_table_file_reaches_pipeline_main(tmp_path, monkeypatch):
    file_path = _write(tmp_path, "good.csv", 'codon,aa,freq_within_aa\nGCT,A,0.5\nGCC,A,0.5\n')

    captured = {}

    def fake_run_pipeline(**kwargs):
        captured.update(kwargs)
        return 'Success!', []

    monkeypatch.setattr(cli, 'run_pipeline', fake_run_pipeline)

    exit_code = cli.main(['--codon-usage-table-file', file_path])

    assert exit_code == 0
    assert captured['codon_usage_table']['A'] == {'GCT': 0.5, 'GCC': 0.5}


def test_codon_usage_table_together_with_custom_score_is_rejected(tmp_path, monkeypatch, capsys):
    # Regression test for a real gap, upgraded from a warning to a hard
    # error (matching this CLI's other mutual-exclusivity checks, e.g.
    # --codon-usage-table-file + --derive-codon-usage-table-from-assembly):
    # optimize.py ignores codon_usage_table entirely whenever custom_score_fn
    # is also given (see eso.optimize._codon_optimization_objectives) - a
    # warning alone still let a user waste a real network fetch/CPU work
    # deriving a codon-usage table that would then be silently discarded.
    # Now rejected upfront, before either side's derivation runs at all.
    table_path = _write(tmp_path, "good.csv", 'codon,aa,freq_within_aa\nGCT,A,0.5\nGCC,A,0.5\n')
    score_path = _write(tmp_path, "score.py", "def score(seq):\n    return len(seq)\n")

    monkeypatch.setattr(cli, 'run_pipeline', lambda **kwargs: ('Success!', []))

    exit_code = cli.main(['--codon-usage-table-file', table_path, '--custom-score-file', score_path])

    assert exit_code == 1
    assert "can't both be given" in capsys.readouterr().err


def test_codon_usage_table_conflict_is_rejected_before_any_network_fetch(monkeypatch, capsys):
    # Same rejection, via the --derive-*-from-assembly flags - must fail
    # before attempting either network fetch, not just before pipeline.main.
    def _unexpected_fetch(*args, **kwargs):
        raise AssertionError("should not fetch anything - rejected before either derivation runs")

    monkeypatch.setattr(cli, 'fetch_genome_package_for', _unexpected_fetch)

    exit_code = cli.main([
        '--derive-codon-usage-table-from-assembly', 'GCF_000005845.2',
        '--derive-tai-score-from-assembly', 'GCF_000005845.2', '--tai-kingdom', 'prokaryote',
    ])

    assert exit_code == 1
    assert "can't both be given" in capsys.readouterr().err


def test_no_codon_usage_table_file_leaves_it_none(monkeypatch):
    captured = {}

    def fake_run_pipeline(**kwargs):
        captured.update(kwargs)
        return 'Success!', []

    monkeypatch.setattr(cli, 'run_pipeline', fake_run_pipeline)

    exit_code = cli.main(['--organism-name', 'kompas'])

    assert exit_code == 0
    assert captured['codon_usage_table'] is None
