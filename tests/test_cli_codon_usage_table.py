"""Tests that --codon-usage-table-file is wired correctly through
eso.cli.main: a bad file must fail fast with a friendly message (no
pipeline run), and a good file's table must reach eso.pipeline.main
unchanged.
"""

import pytest

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


def test_codon_usage_table_together_with_custom_score_warns_it_will_be_ignored(tmp_path, monkeypatch):
    # Regression test for a real, silent-waste gap: optimize.py ignores
    # codon_usage_table entirely whenever custom_score_fn is also given (see
    # eso.optimize._codon_optimization_objectives) - a user deriving BOTH
    # used to get no indication the codon-usage table's real fetch/CPU work
    # was thrown away. Confirmed directly before this fix.
    table_path = _write(tmp_path, "good.csv", 'codon,aa,freq_within_aa\nGCT,A,0.5\nGCC,A,0.5\n')
    score_path = _write(tmp_path, "score.py", "def score(seq):\n    return len(seq)\n")

    monkeypatch.setattr(cli, 'run_pipeline', lambda **kwargs: ('Success!', []))

    with pytest.warns(UserWarning, match="codon-usage table will be ignored"):
        exit_code = cli.main(['--codon-usage-table-file', table_path, '--custom-score-file', score_path])

    assert exit_code == 0


def test_no_codon_usage_table_file_leaves_it_none(monkeypatch):
    captured = {}

    def fake_run_pipeline(**kwargs):
        captured.update(kwargs)
        return 'Success!', []

    monkeypatch.setattr(cli, 'run_pipeline', fake_run_pipeline)

    exit_code = cli.main(['--organism-name', 'kompas'])

    assert exit_code == 0
    assert captured['codon_usage_table'] is None
