"""Tests that --custom-score-file is wired correctly through eso.cli.main:
a bad file must fail fast with a friendly message (no pipeline run), and a
good file's function/window must reach eso.pipeline.main unchanged.
"""

import eso.cli as cli


def _write(tmp_path, name, content):
    file_path = tmp_path / name
    file_path.write_text(content)
    return str(file_path)


def test_bad_custom_score_file_fails_fast_with_friendly_message(tmp_path, capsys):
    file_path = _write(tmp_path, "bad.py", "def not_score(seq):\n    return 1\n")

    exit_code = cli.main(['--custom-score-file', file_path])

    assert exit_code == 1
    assert "doesn't define a function called `score`" in capsys.readouterr().err


def test_good_custom_score_file_reaches_pipeline_main(tmp_path, monkeypatch):
    file_path = _write(tmp_path, "good.py", "WINDOW = 3\ndef score(seq):\n    return seq.count('G')\n")

    captured = {}

    def fake_run_pipeline(**kwargs):
        captured.update(kwargs)
        return 'Success!', []

    monkeypatch.setattr(cli, 'run_pipeline', fake_run_pipeline)

    exit_code = cli.main(['--custom-score-file', file_path, '--custom-score-minimize'])

    assert exit_code == 0
    assert captured['custom_score_window'] == 3
    assert captured['custom_score_minimize'] is True
    assert captured['custom_score_fn']("GGG") == 3


def test_no_custom_score_file_leaves_it_none(tmp_path, monkeypatch):
    captured = {}

    def fake_run_pipeline(**kwargs):
        captured.update(kwargs)
        return 'Success!', []

    monkeypatch.setattr(cli, 'run_pipeline', fake_run_pipeline)

    exit_code = cli.main(['--organism-name', 'not_specified'])

    assert exit_code == 0
    assert captured['custom_score_fn'] is None


def test_common_motifs_flag_parses_to_a_list(monkeypatch):
    captured = {}

    def fake_run_pipeline(**kwargs):
        captured.update(kwargs)
        return 'Success!', []

    monkeypatch.setattr(cli, 'run_pipeline', fake_run_pipeline)

    exit_code = cli.main(['--compute-motifs', '--common-motifs', 'dam, dcm'])

    assert exit_code == 0
    assert captured['common_motifs'] == ['dam', 'dcm']


def test_no_common_motifs_flag_leaves_it_none(monkeypatch):
    captured = {}

    def fake_run_pipeline(**kwargs):
        captured.update(kwargs)
        return 'Success!', []

    monkeypatch.setattr(cli, 'run_pipeline', fake_run_pipeline)

    exit_code = cli.main([])

    assert exit_code == 0
    assert captured['common_motifs'] is None
