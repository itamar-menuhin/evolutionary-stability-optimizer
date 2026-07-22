"""Tests for eso.custom_score.load_custom_score_from_file - the file-based
loader behind `eso-optimize --custom-score-file`, aimed at non-programmer
users authoring a scoring function from examples/custom_score_template.py.
"""

from os import path

import pytest

from eso.custom_score import CustomScoreFileError, load_custom_score_from_file


def _write(tmp_path, name, content):
    file_path = tmp_path / name
    file_path.write_text(content)
    return str(file_path)


def test_loads_the_actual_template_file():
    template = path.join(path.dirname(__file__), '..', 'examples', 'custom_score_template.py')
    score_fn, window = load_custom_score_from_file(template)
    assert window == 3
    assert score_fn("GCG") == 3
    assert score_fn("ATA") == 0


def test_missing_file_gives_friendly_message():
    with pytest.raises(CustomScoreFileError, match="Can't find"):
        load_custom_score_from_file("no_such_file_here.py")


def test_file_with_syntax_error_gives_friendly_message(tmp_path):
    file_path = _write(tmp_path, "bad.py", "def score(seq)\n    return 1\n")
    with pytest.raises(CustomScoreFileError, match="error in it"):
        load_custom_score_from_file(file_path)


def test_missing_function_gives_friendly_message(tmp_path):
    file_path = _write(tmp_path, "no_score.py", "WINDOW = 3\ndef banana(seq):\n    return 1\n")
    with pytest.raises(CustomScoreFileError, match="doesn't define a function called `score`"):
        load_custom_score_from_file(file_path)


def test_custom_function_name_is_respected(tmp_path):
    file_path = _write(tmp_path, "renamed.py", "WINDOW = 3\ndef my_gc_score(seq):\n    return seq.count('G')\n")
    score_fn, window = load_custom_score_from_file(file_path, function_name='my_gc_score')
    assert window == 3
    assert score_fn("GGG") == 3


def test_non_callable_score_gives_friendly_message(tmp_path):
    file_path = _write(tmp_path, "not_callable.py", "score = 42\n")
    with pytest.raises(CustomScoreFileError, match="isn't a function"):
        load_custom_score_from_file(file_path)


def test_score_function_raising_an_error_is_caught_with_friendly_message(tmp_path):
    file_path = _write(tmp_path, "raises.py", "WINDOW = 3\ndef score(seq):\n    return 1 / 0\n")
    with pytest.raises(CustomScoreFileError, match="raised an error when tested"):
        load_custom_score_from_file(file_path)


def test_score_function_returning_non_number_gives_friendly_message(tmp_path):
    file_path = _write(tmp_path, "wrong_type.py", "WINDOW = 3\ndef score(seq):\n    return 'not a number'\n")
    with pytest.raises(CustomScoreFileError, match="instead of a number"):
        load_custom_score_from_file(file_path)


def test_invalid_window_value_gives_friendly_message(tmp_path):
    file_path = _write(tmp_path, "bad_window.py", "WINDOW = -1\ndef score(seq):\n    return 1\n")
    with pytest.raises(CustomScoreFileError, match="positive whole number"):
        load_custom_score_from_file(file_path)


def test_missing_window_defaults_to_none_global_mode(tmp_path):
    file_path = _write(tmp_path, "no_window.py", "def score(seq):\n    return len(seq)\n")
    score_fn, window = load_custom_score_from_file(file_path)
    assert window is None
    assert score_fn("ATGATG") == 6
