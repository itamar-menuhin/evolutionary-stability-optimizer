"""Wrap an arbitrary user-supplied scoring function as a DNAChisel objective,
as an alternative to eso.optimize's built-in CodonOptimize (CAI/tAI-style)
objective.
"""

import importlib.util
import warnings
from os import path

import dnachisel

#: DNA alphabet used to build the dummy sequence that a freshly-loaded custom
#: score file is test-run against, so authoring mistakes surface immediately
#: with a plain-English message instead of mid-optimization, deep inside
#: DNAChisel's internals.
_VALIDATION_TEST_SEQUENCE = "ATGATGATGATGATGATGATGATGATGATG"  # 30nt, multiple of 3


class CustomScore(dnachisel.Specification):
    """DNAChisel objective that maximizes an arbitrary Python function of the
    sequence, instead of a codon-usage table.

    Two modes, chosen by whether `window` is given:

    - **windowed** (`window=N`): `score_fn` is called on each successive
      non-overlapping N-nt chunk of the sequence and the results are summed.
      Only correct if your true score genuinely decomposes as a sum over
      fixed-size windows (this is what CAI/tAI-style per-codon scoring does -
      e.g. N=3 for a per-codon function). In exchange, this is fast: DNAChisel
      only re-evaluates the windows touched by a given trial mutation
      (see `localized()` below), not the whole sequence, on every step.

    - **global** (`window=None`, the default): `score_fn` is called once on
      the whole sequence (or `location`, if given). Always correct - no
      assumption about how the score decomposes - but every trial mutation
      during `problem.optimize()` re-evaluates the entire sequence from
      scratch, since a global score can't be localized. This can be very
      slow for long sequences or an expensive `score_fn`; a warning is
      raised at construction time to make this cost explicit up front.

    Parameters
    ----------
    score_fn
        Callable taking a DNA sequence (str) and returning a float, higher is
        better (pass `minimize=True` if lower is better).
    window
        If given, enables windowed mode (see above) with this window size, in
        nucleotides. If None (default), global mode is used.
    location
        Restrict the objective to a sub-region of the full sequence. Defaults
        to the whole sequence.
    minimize
        If True, `score_fn`'s return value is negated before use, so that a
        *lower* raw score is treated as better.
    """

    best_possible_score = None

    def __init__(self, score_fn, window=None, location=None, minimize=False, boost=1.0):
        self.score_fn = score_fn
        self.window = window
        self.location = dnachisel.Location.from_data(location)
        self.minimize = minimize
        self.boost = boost
        if window is None:
            warnings.warn(
                "CustomScore is running in global mode (no `window` given): "
                "score_fn will be re-evaluated on the full sequence for "
                "every trial mutation during optimize(), which can be very "
                "slow on long sequences or an expensive score_fn. If your "
                "score can be computed as a sum over fixed-size windows "
                "(like per-codon scoring), pass `window=N` for a much "
                "faster, localized evaluation.",
                stacklevel=2,
            )

    def initialized_on_problem(self, problem, role=None):
        return self._copy_with_full_span_if_no_location(problem)

    def _score_sequence(self, sequence):
        if self.window is None:
            score = self.score_fn(sequence)
        else:
            score = sum(
                self.score_fn(sequence[start:start + self.window])
                for start in range(0, len(sequence) - self.window + 1, self.window)
            )
        return -score if self.minimize else score

    def evaluate(self, problem):
        sequence = self.location.extract_sequence(problem.sequence)
        score = self._score_sequence(sequence)
        return dnachisel.SpecEvaluation(self, problem, score, locations=[self.location])

    def localized(self, location, problem=None):
        if self.window is None:
            # A global score can't be restricted to a sub-region without
            # changing its meaning - always re-evaluate the whole thing.
            return self
        extended_location = location.extended(self.window - 1)
        new_location = self.location.overlap_region(extended_location)
        if new_location is None:
            return None
        return self.copy_with_changes(location=new_location)

    def label_parameters(self):
        params = [("window", str(self.window) if self.window else "global")]
        if self.minimize:
            params.append(("minimize", "True"))
        return params

    def short_label(self):
        return "custom score (%s)" % ("global" if self.window is None else f"window={self.window}")


class CustomScoreFileError(Exception):
    """A custom-score file failed to load or didn't behave as expected.

    Raised with a plain-English message aimed at someone who wrote a scoring
    function but doesn't necessarily know Python packaging or DNAChisel -
    the goal is that this exception's message alone is enough to fix the
    problem, without needing to read a traceback through this codebase.
    """


def load_custom_score_from_file(file_path, function_name='score', window_variable='WINDOW'):
    """Load a user-authored scoring function from a plain Python file.

    The file is expected to define:
    - a function named `function_name` (default: `score`) taking a DNA
      sequence string and returning a number (higher = better).
    - optionally, a module-level variable named `window_variable`
      (default: `WINDOW`) set to either an integer (the function will be
      called on each successive chunk of that many nucleotides and the
      results summed - fast, use this if your score is naturally per-codon
      or per-window) or `None` (the function is called once on the whole
      sequence each time - always correct, but much slower to optimize).
      If the file doesn't define `WINDOW` at all, `None` (whole-sequence
      mode) is assumed.

    This is the mechanism behind `eso-optimize --custom-score-file`, and is
    also usable directly from Python. Validates the file eagerly (missing
    function, wrong type, or an error raised on a short test sequence) and
    raises `CustomScoreFileError` with a message meant to be read and acted
    on directly by whoever wrote the scoring file - not a Python expert.

    Returns
    -------
    (score_fn, window)
    """
    if not path.isfile(file_path):
        raise CustomScoreFileError(
            f"Can't find the custom score file '{file_path}'. Check the path is correct.")

    module_name = f"eso_custom_score_{abs(hash(file_path))}"
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
    except Exception as e:
        raise CustomScoreFileError(
            f"'{file_path}' could not be loaded - it has an error in it. "
            f"The error was: {e!r}. Open the file and fix that error, then try again."
        ) from e

    if not hasattr(module, function_name):
        raise CustomScoreFileError(
            f"'{file_path}' doesn't define a function called `{function_name}`. "
            f"Add a function like:\n\n    def {function_name}(seq):\n        "
            f"return ...  # a number, higher = better\n\nat the top level of the file."
        )

    score_fn = getattr(module, function_name)
    if not callable(score_fn):
        raise CustomScoreFileError(
            f"'{function_name}' in '{file_path}' isn't a function - it's a {type(score_fn).__name__}. "
            f"It needs to be defined with `def {function_name}(seq): ...`."
        )

    window = getattr(module, window_variable, None)
    if window is not None and (not isinstance(window, int) or window <= 0):
        raise CustomScoreFileError(
            f"'{window_variable}' in '{file_path}' should be a positive whole number "
            f"(like 3) or `None`, not {window!r}.")

    test_sequence = _VALIDATION_TEST_SEQUENCE if window is None else (
        _VALIDATION_TEST_SEQUENCE[:window] if window <= len(_VALIDATION_TEST_SEQUENCE)
        else _VALIDATION_TEST_SEQUENCE * (window // len(_VALIDATION_TEST_SEQUENCE) + 1)
    )
    try:
        result = score_fn(test_sequence)
    except Exception as e:
        raise CustomScoreFileError(
            f"Your `{function_name}` function raised an error when tested on the "
            f"sequence '{test_sequence}': {e!r}. Please fix `{function_name}` in "
            f"'{file_path}' and try again."
        ) from e

    if not isinstance(result, (int, float)):
        raise CustomScoreFileError(
            f"Your `{function_name}` function returned a {type(result).__name__} "
            f"({result!r}) instead of a number, when tested on '{test_sequence}'. "
            f"Make sure it ends with `return <a number>`."
        )

    return score_fn, window
