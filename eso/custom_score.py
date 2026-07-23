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

    `score_fn` is called once on the whole scored region (the full sequence,
    or `location`, if given) on every trial mutation during `optimize()`.
    This is always correct, with no assumption about how the score behaves.

    An earlier version of this class also supported a "windowed" mode
    (calling score_fn per fixed-size chunk and summing), matching how the
    built-in CAI/tAI codon-usage scoring works. It was removed: benchmarking
    found no case where it was actually faster than this whole-sequence
    approach (comparable at best, meaningfully slower at worst, since
    DNAChisel's own optimizer ends up calling score_fn considerably more
    often when the objective is chunk-localizable) - see
    docs/detector-comparisons.md for the full investigation, including an
    initial, incorrect benchmark that was itself corrected. Windowed mode
    also carried a real, unpreventable correctness risk: a user-supplied
    score_fn that doesn't genuinely decompose as a sum over independent
    chunks (true of most real external/ML models) would silently compute a
    different, structurally unrelated quantity, with no reliable way to
    detect this automatically. Given it offered no confirmed benefit and a
    real risk, it was removed rather than kept as an unverified "maybe
    faster sometimes" option.

    Parameters
    ----------
    score_fn
        Callable taking a DNA sequence (str) and returning a float, higher is
        better (pass `minimize=True` if lower is better).
    location
        Restrict the objective to a sub-region of the full sequence. Defaults
        to the whole sequence.
    minimize
        If True, `score_fn`'s return value is negated before use, so that a
        *lower* raw score is treated as better.
    """

    best_possible_score = None

    def __init__(self, score_fn, location=None, minimize=False, boost=1.0):
        self.score_fn = score_fn
        self.location = dnachisel.Location.from_data(location)
        self.minimize = minimize
        self.boost = boost
        warnings.warn(
            "CustomScore re-evaluates score_fn on the full scored region for "
            "every trial mutation during optimize(), which can be slow for "
            "an expensive score_fn or a long sequence.",
            stacklevel=2,
        )

    def initialized_on_problem(self, problem, role=None):
        return self._copy_with_full_span_if_no_location(problem)

    def _score_sequence(self, sequence):
        score = self.score_fn(sequence)
        return -score if self.minimize else score

    def evaluate(self, problem):
        sequence = self.location.extract_sequence(problem.sequence)
        score = self._score_sequence(sequence)
        return dnachisel.SpecEvaluation(self, problem, score, locations=[self.location])

    def localized(self, location, problem=None):
        # The score can't be restricted to a sub-region without changing its
        # meaning, so always re-evaluate the whole thing.
        return self

    def label_parameters(self):
        params = []
        if self.minimize:
            params.append(("minimize", "True"))
        return params

    def short_label(self):
        return "custom score"


class CustomScoreFileError(Exception):
    """A custom-score file failed to load or didn't behave as expected.

    Raised with a plain-English message aimed at someone who wrote a scoring
    function but doesn't necessarily know Python packaging or DNAChisel -
    the goal is that this exception's message alone is enough to fix the
    problem, without needing to read a traceback through this codebase.
    """


def load_custom_score_from_file(file_path, function_name='score'):
    """Load a user-authored scoring function from a plain Python file.

    The file is expected to define a function named `function_name`
    (default: `score`) taking a DNA sequence string and returning a number
    (higher = better).

    This is the mechanism behind `eso-optimize --custom-score-file`, and is
    also usable directly from Python. Validates the file eagerly (missing
    function, wrong type, or an error raised on a short test sequence) and
    raises `CustomScoreFileError` with a message meant to be read and acted
    on directly by whoever wrote the scoring file - not a Python expert.

    Returns
    -------
    score_fn
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

    try:
        result = score_fn(_VALIDATION_TEST_SEQUENCE)
    except Exception as e:
        raise CustomScoreFileError(
            f"Your `{function_name}` function raised an error when tested on the "
            f"sequence '{_VALIDATION_TEST_SEQUENCE}': {e!r}. Please fix `{function_name}` in "
            f"'{file_path}' and try again."
        ) from e

    if not isinstance(result, (int, float)):
        raise CustomScoreFileError(
            f"Your `{function_name}` function returned a {type(result).__name__} "
            f"({result!r}) instead of a number, when tested on '{_VALIDATION_TEST_SEQUENCE}'. "
            f"Make sure it ends with `return <a number>`."
        )

    return score_fn
