"""Tests for eso.custom_score.CustomScore - the user-supplied-score DNAChisel
objective, as an alternative to CodonOptimize (CAI/tAI).

CustomScore previously also supported a "windowed" mode (score_fn called on
fixed-size chunks, summed) mirroring how the built-in CAI/tAI codon-usage
scoring works. It was removed after benchmarking found no case where it was
actually faster than whole-sequence evaluation (comparable at best,
meaningfully slower at worst - see docs/detector-comparisons.md for the
full investigation, including an initial benchmark that was itself found to
be flawed and corrected), while carrying a real, unpreventable correctness
risk (a score_fn that doesn't genuinely decompose per-chunk would silently
compute something structurally unrelated to the intended score).
"""

import warnings

import dnachisel
import pytest

from eso.custom_score import CustomScore
from eso.optimize import optimization_engine


def _gc_count(seq):
    return sum(1 for nt in seq if nt in "GC")


def test_scores_whole_sequence_and_warns():
    with pytest.warns(UserWarning, match="full scored region"):
        spec = CustomScore(_gc_count)

    problem = dnachisel.DnaOptimizationProblem(sequence="ATGCATGCAT", objectives=[spec])
    spec = spec.initialized_on_problem(problem)
    evaluation = spec.evaluate(problem)
    assert evaluation.score == _gc_count("ATGCATGCAT")


def test_localized_returns_self_unchanged():
    # The score can't be restricted to a sub-region without changing its
    # meaning, so localized() must always return the same (whole-region)
    # specification, never None and never a narrowed location.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        spec = CustomScore(_gc_count)

    problem = dnachisel.DnaOptimizationProblem(sequence="ATGCATGCAT" * 5, objectives=[spec])
    spec = spec.initialized_on_problem(problem)
    localized = spec.localized(dnachisel.Location(5, 8), problem=problem)
    assert localized.location == spec.location


def test_minimize_negates_score():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        maximize_spec = CustomScore(_gc_count)
        minimize_spec = CustomScore(_gc_count, minimize=True)

    seq = "ATGCATGCAT"
    problem = dnachisel.DnaOptimizationProblem(sequence=seq, objectives=[maximize_spec])
    maximize_spec = maximize_spec.initialized_on_problem(problem)
    minimize_spec = minimize_spec.initialized_on_problem(problem)

    assert minimize_spec.evaluate(problem).score == -maximize_spec.evaluate(problem).score


def test_location_restricts_scoring_to_a_sub_region():
    seq = "AAA" + "GGG" + "AAA"  # only the middle 3nt should be scored
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        spec = CustomScore(_gc_count, location=(3, 6))

    problem = dnachisel.DnaOptimizationProblem(sequence=seq, objectives=[spec])
    spec = spec.initialized_on_problem(problem)
    assert spec.evaluate(problem).score == _gc_count("GGG")


def test_optimization_engine_accepts_custom_score_and_moves_toward_higher_gc():
    # End-to-end smoke test: optimization_engine with custom_score_fn instead
    # of organism_name/CAI, on a sequence with room to improve GC content
    # while staying in-frame (synonymous codon choices only).
    seq = "ATG" + "AAA" * 20 + "TAA"  # all-AT codons except start/stop

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        final_seq, _, _ = optimization_engine(
            seq,
            custom_score_fn=_gc_count,
        )

    assert _gc_count(final_seq) > _gc_count(seq)
