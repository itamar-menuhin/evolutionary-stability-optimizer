"""Tests for eso.custom_score.CustomScore - the user-supplied-score DNAChisel
objective, as an alternative to CodonOptimize (CAI/tAI).
"""

import warnings

import dnachisel
import pytest

from eso.custom_score import CustomScore
from eso.optimize import optimization_engine


def _gc_count(seq):
    return sum(1 for nt in seq if nt in "GC")


def test_global_mode_scores_whole_sequence_and_warns():
    with pytest.warns(UserWarning, match="global mode"):
        spec = CustomScore(_gc_count)

    problem = dnachisel.DnaOptimizationProblem(sequence="ATGCATGCAT", objectives=[spec])
    spec = spec.initialized_on_problem(problem)
    evaluation = spec.evaluate(problem)
    assert evaluation.score == _gc_count("ATGCATGCAT")


def test_global_mode_localized_returns_self_unchanged():
    # A global score can't be restricted to a sub-region without changing
    # its meaning, so localized() must always return the same (whole-sequence)
    # specification, never None and never a narrowed location.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        spec = CustomScore(_gc_count)

    problem = dnachisel.DnaOptimizationProblem(sequence="ATGCATGCAT" * 5, objectives=[spec])
    spec = spec.initialized_on_problem(problem)
    localized = spec.localized(dnachisel.Location(5, 8), problem=problem)
    assert localized.location == spec.location


def test_windowed_mode_sums_per_window_scores_no_warning_at_construction():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        spec = CustomScore(_gc_count, window=3)

    # 10nt -> 3 full windows of 3, last nt dropped - initialized_on_problem (not
    # construction) is where this now warns, see test below.
    seq = "ATGCATGCAT"
    problem = dnachisel.DnaOptimizationProblem(sequence=seq, objectives=[spec])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        spec = spec.initialized_on_problem(problem)
    evaluation = spec.evaluate(problem)
    expected = sum(_gc_count(seq[i:i + 3]) for i in range(0, 9, 3))
    assert evaluation.score == expected


def test_windowed_mode_warns_when_region_length_not_a_multiple_of_window():
    # regression test: a scored region whose length isn't a multiple of `window`
    # silently drops its trailing remainder from every score_fn call (only whole
    # windows are ever scored) - confirmed directly this was previously silent.
    # Must now warn so a user isn't surprised by nucleotides never reaching their
    # score_fn.
    spec = CustomScore(_gc_count, window=3)
    seq = "ATGCATGCAT"  # 10nt, remainder 1
    problem = dnachisel.DnaOptimizationProblem(sequence=seq, objectives=[spec])

    with pytest.warns(UserWarning, match="not a multiple of window"):
        spec.initialized_on_problem(problem)


def test_windowed_mode_no_warning_when_region_length_is_a_multiple_of_window():
    spec = CustomScore(_gc_count, window=3)
    seq = "ATGCATGCA" * 3  # 27nt, a clean multiple of 3
    problem = dnachisel.DnaOptimizationProblem(sequence=seq, objectives=[spec])

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        spec.initialized_on_problem(problem)


def test_windowed_mode_localized_restricts_to_nearby_windows():
    spec = CustomScore(_gc_count, window=3)
    seq = "A" * 300
    problem = dnachisel.DnaOptimizationProblem(sequence=seq, objectives=[spec])
    spec = spec.initialized_on_problem(problem)

    localized = spec.localized(dnachisel.Location(150, 153), problem=problem)
    assert localized is not None
    # restricted to near the mutation, not the whole 300nt sequence
    assert localized.location.end - localized.location.start < 300


def test_windowed_mode_localized_returns_none_for_disjoint_region():
    spec = CustomScore(_gc_count, window=3, location=dnachisel.Location(0, 30))
    seq = "A" * 300
    problem = dnachisel.DnaOptimizationProblem(sequence=seq, objectives=[spec])
    spec = spec.initialized_on_problem(problem)

    localized = spec.localized(dnachisel.Location(200, 203), problem=problem)
    assert localized is None


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
            custom_score_window=3,
        )

    assert _gc_count(final_seq) > _gc_count(seq)
