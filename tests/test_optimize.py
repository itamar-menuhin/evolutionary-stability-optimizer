"""Tests for eso.optimize.optimization_engine - previously untested beyond
one pipeline smoke test and the custom-score-specific tests.

Includes a regression test for a real, previously-undiscovered bug: motif
avoidance constraints (from df_motifs) never actually did anything, because
eso.detection.methylation's `end_index` is inclusive (the index of the
motif's last nucleotide) but everything downstream of optimize.py's df_mot
construction - exclusion_site_correcter, convert_df_to_constraints, and
DNAChisel's own Location - expects exclusive end (matching Python slicing).
Every motif's AvoidPattern location was exactly one nucleotide too short to
ever contain its own pattern, so DNAChisel could never find a match there and
always reported the constraint as trivially satisfied. Confirmed directly
before fixing: a GATC motif passed via df_motifs survived optimization
completely untouched (0 edits, exact original sequence at that position).
"""

import warnings

import dnachisel
import pandas as pd
import pytest
from dnachisel.biotools import reverse_complement

from dnachisel.DnaOptimizationProblem import NoSolutionError

from eso.optimize import _constraints_to_drop, optimization_engine
from eso.sequence_utils import InvalidSequenceError


def test_ambiguity_code_raises_a_clear_error_instead_of_crashing_dnachisel():
    # Regression test for a real gap: nothing validated the sequence
    # alphabet before this, so an IUPAC ambiguity code (here "N") crashed
    # deep inside DNAChisel with an unhandled TranslationError - confirmed
    # directly before this fix. Now raises a clear, catchable error instead.
    seq = "ATG" + "NNN" * 3 + "TAA"
    with pytest.raises(InvalidSequenceError, match=r"position 4 \('N'\)"):
        optimization_engine(seq, organism_name="not_specified")


def test_empty_sequence_raises_a_clear_error_instead_of_crashing_dnachisel():
    # Regression test for a real bug: an empty sequence produced a
    # zero-length default ORF region, which crashed with a bare
    # `IndexError: string index out of range` deep inside DNAChisel's
    # EnforceTranslation - confirmed directly before this fix.
    with pytest.raises(InvalidSequenceError, match="is empty"):
        optimization_engine("", organism_name="not_specified")


def test_translation_is_preserved():
    seq = "ATG" + "TTT" * 10 + "TAA"  # all-Phe, a real slippage/CAI target
    final_seq, _, _ = optimization_engine(seq, organism_name="kompas")

    assert len(final_seq) == len(seq)
    assert final_seq[:3] == "ATG"
    assert final_seq[-3:] in ("TAA", "TAG", "TGA")


def test_last_codon_is_translation_protected_when_length_is_a_multiple_of_3():
    # Regression test for a real bug: the default orf_regions computation was
    # `((len(seq) - 1) // 3) * 3`, which for any sequence whose length is a
    # clean multiple of 3 (the normal case - a real ORF from start to stop)
    # excluded the final codon from EnforceTranslation entirely. Force an
    # AvoidPattern constraint that requires the last (stop) codon itself to
    # change: without EnforceTranslation covering it, DNAChisel is free to
    # pick ANY replacement satisfying "not TAA" (e.g. "TAC", not a stop codon
    # at all - a silently corrupted stop). With the fix, EnforceTranslation
    # forces a synonymous stop (TAG or TGA), never a non-stop codon.
    seq = "ATG" + "ATC" * 4 + "TAA"  # len=18, a multiple of 3
    assert len(seq) % 3 == 0

    df_motifs = pd.DataFrame([{
        "start_index": 15, "end_index": 17, "matching_motif": "stop_codon",
        "PSSM_score": 5.0, "actual_site": "TAA", "actual_site_reverse_conjugate": "TAA",
    }])

    final_seq, _, num_edits = optimization_engine(seq, df_motifs=df_motifs, organism_name="not_specified")

    assert len(final_seq) == len(seq)
    assert num_edits > 0
    assert final_seq[-3:] != "TAA"
    assert final_seq[-3:] in ("TAG", "TGA")


def test_gc_content_is_enforced():
    seq = "ATG" + "AAA" * 20 + "TAA"  # far below the default 30% GC floor
    final_seq, _, _ = optimization_engine(seq, mini_gc=0.3, maxi_gc=0.7, window_size_gc=20)

    gc_fraction = (final_seq.count("G") + final_seq.count("C")) / len(final_seq)
    assert gc_fraction >= 0.3


def test_exclusion_regions_are_never_modified():
    seq = "ATG" + "TTT" * 10 + "TAA"
    # lock everything except the last 2 codons before the stop
    exclusion_regions = [(3, 27)]

    final_seq, _, _ = optimization_engine(
        seq, organism_name="kompas", exclusion_regions=exclusion_regions)

    assert final_seq[3:27] == seq[3:27]


def test_unknown_organism_skips_codon_optimization_without_crashing():
    # Regression test for a real gap: an unrecognized organism_name silently
    # produced the exact same "no codon optimization happened" outcome as
    # the documented not_specified sentinel, printed (not warned) with no
    # way for a caller to tell "you opted out" apart from "your organism
    # name was a typo" - now a real, catchable UserWarning specifically for
    # the genuinely-unrecognized case.
    seq = "ATG" + "TTT" * 5 + "TAA"

    with pytest.warns(UserWarning, match="not_a_real_organism_xyz.*isn't a recognized"):
        final_seq, obj_description, _ = optimization_engine(seq, organism_name="not_a_real_organism_xyz")

    assert len(final_seq) == len(seq)
    assert final_seq[:3] == "ATG"


def test_not_specified_organism_skips_codon_optimization_silently():
    # The documented "skip codon optimization" sentinel is normal, expected
    # behavior, not a caller mistake - must NOT warn (unlike the genuinely-
    # unrecognized-organism case above), since every call using the default
    # organism_name would otherwise spam a warning for entirely intended
    # behavior.
    seq = "ATG" + "TTT" * 5 + "TAA"

    with warnings.catch_warnings():
        warnings.simplefilter("error")  # any UserWarning here fails the test
        final_seq, _, _ = optimization_engine(seq, organism_name="not_specified")

    assert len(final_seq) == len(seq)


def test_codon_usage_table_takes_precedence_over_organism_name():
    from eso.codon_usage import cub_c1

    table = cub_c1()
    # Force a single, deterministic "best" codon for Leucine - everything
    # else in this table is a real, working table (cub_c1), so this only
    # tests that the explicit table (not organism_name="kompas", a
    # deliberately different table) is what actually gets used.
    table['L'] = {'CTA': 0.0, 'CTC': 0.0, 'CTG': 1.0, 'CTT': 0.0, 'TTA': 0.0, 'TTG': 0.0}
    seq = "ATG" + "TTA" * 5 + "TAA"

    final_seq, _, num_edits = optimization_engine(seq, organism_name="kompas", codon_usage_table=table)

    # Translation preserved (EnforceTranslation) - not asserting the exact
    # stop codon, since 'use_best_codon' optimizes that position too, and
    # cub_c1's own real stop-codon frequencies (unmodified by this test)
    # may legitimately prefer a different stop codon than the original.
    from Bio.Seq import Seq
    assert str(Seq(final_seq).translate()) == str(Seq(seq).translate())
    leu_codons = [final_seq[i:i + 3] for i in range(3, len(final_seq) - 3, 3)]
    assert all(codon == "CTG" for codon in leu_codons)
    assert num_edits > 0


def test_codon_usage_table_is_ignored_when_custom_score_fn_is_given():
    # A deliberately broken table (missing Met/stop entries entirely) would
    # crash immediately if it were actually used to score the sequence -
    # completing successfully is itself the proof that custom_score_fn took
    # precedence and codon_usage_table was ignored, matching organism_name's
    # own documented precedence rule.
    broken_table = {'A': {'GCT': 1.0}}
    seq = "ATG" + "TTA" * 5 + "TAA"

    final_seq, _, _ = optimization_engine(
        seq, codon_usage_table=broken_table, custom_score_fn=lambda s: s.count("G"))

    assert len(final_seq) == len(seq)


def test_recombination_avoidance_breaks_the_near_duplicate_pair():
    site = "ACGTGGCATTAGCTAGCCTA"
    spacer = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
    seq = "ATG" + site + spacer + site + "TGA"
    # translate-preserving frame isn't guaranteed here since this isn't a
    # real ORF-aligned repeat; use exclusion_regions=() and no orf to just
    # check the constraint mechanics via direct df_recombination injection
    start_1, end_1 = 3, 3 + len(site)
    start_2 = 3 + len(site) + len(spacer)
    end_2 = start_2 + len(site)

    df_recombination = pd.DataFrame([{
        "sequence_1": site, "start_1": start_1, "end_1": end_1,
        "sequence_2": site, "start_2": start_2, "end_2": end_2,
        "location_delta": start_2 - end_1, "site_length": len(site),
        "log10_prob_recombination_ecoli": -1.0,
    }])

    final_seq, _, num_edits = optimization_engine(seq, df_recombination=df_recombination)

    assert num_edits > 0
    assert final_seq[start_1:end_1] != final_seq[start_2:end_2]


def test_slippage_avoidance_disrupts_the_repeat():
    seq = "ATG" + "GCT" * 12 + "TGA"  # 12x the same codon - a severe SSR hotspot
    df_slippage = pd.DataFrame([{
        "start": 3, "end": 3 + 36, "length_base_unit": 3, "sequence": "GCT" * 12,
        "num_base_units": 12, "log10_prob_slippage_ecoli": -1.0,
    }])

    final_seq, _, num_edits = optimization_engine(
        seq, df_slippage=df_slippage, orf_regions=[(0, len(seq))])

    assert num_edits > 0
    assert final_seq[3:39] != "GCT" * 12


def test_motif_avoidance_actually_removes_the_motif():
    # regression test for the end_index off-by-one bug (see module docstring) -
    # GC-balanced flanks so EnforceGCContent doesn't force unrelated edits
    # that could mask whether the motif constraint itself does anything.
    seq = "ATG" + "ATC" * 10 + "GATC" + "ATC" * 10 + "TAA"
    start = seq.index("GATC", 3)
    df_motifs = pd.DataFrame([{
        "start_index": start, "end_index": start + 4 - 1, "matching_motif": "dam",
        "PSSM_score": 5.0, "actual_site": "GATC", "actual_site_reverse_conjugate": "GATC",
    }])

    final_seq, _, num_edits = optimization_engine(seq, df_motifs=df_motifs)

    assert num_edits > 0
    assert final_seq[start:start + 4] != "GATC"


def test_motif_avoidance_respects_exclusion_regions():
    # the motif sits fully inside a locked region - it must survive untouched.
    seq = "ATG" + "ATC" * 10 + "GATC" + "ATC" * 10 + "TAA"
    start = seq.index("GATC", 3)
    df_motifs = pd.DataFrame([{
        "start_index": start, "end_index": start + 4 - 1, "matching_motif": "dam",
        "PSSM_score": 5.0, "actual_site": "GATC", "actual_site_reverse_conjugate": "GATC",
    }])

    final_seq, _, _ = optimization_engine(
        seq, df_motifs=df_motifs, exclusion_regions=[(start, start + 4)])

    assert final_seq[start:start + 4] == "GATC"


def test_non_codon_aligned_slippage_site_no_longer_crashes_dnachisel():
    # regression test for a real DNAChisel-internal crash: an AvoidPattern
    # window that doesn't align to the codon grid under EnforceTranslation
    # (e.g. a 2nt "GC" repeat starting mid-codon) could make
    # DnaOptimizationProblem.resolve_constraint compute a mutation-space
    # "choices span" that doesn't overlap the constraint's own location at
    # all - AvoidPattern.localized() correctly returns None for that, but
    # DNAChisel's caller doesn't check before calling .evaluate() on it,
    # crashing with "AttributeError: 'NoneType' object has no attribute
    # 'evaluate'". Confirmed this reproduces with codon-unaligned length-1
    # (homopolymer) and length-2 repeats alike, but NOT with a codon-aligned
    # length-3 repeat (see test_slippage_avoidance_disrupts_the_repeat, which
    # never crashed).
    #
    # Fixed not by widening the avoid-window (an earlier attempt at that
    # broadened "avoid this pattern at this spot" into "avoid this pattern
    # ANYWHERE in the whole codon", which for a single-nucleotide pattern in a
    # homopolymer run can be unsatisfiable for every codon in the run and
    # silently drop every constraint with no effect - see
    # test_non_codon_aligned_homopolymer_falls_back_to_dropping_the_site
    # below), but by catching this exact crash in optimize.py's retry loop
    # and treating it like any other unsatisfiable constraint: re-evaluate
    # everything and drop whichever constraints are still actually failing.
    seq = "ATG" + "GC" * 10 + "TAA"  # "GC" windows start at position 3, mid-codon
    df_slippage = pd.DataFrame([
        {"start": s, "end": s + 2, "length_base_unit": 2, "sequence": "GC",
         "num_base_units": 2, "log10_prob_slippage_ecoli": -1.0}
        for s in range(3, 23, 4)
    ])

    final_seq, _, num_edits = optimization_engine(
        seq, df_slippage=df_slippage, organism_name="kompas")

    assert len(final_seq) == len(seq)
    assert final_seq[:3] == "ATG" and final_seq[-3:] in ("TAA", "TAG", "TGA")
    assert num_edits > 0


def test_non_codon_aligned_homopolymer_falls_back_to_dropping_the_site():
    # same crash as above, but for a 1nt homopolymer run of "T" landing on
    # all-Phe codons (TTT/TTC, both of which always contain a T). There is no
    # way to satisfy "avoid a T here" while preserving translation over every
    # position of a pure poly-T run under this codon usage table, so - unlike
    # the "GC" case above, where other synonymous codons exist - every one of
    # these per-position constraints is expected to end up dropped rather
    # than resolved. The key regression this guards against is the crash
    # itself (and cnst.remove(None) on DNAChisel's own final-check error) -
    # not that every site gets fixed.
    seq = "ATG" + "T" * 12 + "TAA"
    df_slippage = pd.DataFrame([{
        "start": 3, "end": 15, "length_base_unit": 1, "sequence": "T" * 12,
        "num_base_units": 12, "log10_prob_slippage_ecoli": -1.0,
    }])

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        final_seq, _, _ = optimization_engine(
            seq, df_slippage=df_slippage, organism_name="kompas")

    assert len(final_seq) == len(seq)
    assert final_seq[:3] == "ATG" and final_seq[-3:] in ("TAA", "TAG", "TGA")

    # the user must be told this happened, and roughly why - not just left to
    # infer it from a lower-than-expected num_edits.
    messages = [str(w.message) for w in caught]
    assert any("Could not satisfy constraint" in m for m in messages)
    assert any("translation preservation" in m for m in messages)


def test_constraints_to_drop_picks_one_least_severe_constraint_when_a_protected_one_is_at_risk():
    # Core selection rule, protected-conflict branch: when a hard, protected
    # constraint (EnforceGCContent) is among what's currently failing, only
    # ONE non-protected candidate comes back - the least severe (ascending
    # risk score) currently-failing one - never the protected constraint
    # itself, and never a constraint that isn't currently failing at all (the
    # absent-pattern one below, despite being the "least risky" of everything
    # present, is never touched since it always passes).
    seq = "ATG" + "AAA" * 5 + "TAA"
    present_low_risk = dnachisel.AvoidPattern("AAA", location=(3, 6))
    present_low_risk.eso_severity = -5.0  # more negative log10(prob) = less risky
    present_high_risk = dnachisel.AvoidPattern("AAA", location=(6, 9))
    present_high_risk.eso_severity = -1.0
    absent_lowest_risk = dnachisel.AvoidPattern("GGGGGG", location=(9, 15))
    absent_lowest_risk.eso_severity = -9.0  # never present in seq -> always passes
    gc_constraint = dnachisel.EnforceGCContent(mini=1.0, maxi=1.0)  # impossible -> always fails
    cnst = [gc_constraint, present_low_risk, present_high_risk, absent_lowest_risk]
    problem = dnachisel.DnaOptimizationProblem(sequence=seq, constraints=cnst)

    assert _constraints_to_drop(cnst, problem) == [present_low_risk]


def test_constraints_to_drop_batches_when_no_protected_constraint_is_at_risk():
    # When nothing protected is in conflict - just several mutually-failing
    # hotspot constraints - every currently-failing one comes back at once
    # (not one-at-a-time): there's no hard requirement to protect here, and
    # batching keeps a dense conflict from burning one retry-budget round per
    # constraint (a real regression risk found while building this: a ~100-
    # constraint crash case blew through the scaled retry budget once dropping
    # became strictly one-at-a-time).
    seq = "ATG" + "AAA" * 5 + "TAA"
    present_1 = dnachisel.AvoidPattern("AAA", location=(3, 6))
    present_2 = dnachisel.AvoidPattern("AAA", location=(6, 9))
    translation_constraint = dnachisel.EnforceTranslation(location=(0, len(seq)))  # passes, not in conflict
    cnst = [present_1, present_2, translation_constraint]
    problem = dnachisel.DnaOptimizationProblem(sequence=seq, constraints=cnst)

    dropped = _constraints_to_drop(cnst, problem)
    assert len(dropped) == 2 and present_1 in dropped and present_2 in dropped


def test_constraints_to_drop_returns_empty_if_only_protected_constraints_fail():
    seq = "ATG" + "AAA" * 5 + "TAA"
    gc_constraint = dnachisel.EnforceGCContent(mini=1.0, maxi=1.0)  # impossible -> always fails
    translation_constraint = dnachisel.EnforceTranslation(location=(0, len(seq)))
    cnst = [gc_constraint, translation_constraint]
    problem = dnachisel.DnaOptimizationProblem(sequence=seq, constraints=cnst)

    assert _constraints_to_drop(cnst, problem) == []


def test_enforce_gc_content_conflict_is_resolved_by_dropping_the_lower_severity_site_first(monkeypatch):
    # The redesigned behavior (per explicit feedback: don't fail loudly the
    # moment EnforceGCContent is implicated - try resolving the clash by
    # dropping some other, less-critical constraint first, in order of
    # severity, and only remove what's actually necessary). Simulates a
    # GC-content conflict that resolves as soon as both real, currently-
    # present hotspot sites are gone, and asserts the LOWER-severity one
    # (-5.0, less risky) is tried and dropped before the HIGHER-severity one
    # (-1.0, riskier) - not raised immediately, and not dropped in an
    # arbitrary/construction order.
    seq = "ATG" + "AAA" * 5 + "TAA"
    dropped = []
    monkeypatch.setattr(
        "eso.optimize._warn_dropped_constraint",
        lambda constraint: dropped.append(constraint),
    )

    def fake_resolve_constraints(self, *args, **kwargs):
        if any(isinstance(c, dnachisel.AvoidPattern) for c in self.constraints):
            gc_constraint = next(c for c in self.constraints if isinstance(c, dnachisel.EnforceGCContent))
            raise NoSolutionError("simulated GC-content conflict", problem=self, constraint=gc_constraint)
        # both hotspot sites gone - the simulated conflict is resolved.

    monkeypatch.setattr(dnachisel.DnaOptimizationProblem, "resolve_constraints", fake_resolve_constraints)

    df_slippage = pd.DataFrame([
        {"start": 3, "end": 6, "length_base_unit": 3, "sequence": "AAA",
         "num_base_units": 2, "log10_prob_slippage_ecoli": -1.0},  # higher risk (kept longer)
        {"start": 6, "end": 9, "length_base_unit": 3, "sequence": "AAA",
         "num_base_units": 2, "log10_prob_slippage_ecoli": -5.0},  # lower risk (dropped first)
    ])

    # window_size_gc must actually fit inside this short test sequence - the
    # default (50) is longer than the whole sequence, so DNAChisel can't slide
    # any window and EnforceGCContent's own evaluate() vacuously passes,
    # which would hide the real conflict this test needs to simulate.
    optimization_engine(
        seq, mini_gc=0.3, maxi_gc=0.7, window_size_gc=10,
        df_slippage=df_slippage, organism_name="not_specified",
    )

    # In this scenario BOTH sites are genuinely necessary (the fake resolver
    # only succeeds once every AvoidPattern is gone), so the shrink/restore
    # pass below can't un-drop either one - both stay dropped, but restore is
    # attempted highest-severity-first, so the warning order is reversed
    # from drop order: -1.0 (tried to restore first, fails) then -5.0.
    assert len(dropped) == 2
    assert [round(c.eso_severity, 1) for c in dropped] == [-1.0, -5.0]


def test_gc_content_shrink_pass_restores_a_constraint_that_turned_out_unnecessary(monkeypatch):
    # The real gap the growth-only design (above) left open: the greedy walk
    # drops constraints in ascending-severity order as it goes, but never
    # re-checks whether an EARLIER drop is still needed once LATER drops
    # have also happened - so a constraint can end up permanently dropped
    # even though, given everything else that ended up dropped too, it was
    # never actually necessary on its own. Simulates exactly that: dropping
    # A (severity -9, tried first - least severe) doesn't help on its own,
    # so the walk also drops B (severity -1, the real blocker) - but once B
    # is gone, A was never actually needed. The shrink pass must restore A
    # and leave only B dropped (and warned about).
    seq = "ATG" + "AAA" * 5 + "TAA"

    a_constraint = dnachisel.AvoidPattern("GATC")
    a_constraint.eso_severity = -9.0
    b_constraint = dnachisel.AvoidPattern("GATC")
    b_constraint.eso_severity = -1.0

    class _FakeFailingEvaluation:
        passes = False

    # GC content and every AvoidPattern candidate always evaluate as
    # failing, regardless of real sequence content - _constraints_to_drop
    # must always see both A and B as live candidates each round.
    monkeypatch.setattr(dnachisel.EnforceGCContent, "evaluate", lambda self, problem: _FakeFailingEvaluation())
    monkeypatch.setattr(dnachisel.AvoidPattern, "evaluate", lambda self, problem: _FakeFailingEvaluation())

    def fake_resolve_constraints(self, *args, **kwargs):
        # "Resolved" iff B specifically is absent - A's presence or absence
        # never actually mattered, only the greedy walk doesn't know that.
        # DnaOptimizationProblem copies each constraint on construction (
        # confirmed directly: `problem.constraints[0] is original_constraint`
        # is False), so identity/equality checks against the original object
        # don't survive - checking the custom eso_severity attribute instead,
        # which DOES survive the copy (also confirmed directly).
        if any(getattr(c, 'eso_severity', None) == b_constraint.eso_severity for c in self.constraints):
            raise NoSolutionError("simulated ambiguous failure", problem=self, constraint=None)

    monkeypatch.setattr(dnachisel.DnaOptimizationProblem, "resolve_constraints", fake_resolve_constraints)

    dropped = []
    monkeypatch.setattr("eso.optimize._warn_dropped_constraint", lambda constraint: dropped.append(constraint))

    def fake_convert_df_to_constraints(df):
        return [a_constraint, b_constraint]

    monkeypatch.setattr("eso.optimize.convert_df_to_constraints", fake_convert_df_to_constraints)

    df_slippage = pd.DataFrame([{
        "start": 3, "end": 6, "length_base_unit": 3, "sequence": "AAA",
        "num_base_units": 2, "log10_prob_slippage_ecoli": -1.0,
    }])

    final_seq, _, _ = optimization_engine(
        seq, mini_gc=0.3, maxi_gc=0.7, df_slippage=df_slippage, organism_name="not_specified")

    assert final_seq == seq  # fake resolver never mutates anything
    # A was dropped during growth (tried first, least severe) but restored
    # once the shrink pass confirmed B alone was the actual blocker.
    assert dropped == [b_constraint]


def test_shrink_pass_is_skipped_above_the_verification_cap(monkeypatch):
    # Restoring each dropped constraint costs one full resolve_constraints()
    # call, and that compounds with the retry loop's own already-quadratic
    # worst case for a genuinely dense conflict (benchmarked directly: a
    # near-0%-GC poly-T stress sequence with ~475 dropped constraints took
    # ~180s with unconditional shrink verification). Above
    # _SHRINK_VERIFICATION_CAP, the shrink pass must be skipped entirely -
    # every dropped constraint stays dropped and warned about, with zero
    # _try_restore calls, exactly like before the shrink pass existed.
    from eso.optimize import _SHRINK_VERIFICATION_CAP

    seq = "ATG" + "AAA" * 5 + "TAA"
    constraints = []
    for i in range(_SHRINK_VERIFICATION_CAP + 5):
        c = dnachisel.AvoidPattern("GATC")
        c.eso_severity = float(i)
        constraints.append(c)

    class _FakeFailingEvaluation:
        passes = False

    monkeypatch.setattr(dnachisel.AvoidPattern, "evaluate", lambda self, problem: _FakeFailingEvaluation())

    def fake_resolve_constraints(self, *args, **kwargs):
        # Resolved only once every AvoidPattern is gone - none of them
        # matters individually here, only the total count relative to the cap.
        if any(isinstance(c, dnachisel.AvoidPattern) for c in self.constraints):
            raise NoSolutionError("simulated ambiguous failure", problem=self, constraint=None)

    monkeypatch.setattr(dnachisel.DnaOptimizationProblem, "resolve_constraints", fake_resolve_constraints)
    monkeypatch.setattr("eso.optimize.convert_df_to_constraints", lambda df: constraints)

    restore_calls = []
    monkeypatch.setattr(
        "eso.optimize._try_restore",
        lambda cnst, constraint, obj, seq: restore_calls.append(constraint) or None,
    )
    dropped = []
    monkeypatch.setattr("eso.optimize._warn_dropped_constraint", lambda c: dropped.append(c))

    df_slippage = pd.DataFrame([{
        "start": 3, "end": 6, "length_base_unit": 3, "sequence": "AAA",
        "num_base_units": 2, "log10_prob_slippage_ecoli": -1.0,
    }])

    optimization_engine(seq, mini_gc=0.3, maxi_gc=0.7, df_slippage=df_slippage, organism_name="not_specified")

    assert len(dropped) == _SHRINK_VERIFICATION_CAP + 5
    assert restore_calls == []


def test_enforce_gc_content_is_never_silently_dropped_when_no_alternative_helps(monkeypatch):
    # The narrower, still-real gap this whole fix is ultimately about: if
    # dropping every other, less-critical constraint still can't rescue
    # EnforceGCContent (or nothing else is even in conflict), optimize.py must
    # still refuse to silently drop it and ship a sequence outside the
    # requested GC range - it raises instead. Forces EnforceGCContent's own
    # evaluate() to always report failure, with no hotspot/enzyme constraint
    # present at all to try dropping instead.
    seq = "ATG" + "AAA" * 5 + "TAA"

    def fake_resolve_constraints(self, *args, **kwargs):
        raise NoSolutionError("simulated ambiguous failure", problem=self, constraint=None)

    monkeypatch.setattr(dnachisel.DnaOptimizationProblem, "resolve_constraints", fake_resolve_constraints)

    class _FakeFailingEvaluation:
        passes = False

    monkeypatch.setattr(dnachisel.EnforceGCContent, "evaluate", lambda self, problem: _FakeFailingEvaluation())

    with pytest.raises(NoSolutionError, match="hard, user-specified requirement"):
        optimization_engine(seq, mini_gc=0.3, maxi_gc=0.7, organism_name="not_specified")


def test_retry_budget_scales_past_the_old_fixed_60_round_cap():
    # Smoke test, not a proven regression reproduction: I could not construct
    # a case that actually failed under the old fixed `flag < 60` (checked
    # directly - DNAChisel's own "drop every currently-failing constraint in
    # one round" path, eso.optimize's `drop_and_retry` branch, handles this
    # exact homopolymer scenario efficiently regardless of how many point
    # constraints exist, so this scenario alone doesn't need more than a
    # handful of rounds either way). The scaling fix (`retry_budget =
    # max(60, len(cnst))`) is still a real, deliberate defensive improvement
    # for the OTHER retry path this loop has - DNAChisel surfacing one
    # specific named unsatisfiable constraint per top-level
    # resolve_constraints() call (the `e.constraint is not None` branch),
    # which drops exactly one constraint per round - since that path scales
    # with constraint count, not with rounds-per-batch, if it's ever the one
    # DNAChisel takes for a dense-enough sequence. I did not find a
    # deterministic way to force DNAChisel down that specific path with more
    # than 60 named failures from the outside. This test just confirms a
    # large (~100 point constraints), all-unsatisfiable case still resolves
    # correctly with the new budget in place - a coverage-at-scale check, not
    # proof the old cap was reachable.
    # num_t must be a multiple of 3, so the trailing "TAA" actually lands
    # in-frame as a real, protected stop codon. A previous version of this
    # test used num_t=200 (not a multiple of 3): the resulting 206nt
    # sequence isn't a multiple of 3 either, so optimization_engine's own
    # default orf_regions computation (`(len(seq) // 3) * 3`) correctly
    # truncated the "ORF" to 204nt, leaving the last 2 characters of the
    # string - part of what looked like "TAA" - genuinely OUTSIDE
    # EnforceTranslation's protection, free for GC-content optimization to
    # mutate. That's a real bug in the TEST FIXTURE (an accidentally
    # frame-shifted, not-really-protected "stop codon"), not in the retry
    # loop: confirmed directly that with a properly frame-aligned sequence,
    # the real default GC bounds (0.3-0.7, restored below - a previous fix
    # here had incorrectly widened them to 0.0-1.0 to dodge this same
    # flakiness instead of finding its actual cause) are stable across many
    # runs.
    num_t = 201  # ~100 point constraints via modify_df_slippage's "every other unit"
    seq = "ATG" + "T" * num_t + "TAA"
    assert len(seq) % 3 == 0
    df_slippage = pd.DataFrame([{
        "start": 3, "end": 3 + num_t, "length_base_unit": 1, "sequence": "T" * num_t,
        "num_base_units": num_t, "log10_prob_slippage_ecoli": -1.0,
    }])

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        final_seq, _, _ = optimization_engine(seq, df_slippage=df_slippage, organism_name="kompas")

    assert len(final_seq) == len(seq)
    assert final_seq[:3] == "ATG" and final_seq[-3:] in ("TAA", "TAG", "TGA")


def test_custom_score_is_scoped_to_orf_regions_not_the_whole_sequence():
    # regression test for a real bug: custom scoring applied across the WHOLE
    # sequence regardless of orf_regions, unlike the built-in CAI/tAI codon
    # scoring (which is correctly scoped via CodonOptimize(location=orf, ...)).
    # Confirmed directly before fixing: with orf_regions=[(0, 36)] on a 41nt
    # sequence, CustomScore's own location was 0-41, not 0-36 - a custom score
    # would get called on non-ORF flanking nucleotides too.
    import dnachisel
    from eso.custom_score import CustomScore

    seq = "ATG" + "GCT" * 10 + "TAA" + "AAAAA"  # 41nt total; ORF is 0-36
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obj = [CustomScore(lambda s: 0, location=orf) for orf in [(0, 36)]]
    problem = dnachisel.DnaOptimizationProblem(sequence=seq, objectives=obj)
    initialized = obj[0].initialized_on_problem(problem, role="objective")

    assert (initialized.location.start, initialized.location.end) == (0, 36)

    seen = []
    initialized.score_fn = lambda c: seen.append(c) or 0
    initialized.evaluate(problem)

    # score_fn saw exactly the ORF, nothing from the flanking "AAAAA" tail.
    assert seen == ["ATG" + "GCT" * 10 + "TAA"]


# --- avoid_hairpins / avoid_enzymes (previously unused DNAChisel features) -

def test_avoid_hairpins_removes_a_real_hairpin():
    # DNAChisel's own AvoidHairpins/EnzymeSitePattern were installed,
    # unused dependencies before this - a real gap for standard gene-design
    # practice (avoiding secondary structure, and specific restriction
    # sites for cloning). A stem immediately followed by its own reverse
    # complement is a textbook hairpin - confirmed directly (before adding
    # avoid_hairpins support) that DNAChisel's own AvoidHairpins reports
    # this exact construction as a real, detected hairpin.
    stem = "ACGTACGTACGTACG"  # 15nt
    spacer = "T" * 9
    hairpin_region = stem + spacer + reverse_complement(stem)
    seq = "ATG" + hairpin_region + "TAA"
    assert len(seq) % 3 == 0

    # only the start codon is translation-protected, so AvoidHairpins is the
    # only real constraint on the hairpin region itself - isolates exactly
    # what's being tested from any translation-preservation interaction.
    final_seq, _, num_edits = optimization_engine(
        seq, organism_name='not_specified', orf_regions=((0, 3),),
        avoid_hairpins=True, hairpin_stem_size=15, hairpin_window=60,
        mini_gc=0.0, maxi_gc=1.0)

    assert final_seq[:3] == "ATG"
    assert num_edits > 0
    problem_after = dnachisel.DnaOptimizationProblem(sequence=final_seq, constraints=[])
    hairpin_check = dnachisel.AvoidHairpins(
        stem_size=15, hairpin_window=60).initialized_on_problem(problem_after, role='constraint')
    assert hairpin_check.evaluate(problem_after).passes


def test_avoid_enzymes_removes_the_recognition_site():
    seq = "ATG" + "GAATTC" * 3 + "TAA"  # a real EcoRI site (GAATTC)

    final_seq, _, num_edits = optimization_engine(
        seq, organism_name='not_specified', orf_regions=((0, 3),),
        avoid_enzymes=['EcoRI'], mini_gc=0.0, maxi_gc=1.0)

    assert final_seq[:3] == "ATG"
    assert num_edits > 0
    assert "GAATTC" not in final_seq


def test_avoid_enzymes_rejects_an_unrecognized_name():
    seq = "ATG" + "AAA" * 5 + "TAA"
    with pytest.raises(ValueError, match="NotARealEnzyme"):
        optimization_engine(seq, avoid_enzymes=['NotARealEnzyme'])
