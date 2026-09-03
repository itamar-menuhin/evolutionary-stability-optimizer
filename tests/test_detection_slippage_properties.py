"""Property-based (fuzz) tests for eso.detection.slippage.

Generates random DNA sequences (not hand-picked examples) and checks
internal-consistency invariants that must hold for every detected candidate,
regardless of what the sequence happens to contain.
"""

from hypothesis import given, settings, strategies as st

from eso.detection.slippage import find_slippage_candidates, modify_df_slippage

_dna_sequence = st.text(alphabet='ACGT', min_size=0, max_size=300)


@settings(deadline=None)  # detection itself is the thing under test, not its speed
@given(seq=_dna_sequence)
def test_every_candidate_start_end_and_sequence_are_mutually_consistent(seq):
    candidates = find_slippage_candidates(seq)

    for row in candidates.itertuples():
        # the reported `sequence` must be the real substring at (start, end) -
        # not a stale or miscomputed copy.
        assert row.sequence == seq[row.start:row.end]
        # a candidate's span is always exactly num_base_units whole repeat
        # units of length_base_unit, never a partial unit.
        assert len(row.sequence) == row.num_base_units * row.length_base_unit
        assert 0 <= row.start < row.end <= len(seq)


@settings(deadline=None)
@given(seq=_dna_sequence)
def test_every_candidate_passes_the_declared_score_filter(seq):
    # find_slippage_candidates documents filtering to log10_prob > -9 - a
    # fuzz test is exactly the way to catch a candidate that slipped through
    # unfiltered (e.g. a boundary/rounding mistake in the formula, or a
    # length-1 vs. length>1 branch mix-up).
    candidates = find_slippage_candidates(seq)
    assert (candidates.log10_prob_slippage_ecoli > -9).all()


@settings(deadline=None)
@given(seq=_dna_sequence)
def test_modify_df_slippage_only_ever_targets_positions_inside_the_original_site(seq):
    candidates = find_slippage_candidates(seq)
    modified = modify_df_slippage(candidates)

    # every scattered avoidance site must itself be a real substring at its
    # own (start, end) - modify_df_slippage slices `sequence` by index
    # arithmetic (i*length : (i+1)*length), exactly the shape of bug this
    # codebase has hit before (off-by-one range conventions).
    for row in modified.itertuples():
        assert row.sequence == seq[row.start:row.end]
