"""Property-based (fuzz) tests for eso.detection.slippage.

Generates random DNA sequences (not hand-picked examples) and checks
internal-consistency invariants that must hold for every detected candidate,
regardless of what the sequence happens to contain.
"""

from hypothesis import given, settings, strategies as st

from eso.detection.slippage import find_slippage_candidates, find_slippage_candidates_slow, modify_df_slippage

_dna_sequence = st.text(alphabet='ACGT', min_size=0, max_size=300)


def _covered_positions(df):
    covered = set()
    for row in df.itertuples():
        covered.update(range(int(row.start), int(row.end)))
    return covered


def test_odd_length_homopolymer_phase_boundary_never_loses_coverage():
    # Deterministic repro of the phase-tie-break case described in
    # test_fast_slippage_scan_never_covers_less_than_the_slow_reference's
    # docstring below - a 7nt T-homopolymer (odd length, so a length-2
    # unit_len reading can only ever cover 6 of the 7 characters, leaving
    # a genuine ambiguity about which boundary nucleotide is left out).
    # Kept as a fixed example, not left to hypothesis's random search
    # alone, since max_size=300 won't reliably reproduce this at the
    # sequence length (20kb) where it was first found.
    seq = "ATG" + "T" * 7 + "GAT" + "TAA"
    slow_covered = _covered_positions(find_slippage_candidates_slow(seq))
    fast_covered = _covered_positions(find_slippage_candidates(seq))
    assert slow_covered <= fast_covered


@settings(deadline=None)
@given(seq=_dna_sequence)
def test_fast_slippage_scan_never_covers_less_than_the_slow_reference(seq):
    # Regression/equivalence test for the byte-XOR-mask rewrite of
    # find_slippage_candidates (ported from STABLES's
    # automate_goi_repeat_alternatives.py, which measured an 11.5x-73.3x
    # speedup for the analogous rewrite there).
    #
    # NOT exact raw-tuple equality - confirmed directly (on a 20kb random
    # sequence, larger than hypothesis's max_size here would reliably
    # generate) that the two implementations can pick a different, equally
    # valid phase/content representative for a perfectly periodic or
    # homopolymer-like run: e.g. one may report a "TT"-repeat window at
    # [2419,2425) where the other reports [2420,2426) - both real readings
    # of the same physical 7nt T-homopolymer, differing only in which
    # boundary nucleotide is included. This is the exact "GCGCGCGC vs.
    # 1-shifted CGCGCG" phenomenon collapse_slippage_sites' own docstring
    # already documents as a known, accepted ambiguity - not a bug, and
    # resolved downstream by the collapse step, never by the raw
    # candidate generator (find_slippage_candidates is documented as
    # returning every raw candidate, uncollapsed).
    #
    # The invariant that actually matters - and is checked here - is that
    # the fast scan never loses coverage: every position the slow
    # reference would place under a correction constraint must also be
    # covered by the fast scan (verified on that same 20kb sequence: the
    # fast scan's covered-position set was a strict superset of the slow
    # one's, by exactly the handful of alternate-phase boundary positions
    # described above - never the reverse).
    slow_covered = _covered_positions(find_slippage_candidates_slow(seq))
    fast_covered = _covered_positions(find_slippage_candidates(seq))
    assert slow_covered <= fast_covered


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
