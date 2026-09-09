"""Property-based (fuzz) tests for eso.detection.recombination's fast
candidate-pair generator, checked against the slow reference implementation.

Generates random DNA sequences (not hand-picked examples) and checks that
the fast (deletion-signature-based) path and the slow (exact edit-distance-1
variant) path converge on the exact same final result, after elongation and
scoring - the two are allowed to disagree on raw, pre-elongation seeds (both
can find a different sub-window of the same real site first), but not on
what gets reported.
"""

from hypothesis import given, settings, strategies as st

from eso.detection.recombination import (
    _generate_all_recombination_sites_slow,
    _generate_relevant_pairs_fast,
    _generate_relevant_pairs_slow,
    _elongate_sites,
    calc_recombination_score,
)
from eso.sequence_utils import reverse_complement_seq

_dna_sequence = st.text(alphabet='ACGT', min_size=0, max_size=200)


def _finalize(df_pairs, seq):
    # Mirrors find_recombination_candidates' own post-processing (elongate,
    # score, filter) - see that function in eso/detection/recombination.py.
    if df_pairs.shape[0] == 0:
        return df_pairs
    df_pairs = df_pairs.apply(lambda row: _elongate_sites(row, seq), axis=1)
    df_pairs['location_delta'] = df_pairs.start_2 - df_pairs.end_1
    df_pairs['site_length'] = df_pairs.apply(
        lambda row: max(row.end_1 - row.start_1, row.end_2 - row.start_2), axis=1)
    df_pairs['log10_prob_recombination_ecoli'] = df_pairs.apply(
        lambda x: calc_recombination_score(x.location_delta, x.site_length), axis=1)
    return df_pairs[df_pairs.log10_prob_recombination_ecoli > -9]


def _covered_positions(df):
    covered = set()
    for row in df.itertuples():
        covered.update(range(int(row.start_1), int(row.end_1)))
        covered.update(range(int(row.start_2), int(row.end_2)))
    return covered


def _slow_covered(seq):
    sites = _generate_all_recombination_sites_slow(seq)
    return _covered_positions(_finalize(_generate_relevant_pairs_slow(*sites), seq))


def _fast_covered(seq):
    return _covered_positions(_finalize(_generate_relevant_pairs_fast(seq), seq))


def test_convention_mismatch_no_longer_over_extends_past_the_sequence_end():
    # Regression test for a real bug, caught before shipping: the fast
    # path's window-generating helpers (_all_windows_both_strands, the
    # seed_windows in _generate_relevant_pairs_fast) originally labelled a
    # window's `end` using the exclusive/standard-slicing convention
    # (start + length), but _elongate_sites and the overlap check
    # (`end_a >= start_b`) both assume the slow path's convention of an
    # INCLUSIVE last-character index (start + length - 1). Feeding
    # exclusive-convention values into code written for the inclusive
    # convention silently over-extended every elongated site by one extra
    # position (confirmed: for this exact 40nt homopolymer, the fast path
    # reported a site extending to index 40 - one past the last valid
    # index, 39, of a 40nt sequence). Fixed by making the fast path's
    # window helpers use the same inclusive-end convention as the slow
    # path throughout.
    seq = "A" * 40
    fast_covered = _fast_covered(seq)
    assert fast_covered == _slow_covered(seq)
    assert max(fast_covered) < len(seq)


@settings(deadline=None)
@given(seq=_dna_sequence)
def test_fast_and_slow_agree_on_final_covered_positions(seq):
    assert _fast_covered(seq) == _slow_covered(seq)


def test_fast_and_slow_agree_on_targeted_edge_cases():
    site = "ACGTGGCATTAGCTAGCCTA"  # 20nt
    tests = {
        "exact duplicate": "ATGCATGCAT" + site + "ATCGGATCCAAGCTTGGATCCAAGCTTGGA" + site + "TTGGCCAATT",
        "homopolymer 40": "A" * 40,
        "tandem repeat AT*30": "AT" * 30,
        "single substitution": "ATGCATGCAT" + site + "ATCGGATCCAAGCTTGGATCCAAGCTTGGA"
        + site[:10] + "T" + site[11:] + "TTGGCCAATT",
        "single insertion": "ATGCATGCAT" + site + "ATCGGATCCAAGCTTGGATCCAAGCTTGGA"
        + site[:10] + "C" + site[10:] + "TTGGCCAATT",
        "single deletion": "ATGCATGCAT" + site + "ATCGGATCCAAGCTTGGATCCAAGCTTGGA"
        + site[:10] + site[11:] + "TTGGCCAATT",
        "inverted repeat (RC)": "ATGCATGCAT" + site + "ATCGGATCCAAGCTTGGATCCAAGCTTGGA"
        + reverse_complement_seq(site) + "TTGGCCAATT",
    }
    for name, seq in tests.items():
        assert _fast_covered(seq) == _slow_covered(seq), name
