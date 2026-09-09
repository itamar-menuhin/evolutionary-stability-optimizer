import pandas as pd
import pytest

from eso.detection.recombination import (
    find_recombination_sites,
    find_recombination_candidates,
    calc_recombination_score,
    collapse_recombination_sites,
    recombination_sites_for_constraints,
    _elongate_sites,
)
from eso.sequence_utils import reverse_complement_seq

# non-repetitive spacer: a homopolymer/simple-repeat spacer would itself be a
# genuine (distinct) slippage/recombination hotspot and pollute row counts.
SPACER = "ATCGGATCCAAGCTTGGATCCAAGCTTGGA"


def test_finds_exact_duplicate_site():
    site = "ACGTGGCATTAGCTAGCCTA"  # 20nt
    seq = "ATGCATGCAT" + site + SPACER + site + "TTGGCCAATT"

    df = find_recombination_sites(seq)

    # collapsed to exactly one row: many candidate seeds converge, via
    # elongation, on this same real hotspot and must not be reported separately
    assert df.shape[0] == 1
    row = df.iloc[0]
    assert site in row.sequence_1 or site in row.sequence_2


def test_finds_near_duplicate_with_centered_mutation():
    # A single substitution dead-center (prefix=12nt, suffix=11nt) leaves no
    # 16-consecutive-nt exact match on either side - this is the case the
    # exact-match (STABLES) variant cannot catch, but the Levenshtein-tolerant
    # primary detector should, since the two sites are still edit distance 1 apart.
    site = "ACGTGGCATTAGGCTAGCCTAGGC"
    site_mut = "ACGTGGCATTAGACTAGCCTAGGC"
    seq = "ATGCATGCAT" + site + SPACER + site_mut + "TTGGCCAATT"

    df = find_recombination_sites(seq)

    assert df.shape[0] == 1
    row = df.iloc[0]
    assert {row.sequence_1, row.sequence_2} == {site, site_mut}


def test_no_recombination_sites_in_non_repetitive_sequence():
    # Hand-checked to contain no repeated/near-duplicate 16+nt window.
    seq = "ACGTAGCTTGACCTGAAGCTAGCATTGCA"
    df = find_recombination_sites(seq)
    assert df.empty


def test_distinct_hotspots_are_not_merged_into_each_other():
    # Two unrelated exact-duplicate pairs, far apart, with a non-repetitive
    # sequence between and around them - the overlap-based collapse must not
    # accidentally merge genuinely distinct hotspots. Uses two DIFFERENT
    # spacers: reusing the same spacer text near both pairs would itself
    # create a third, genuine duplicate at the site/spacer boundary.
    site_a = "ACGTGGCATTAGCTAGCCTA"
    site_b = "TTGACCGGAATCCGTTAGCA"
    spacer_2 = "GTAGCTAACGATTGCGATCCGTAACTAGGA"
    seq = (
        "ATGCATGCAT" + site_a + SPACER + site_a
        + "GCGCGCTTAACC" + site_b + spacer_2 + site_b
        + "TTGGCCAATT"
    )

    df = find_recombination_sites(seq)

    assert df.shape[0] == 2
    found_sequences = set(df.sequence_1) | set(df.sequence_2)
    assert any(site_a in s for s in found_sequences)
    assert any(site_b in s for s in found_sequences)


def test_direct_repeat_candidates_are_labeled_direct_strand():
    # find_recombination_candidates' output must carry the new strand field
    # (issue #17) so a caller can audit that only direct repeats were used
    # to build correction constraints.
    site = "ACGTGGCATTAGCTAGCCTA"  # 20nt
    seq = "ATGCATGCAT" + site + SPACER + site + "TTGGCCAATT"

    df = find_recombination_candidates(seq)

    assert not df.empty
    assert set(df.strand) == {'direct'}


def test_inverted_repeat_pair_is_excluded_from_candidates():
    # issue #17: calc_recombination_score is Oliveira et al. 2008's formula
    # for RecA-mediated deletion between DIRECT repeats specifically. An
    # inverted repeat (a site paired with its own reverse complement) is a
    # mechanistically different hazard (hairpin/cruciform formation, not
    # RecA-mediated deletion) that this formula does not model, so it must
    # not surface anywhere in the real detection/scoring/correction output.
    site = "ACGTGGCATTAGCTAGCCTA"  # 20nt
    seq = "ATGCATGCAT" + site + SPACER + reverse_complement_seq(site) + "TTGGCCAATT"

    df = find_recombination_candidates(seq)

    assert df.empty
    # sanity check: this sequence really does contain a detectable
    # inverted-repeat pair pre-filtering, so an empty result here is the
    # strand filter doing its job, not an unrelated absence of any match.
    from eso.detection.recombination import _generate_relevant_pairs_fast
    raw_pairs = _generate_relevant_pairs_fast(seq)
    assert (raw_pairs.strand == 'inverted').any()


def test_direct_repeat_still_detected_alongside_an_unrelated_inverted_repeat():
    # Regression check that the strand filter added for issue #17 does not
    # affect same-strand (direct-repeat) detection: a genuine direct-repeat
    # hotspot is still found and scored even when the sequence also
    # contains an unrelated inverted-repeat pair elsewhere, and the
    # inverted pair must not leak into the output alongside it.
    site = "ACGTGGCATTAGCTAGCCTA"  # 20nt, genuine direct repeat
    site_b = "TTGACCGGAATCCGTTAGCA"  # distinct 20nt site, used as an inverted repeat
    spacer_2 = "GTAGCTAACGATTGCGATCCGTAACTAGGA"
    seq = (
        "ATGCATGCAT" + site + SPACER + site
        + "GCGCGCTTAACC" + site_b + spacer_2 + reverse_complement_seq(site_b)
        + "TTGGCCAATT"
    )

    df = find_recombination_candidates(seq)

    assert not df.empty
    assert set(df.strand) == {'direct'}
    found_sequences = set(df.sequence_1) | set(df.sequence_2)
    assert any(site in s for s in found_sequences)
    assert not any(site_b in s for s in found_sequences)
    assert not any(reverse_complement_seq(site_b) in s for s in found_sequences)


def test_partially_overlapping_pairs_both_kept_for_constraints():
    # Regression test for the same coverage-gap bug as
    # eso.detection.slippage's overlapping-but-not-nested case, applied to
    # recombination pairs: two candidate pairs whose site_1/site_2 ranges
    # each overlap but neither fully contains the other. Plain non-max
    # suppression (collapse_recombination_sites, the report view) correctly
    # keeps only the higher-scoring one - but recombination_sites_for_constraints
    # must keep both, since the lower-scoring pair's site_2 sticks out beyond
    # what the higher-scoring one covers.
    df = pd.DataFrame([
        {
            'sequence_1': 'A' * 16, 'start_1': 0, 'end_1': 16,
            'sequence_2': 'A' * 16, 'start_2': 100, 'end_2': 116,
            'log10_prob_recombination_ecoli': -1.0,
        },
        {
            'sequence_1': 'A' * 16, 'start_1': 0, 'end_1': 16,
            'sequence_2': 'A' * 20, 'start_2': 100, 'end_2': 120,  # sticks out past 116
            'log10_prob_recombination_ecoli': -2.0,
        },
    ])

    reported = collapse_recombination_sites(df)
    assert reported.shape[0] == 1
    assert reported.iloc[0].end_2 == 116  # higher-scoring pair only

    for_constraints = recombination_sites_for_constraints(df)
    assert set(for_constraints.end_2) == {116, 120}  # both pairs, [116, 120) not lost


def test_calc_recombination_score_decreases_with_distance():
    close = calc_recombination_score(location_delta=5, site_length=20)
    far = calc_recombination_score(location_delta=500, site_length=20)
    assert close > far


def test_calc_recombination_score_matches_efm_calculator_reference():
    # Pins the formula to the primary source (Oliveira et al. 2008, Plasmid
    # 60:159-165, Table 3, recA+ row: A=5.8, B=1465.6, a=29.0) - not just our
    # own internal consistency. See calc_recombination_score's docstring for
    # why the EFM Calculator tool's own hardcoded A=8.8 is a mix-up with a
    # different row/parameter of that same table, not a correction.
    import math

    def efm_reference(location_delta, site_length):
        return math.log10(
            ((5.8 + location_delta) ** (-29.0 / site_length))
            * (site_length / (1 + 1465.6 * site_length))
        )

    for location_delta, site_length in [(1, 16), (5, 20), (100, 30), (500, 50)]:
        ours = calc_recombination_score(location_delta, site_length)
        reference = efm_reference(location_delta, site_length)
        assert ours == pytest.approx(reference, abs=1e-9)


def test_elongate_sites_does_not_wrap_past_either_sequence_boundary():
    # Regression test for a real gap: with no explicit bounds guard, a site
    # pair starting at index 0 (nothing left to extend into) or ending at
    # the sequence's last index (nothing right to extend into) could have
    # start_1 or end_2 walked past the sequence's real boundaries. Python
    # doesn't raise for a negative/out-of-range string index in a slice - it
    # silently wraps to slice from the OTHER end of the string (or returns
    # ''), which would corrupt the reported site boundary rather than crash.
    # This site pair sits at both boundaries at once: site 1 starts at index
    # 0, site 2 ends at the sequence's last index.
    full_seq = "GATCGATCGATCGATC" + "AAAA" + "GATCGATCGATCGATC"
    row = pd.Series({
        "sequence_1": full_seq[0:16], "start_1": 0, "end_1": 15,
        "sequence_2": full_seq[20:36], "start_2": 20, "end_2": 35,
    })

    result = _elongate_sites(row, full_seq)

    assert 0 <= result.start_1 and result.end_1 <= len(full_seq)
    assert 0 <= result.start_2 and result.end_2 <= len(full_seq)
    assert full_seq[result.start_1:result.end_1] == result.sequence_1
    assert full_seq[result.start_2:result.end_2] == result.sequence_2
