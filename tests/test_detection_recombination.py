import pytest

from eso.detection.recombination import find_recombination_sites, calc_recombination_score

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


def test_calc_recombination_score_decreases_with_distance():
    close = calc_recombination_score(location_delta=5, site_length=20)
    far = calc_recombination_score(location_delta=500, site_length=20)
    assert close > far


def test_calc_recombination_score_matches_efm_calculator_reference():
    # Pins the formula to the authoritative reference implementation
    # (github.com/barricklab/efm-calculator, get_recombo_rate, E. coli/yeast
    # case) - not just our own internal consistency. Confirmed via that
    # repo's git history that a=8.8 has been the value since its first
    # commit (2015-03-31); this codebase had inherited a wrong a=5.8 from
    # ESO_curr/STABLES (see calc_recombination_score's docstring).
    import math

    def efm_reference(location_delta, site_length):
        return math.log10(
            ((8.8 + location_delta) ** (-29.0 / site_length))
            * (site_length / (1 + 1465.6 * site_length))
        )

    for location_delta, site_length in [(1, 16), (5, 20), (100, 30), (500, 50)]:
        ours = calc_recombination_score(location_delta, site_length)
        reference = efm_reference(location_delta, site_length)
        assert ours == pytest.approx(reference, abs=1e-9)
