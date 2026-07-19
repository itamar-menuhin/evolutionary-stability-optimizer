from eso.detection.recombination import find_recombination_sites, calc_recombination_score


def test_finds_exact_duplicate_site():
    site = "ACGTGGCATTAGCTAG"  # 16nt
    seq = "ATGCATGCAT" + site + "CCCCCCCCCCCCCCCCCCCC" + site + "TTGGCCAATT"

    df = find_recombination_sites(seq)

    assert not df.empty
    # results use a sequence_1/sequence_2 pair schema, not a single 'sequence' column
    assert (df.sequence_1.str.contains(site) | df.sequence_2.str.contains(site)).any()


def test_no_recombination_sites_in_non_repetitive_sequence():
    # Hand-checked to contain no repeated/near-duplicate 16+nt window.
    seq = "ACGTAGCTTGACCTGAAGCTAGCATTGCA"
    df = find_recombination_sites(seq)
    assert df.empty


def test_calc_recombination_score_decreases_with_distance():
    close = calc_recombination_score(location_delta=5, site_length=20)
    far = calc_recombination_score(location_delta=500, site_length=20)
    assert close > far
