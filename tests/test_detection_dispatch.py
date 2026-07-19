import pytest

from eso.detection.dispatch import find_recombination_sites, find_slippage_sites

SPACER = "ATCGGATCCAAGCTTGGATCCAAGCTTGGA"


def test_thorough_mode_is_default_and_matches_primary_module():
    from eso.detection.recombination import find_recombination_sites as primary_impl

    site = "ACGTGGCATTAGCTAGCCTA"
    seq = "ATGCATGCAT" + site + SPACER + site + "TTGGCCAATT"

    df_default = find_recombination_sites(seq)
    df_explicit = find_recombination_sites(seq, mode="thorough")
    df_primary = primary_impl(seq)

    assert df_default.shape[0] == df_explicit.shape[0] == df_primary.shape[0] == 1


def test_fast_mode_routes_to_staubility_variant():
    from eso.detection.staubility_variant import find_recombination_sites as fast_impl

    site = "ACGTGGCATTAGCTAGCCTA"
    seq = "ATGCATGCAT" + site + SPACER + site + "TTGGCCAATT"

    df_dispatch = find_recombination_sites(seq, mode="fast")
    df_fast = fast_impl(seq)

    assert df_dispatch.shape[0] == df_fast.shape[0] == 1


def test_thorough_catches_centered_near_duplicate_that_fast_misses():
    # The central case where the two modes genuinely disagree: no 16nt exact
    # window survives a dead-center substitution, so "fast" (exact-match)
    # cannot find it, but "thorough" (Levenshtein-tolerant) still can.
    site = "ACGTGGCATTAGGCTAGCCTAGGC"
    site_mut = "ACGTGGCATTAGACTAGCCTAGGC"
    seq = "ATGCATGCAT" + site + SPACER + site_mut + "TTGGCCAATT"

    df_thorough = find_recombination_sites(seq, mode="thorough")
    df_fast = find_recombination_sites(seq, mode="fast")

    assert df_thorough.shape[0] == 1
    assert df_fast.empty


def test_unknown_mode_raises():
    with pytest.raises(ValueError, match="Unknown recombination mode"):
        find_recombination_sites("ACGTACGTACGT", mode="nonexistent")


def test_slippage_default_and_fast_modes_agree():
    # Unlike recombination, both slippage modes are equivalent in sensitivity
    # (verified via a 300-trial fuzz sweep, see docs/detector-comparisons.md) -
    # this is purely a speed choice, so both modes should always agree.
    seq = "ATGCTAAT" + "GCGCGCGC" + "TTAGGCATGCCTAGC"

    df_default = find_slippage_sites(seq)
    df_explicit = find_slippage_sites(seq, mode="default")
    df_fast = find_slippage_sites(seq, mode="fast")

    assert df_default.shape[0] == df_explicit.shape[0] == df_fast.shape[0] == 1


def test_slippage_unknown_mode_raises():
    with pytest.raises(ValueError, match="Unknown slippage mode"):
        find_slippage_sites("ACGTACGTACGT", mode="nonexistent")
