"""Tests for eso.detection.common_motifs.find_sigma70_promoter_pairs -
strand-aware, correctly-spaced sigma70 -35/-10 hexamer pair detection (a
real cryptic promoter needs both hexamers at ~17±1bp apart, which
find_motif_sites' per-position independent scoring has no way to express -
see common_motifs.py's own docstring).
"""

from Bio.Seq import Seq

from eso.detection.common_motifs import find_sigma70_promoter_pairs, load_common_motifs
from eso.detection.methylation import find_motif_sites

MINUS35 = "TTGACA"
MINUS10 = "TATAAT"


def test_finds_a_correctly_spaced_pair_on_the_forward_strand():
    spacer = "A" * 17  # -35 end to -10 start: exactly 17nt, within 16-18
    seq = MINUS35 + spacer + MINUS10

    result = find_sigma70_promoter_pairs(seq)

    assert len(result) == 1
    row = result.iloc[0]
    assert row.strand == '+'
    assert (row.minus35_start, row.minus35_end) == (0, 5)
    assert (row.minus10_start, row.minus10_end) == (23, 28)
    assert row.spacing == 17


def test_finds_a_correctly_spaced_pair_on_the_reverse_strand():
    spacer = "A" * 17
    forward_seq = MINUS35 + spacer + MINUS10
    seq = str(Seq(forward_seq).reverse_complement())

    result = find_sigma70_promoter_pairs(seq)

    assert len(result) == 1
    assert result.iloc[0].strand == '-'
    assert result.iloc[0].spacing == 17


def test_no_pair_when_spacing_is_wrong():
    # 10nt spacer - well outside the 16-18nt window, no real promoter here.
    seq = MINUS35 + ("A" * 10) + MINUS10

    result = find_sigma70_promoter_pairs(seq)

    assert result.empty


def test_no_pair_when_hexamers_are_on_opposite_strands():
    # -35 on the forward strand, -10 only present via its reverse complement
    # (i.e. really on the opposite strand) - not a real promoter, must not
    # be reported as a pair just because both hexamers appear somewhere.
    spacer = "A" * 17
    seq = MINUS35 + spacer + str(Seq(MINUS10).reverse_complement())

    result = find_sigma70_promoter_pairs(seq)

    assert result.empty


def test_isolated_single_hexamer_yields_no_pair():
    seq = MINUS35 + ("A" * 30)  # no -10 box anywhere

    result = find_sigma70_promoter_pairs(seq)

    assert result.empty


def test_multiple_pairs_ranked_by_combined_score_highest_first():
    spacer = "A" * 17
    # two independent, well-separated candidate pairs in one sequence
    seq = (MINUS35 + spacer + MINUS10) + ("A" * 40) + (MINUS35 + spacer + MINUS10)

    result = find_sigma70_promoter_pairs(seq)

    assert len(result) == 2
    # sorted highest combined_score first - both are identical consensus
    # sequences here, so scores tie, but the ordering contract itself
    # (descending, not insertion order) is what's being checked.
    assert list(result.combined_score) == sorted(result.combined_score, reverse=True)


def test_load_common_motifs_pseudocount_loosens_matching():
    # Regression test for a real gap: COMMON_MOTIFS' bundled motifs had no
    # way to adjust match tolerance without editing source, even though the
    # underlying primitive (motif_from_consensus) already supports it via
    # `pseudocount` - load_common_motifs now passes it through. "GACC" is
    # "GATC" (the dam motif) with one mismatch (T->C); a higher pseudocount
    # measurably widens what counts as a match (more marginal positions
    # elsewhere in the sequence cross the scoring threshold too), confirming
    # the knob actually reaches common_motifs, not just custom ones.
    seq = "ATG" + "GACC" + "AAAAAAAAAAAAAAAAAAAA"

    strict = load_common_motifs(['dam'], pseudocount=1)
    loose = load_common_motifs(['dam'], pseudocount=50)

    strict_hits = find_motif_sites(seq, num_sites=50, relevant_motifs=strict)
    loose_hits = find_motif_sites(seq, num_sites=50, relevant_motifs=loose)

    assert len(loose_hits) > len(strict_hits)


def test_load_common_motifs_default_pseudocount_is_unchanged():
    # the default behavior (no pseudocount passed) must stay exactly what
    # it was before this parameter existed.
    default = load_common_motifs(['dam'])
    explicit = load_common_motifs(['dam'], pseudocount=1)
    assert default[0].counts == explicit[0].counts
