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

import pandas as pd
import pytest

from eso.optimize import optimization_engine


def test_translation_is_preserved():
    seq = "ATG" + "TTT" * 10 + "TAA"  # all-Phe, a real slippage/CAI target
    final_seq, _, _ = optimization_engine(seq, organism_name="kompas")

    assert len(final_seq) == len(seq)
    assert final_seq[:3] == "ATG"
    assert final_seq[-3:] in ("TAA", "TAG", "TGA")


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
    seq = "ATG" + "TTT" * 5 + "TAA"
    final_seq, obj_description, _ = optimization_engine(seq, organism_name="not_a_real_organism_xyz")

    assert len(final_seq) == len(seq)
    assert final_seq[:3] == "ATG"


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
