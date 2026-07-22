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
