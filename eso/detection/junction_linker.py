"""Hypermutable-site detection applied to a new domain: checking whether
inserting a linker between a target and host sequence creates new slippage
or recombination hotspots at the junction boundaries.

Developed inside the STABLES project (linker_suspect_utils.py) as a
lightweight, junction-window-scoped variant of eso.detection.slippage /
eso.detection.recombination - it checks only the region spanning the
target-linker-host boundary rather than scanning the full sequence.
"""

from __future__ import annotations

from typing import Literal

Mode = Literal["linker", "target"]


def find_slippage_1(target_seq: str, linker_seq: str, host_seq: str, mode: Mode = "linker") -> str:
    """Check for a run of 4+ identical nucleotides at the junction (or within `target_seq` if mode='target')."""
    test_seq = target_seq[-3:] + linker_seq + host_seq[:3]
    if mode == "target":
        test_seq = target_seq

    subseqs = ["AAAA", "CCCC", "GGGG", "TTTT"]
    found = [subseq for subseq in subseqs if test_seq.find(subseq) > -1]
    return ",".join(found)


def find_slippage_l(
    target_seq: str,
    linker_seq: str,
    host_seq: str,
    length: int,
    mode: Mode = "linker",
) -> str:
    """Check for a base unit of `length` repeated 3x at the junction."""
    border_add = 3 * length - 1
    test_seq = target_seq[-border_add:] + linker_seq + host_seq[:border_add]
    if mode == "target":
        test_seq = target_seq

    all_substrings = [test_seq[index:index + length] for index in range(len(test_seq))]
    all_substrings = [value for value in all_substrings if len(value) == length]
    all_substrings = sorted(set(all_substrings))
    subseqs = [value * 3 for value in all_substrings]

    found = [subseq for subseq in subseqs if test_seq.find(subseq) > -1]
    return ",".join(found)


def find_recombination(
    target_seq: str,
    linker_seq: str,
    host_seq: str,
    length: int,
    mode: Mode = "linker",
) -> str:
    """Check whether any `length`-mer newly introduced by the junction already
    exists elsewhere in the target/host (mode='linker') or linker/host (mode='target').
    """
    border_add = length - 1
    test_seq = target_seq[-border_add:] + linker_seq + host_seq[:border_add]
    if mode == "target":
        test_seq = target_seq

    all_substrings = [test_seq[index:index + length] for index in range(len(test_seq))]
    all_substrings = [value for value in all_substrings if len(value) == length]
    subseqs = sorted(set(all_substrings))

    found = []
    if mode == "linker":
        for subseq in subseqs:
            if target_seq.find(subseq) > -1:
                found.append(subseq)
            if host_seq.find(subseq) > -1:
                found.append(subseq)
    else:
        for subseq in subseqs:
            if linker_seq.find(subseq) > -1:
                found.append(subseq)
            if host_seq.find(subseq) > -1:
                found.append(subseq)

    return ",".join(found)


def find_suspect(
    target_seq: str,
    linker_seq: str,
    host_seq: str,
    mode: Mode = "linker",
    max_l: int = 12,
) -> str:
    """Run all junction-hotspot checks (slippage lengths 1..max_l-1, recombination at max_l)."""
    suspects = [find_slippage_1(target_seq, linker_seq, host_seq, mode)]
    for length in range(2, max_l):
        suspects.append(find_slippage_l(target_seq, linker_seq, host_seq, length, mode))

    suspects.append(find_recombination(target_seq, linker_seq, host_seq, max_l, mode))
    suspects = [value for value in suspects if value != ""]
    return ",".join(suspects)
