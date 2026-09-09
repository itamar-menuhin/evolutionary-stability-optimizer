"""Property-based (fuzz) test for eso.detection.staubility_variant's
_repeated_kmer_positions, checked against an independent re.finditer-based
reference implementation.

Gap found while auditing: the dict-scan replacement for sklearn's
CountVectorizer + re.finditer (see git commit d777761, "remove eso's only
scikit-learn dependency from staubility_variant") was verified during
development against real sklearn output across "300+ random trials" per its
own commit message, but that verification was never preserved as a permanent
test - only two hand-picked example regressions (the homopolymer/tandem-
repeat merge-collapse cases) made it into the checked-in suite. Re-adding
sklearn here to re-run that exact comparison would reintroduce the very
dependency this change was meant to remove; re.finditer is the more apt
reference anyway - it's the actual stdlib mechanism _repeated_kmer_positions'
own docstring says it's reproducing the non-overlapping-position semantics
of, not just a stand-in for sklearn's counting behavior.
"""

import re

from hypothesis import given, settings, strategies as st

from eso.detection.staubility_variant import _repeated_kmer_positions

_dna_sequence = st.text(alphabet='ACGT', min_size=0, max_size=200)
_kmer_length = st.integers(min_value=1, max_value=16)


def _reference_repeated_kmer_positions(seq, k):
    """Every k-mer with 2+ overlapping occurrences, each with its
    NON-overlapping positions via re.finditer - exactly the semantics
    _repeated_kmer_positions' own docstring says it reproduces.
    """
    positions_by_kmer = {}
    for i in range(len(seq) - k + 1):
        positions_by_kmer.setdefault(seq[i:i + k], []).append(i)

    result = {}
    for kmer, positions in positions_by_kmer.items():
        if len(positions) < 2:
            continue
        result[kmer] = [m.start() for m in re.finditer(re.escape(kmer), seq)]
    return result


@settings(deadline=None)
@given(seq=_dna_sequence, k=_kmer_length)
def test_matches_re_finditer_reference(seq, k):
    assert _repeated_kmer_positions(seq, k) == _reference_repeated_kmer_positions(seq, k)


def test_matches_re_finditer_reference_on_targeted_edge_cases():
    tests = {
        "homopolymer 40, k=16": ("A" * 40, 16),
        "homopolymer exactly 2x k": ("A" * 32, 16),
        "tandem repeat AT*30, k=2": ("AT" * 30, 2),
        "tandem repeat AT*30, k=16": ("AT" * 30, 16),
        "k longer than sequence": ("ACGTACGT", 16),
        "k equal to sequence length": ("ACGTACGTACGTACGT", 16),
        "empty sequence": ("", 4),
        "no repeats": ("ACGTAGCTTGACCTGAAGCTAGCA", 4),
    }
    for name, (seq, k) in tests.items():
        assert _repeated_kmer_positions(seq, k) == _reference_repeated_kmer_positions(seq, k), name
