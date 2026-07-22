"""Methylation motif detection via position-specific scoring matrices (PSSMs).

Scores every position of the sequence (forward and reverse complement)
against a set of methylation-enzyme recognition motifs and keeps the
best-scoring match per position above random-chance probability.
"""

import numpy as np
import pandas as pd
from Bio import motifs

from eso.sequence_utils import reverse_complement_seq

SITE_COLUMNS = [
    'start_index', 'end_index', 'matching_motif', 'PSSM_score',
    'actual_site', 'actual_site_reverse_conjugate',
]


def load_motifs(motifs_path):
    """Load PSSM motifs from a MEME-minimal-format file (e.g. topEnriched.*.meme.txt)."""
    with open(motifs_path, "r") as handle:
        return list(motifs.parse(handle, "minimal"))


def find_motif_sites(seq, num_sites, relevant_motifs):
    """Find the best-scoring motif match at every position of `seq` (forward or
    reverse complement, across all of `relevant_motifs`) that scores above
    random-chance probability (PSSM log-odds > 0). Returns at most `num_sites`
    rows, highest-scoring first.

    Builds one score matrix (2*len(relevant_motifs) rows x len(seq) columns,
    -inf where a motif doesn't reach that far) and reduces it with vectorized
    numpy calls, rather than a per-position/per-motif Python loop building an
    intermediate long-format dataframe with 2*len(relevant_motifs)*len(seq)
    rows - profiling showed that intermediate frame (and a df.apply(axis=1) to
    extract each site's sequence, now also replaced with plain zip()) was the
    dominant cost for realistic multi-motif sets on long sequences.
    """
    seq_len = len(seq)
    if not relevant_motifs or seq_len == 0:
        return pd.DataFrame(columns=SITE_COLUMNS)

    num_motifs = len(relevant_motifs)
    scores_matrix = np.full((2 * num_motifs, seq_len), -np.inf)
    motif_lengths = np.empty(num_motifs, dtype=int)
    motif_names = []

    for ii, motif in enumerate(relevant_motifs):
        length = len(motif)
        motif_lengths[ii] = length
        motif_names.append(motif.name)
        valid = seq_len - length + 1
        if valid <= 0:
            continue
        scores_matrix[2 * ii, :valid] = motif.pssm.calculate(seq)
        scores_matrix[2 * ii + 1, :valid] = motif.pssm.reverse_complement().calculate(seq)

    # ties (equal score at the same position) resolve to the lowest motif_number,
    # forward strand before backward - np.argmax returns the first max, and rows
    # are laid out (motif0-fwd, motif0-rev, motif1-fwd, ...), matching this
    # function's historical tie-breaking order.
    best_row = np.argmax(scores_matrix, axis=0)
    best_score = scores_matrix[best_row, np.arange(seq_len)]

    keep_mask = best_score > 0
    start_indices = np.nonzero(keep_mask)[0]
    if start_indices.size == 0:
        return pd.DataFrame(columns=SITE_COLUMNS)

    winning_rows = best_row[keep_mask]
    winning_scores = best_score[keep_mask]
    winning_motif_numbers = winning_rows // 2
    end_indices = start_indices + motif_lengths[winning_motif_numbers] - 1

    # highest-scoring first, matching this function's historical output order
    order = np.argsort(-winning_scores, kind='stable')
    if num_sites < order.size:
        order = order[:int(num_sites)]

    start_indices = start_indices[order]
    end_indices = end_indices[order]
    winning_scores = winning_scores[order]
    matching_motifs = [motif_names[m] for m in winning_motif_numbers[order]]

    actual_sites = [seq[s:e + 1] for s, e in zip(start_indices, end_indices)]
    actual_sites_rc = [reverse_complement_seq(site) for site in actual_sites]

    return pd.DataFrame({
        'start_index': start_indices,
        'end_index': end_indices,
        'matching_motif': matching_motifs,
        'PSSM_score': winning_scores,
        'actual_site': actual_sites,
        'actual_site_reverse_conjugate': actual_sites_rc,
    })
