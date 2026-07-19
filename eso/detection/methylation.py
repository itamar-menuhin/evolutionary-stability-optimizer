"""Methylation motif detection via position-specific scoring matrices (PSSMs).

Scores every position of the sequence (forward and reverse complement)
against a set of methylation-enzyme recognition motifs and keeps the
best-scoring match per position above random-chance probability.
"""

import pandas as pd
from Bio import motifs

from eso.sequence_utils import reverse_complement_seq


def load_motifs(motifs_path):
    """Load PSSM motifs from a MEME-minimal-format file (e.g. topEnriched.*.meme.txt)."""
    with open(motifs_path, "r") as handle:
        return list(motifs.parse(handle, "minimal"))


def find_motif_sites(seq, num_sites, relevant_motifs):
    site_scores_list = []
    motifs_list = []

    for ii, motif in enumerate(relevant_motifs):
        motif_scores = list(motif.pssm.calculate(seq))
        for jj, score in enumerate(motif_scores):
            site_scores_list.append((jj, ii, score, 'forward'))

        motif_rev_scores = list(motif.pssm.reverse_complement().calculate(seq))
        for jj, score in enumerate(motif_rev_scores):
            site_scores_list.append((jj, ii, score, 'backward'))

        motifs_list.append((ii, motif.name, motif.__len__()))

    df_sites = pd.DataFrame.from_records(
        data=site_scores_list, columns=['start_index', 'motif_number', 'PSSM_score', 'strand'])

    # keep only the best match per position, above random-chance probability
    df_sites = df_sites.sort_values('PSSM_score', ascending=False)
    df_sites = df_sites.drop_duplicates(subset=['start_index'])
    df_sites = df_sites[df_sites.PSSM_score > 0]

    df_motifs = pd.DataFrame.from_records(motifs_list, columns=['motif_number', 'matching_motif', 'motif_length'])
    df_sites = df_sites.merge(df_motifs, on='motif_number')

    df_sites.loc[:, 'end_index'] = df_sites.start_index + df_sites.motif_length - 1

    if num_sites < df_sites.shape[0]:
        df_sites = df_sites.head(num_sites)

    df_sites.loc[:, 'actual_site'] = df_sites.apply(lambda x: seq[x.start_index:(x.end_index + 1)], axis=1)
    df_sites.loc[:, 'actual_site_reverse_conjugate'] = df_sites.actual_site.apply(reverse_complement_seq)

    return df_sites[
        ['start_index', 'end_index', 'matching_motif', 'PSSM_score', 'actual_site', 'actual_site_reverse_conjugate']
    ]
