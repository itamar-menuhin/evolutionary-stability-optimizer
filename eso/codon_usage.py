"""Codon usage bias (CUB) tables for hosts not covered by python-codon-tables."""

from importlib import resources

import pandas as pd
from Bio import SeqUtils


def cub_c1():
    codon_usage_table = {
        'A': {'GCA': 0.13718, 'GCC': 0.432097, 'GCG': 0.264336, 'GCT': 0.166387},
        '*': {'TAA': 0.265377, 'TAG': 0.181818, 'TGA': 0.552805}, 'W': {'TGG': 1},
        'R': {'AGA': 0.0973128, 'AGG': 0.153007, 'CGA': 0.128389, 'CGC': 0.304978, 'CGG': 0.209737,
              'CGT': 0.106575}, 'N': {'AAC': 0.792265, 'AAT': 0.207735},
        'D': {'GAC': 0.718222, 'GAT': 0.281778}, 'C': {'TGC': 0.710465, 'TGT': 0.289535},
        'Q': {'CAA': 0.361039, 'CAG': 0.638961}, 'E': {'GAA': 0.250212, 'GAG': 0.749788},
        'G': {'GGA': 0.142656, 'GGC': 0.506165, 'GGG': 0.184803, 'GGT': 0.166376},
        'H': {'CAC': 0.650846, 'CAT': 0.349154},
        'I': {'ATA': 0.0923118, 'ATC': 0.673554, 'ATT': 0.234134},
        'L': {'CTA': 0.0641845, 'CTC': 0.345609, 'CTG': 0.309073, 'CTT': 0.132086, 'TTA': 0.030553,
              'TTG': 0.118494}, 'K': {'AAA': 0.180592, 'AAG': 0.819408}, 'M': {'ATG': 1},
        'F': {'TTC': 0.672764, 'TTT': 0.327236},
        'P': {'CCA': 0.167124, 'CCC': 0.336275, 'CCG': 0.326667, 'CCT': 0.169935},
        'Y': {'TAC': 0.757778, 'TAT': 0.242222},
        'S': {'AGC': 0.234258, 'AGT': 0.0671675, 'TCA': 0.103458, 'TCC': 0.220449, 'TCG': 0.247293,
              'TCT': 0.127375},
        'T': {'ACA': 0.154802, 'ACC': 0.406355, 'ACG': 0.307613, 'ACT': 0.131231},
        'V': {'GTA': 0.0871417, 'GTC': 0.460197, 'GTG': 0.289844, 'GTT': 0.162817},
    }
    return codon_usage_table


def cub_kompas():
    """Codon usage table for Komagataella phaffii (Pichia pastoris)."""
    cub_full = {
        'Ala': {'GCA': 0.275098, 'GCC': 0.244931, 'GCG': 0.0786548, 'GCT': 0.401316},
        'Arg': {'AGA': 0.455639, 'AGG': 0.181255, 'CGA': 0.119597, 'CGC': 0.0506672, 'CGG': 0.0522151,
                'CGT': 0.140626},
        'Asn': {'AAC': 0.465812, 'AAT': 0.534188},
        'Asp': {'GAC': 0.382762, 'GAT': 0.617238},
        'Cys': {'TGC': 0.377246, 'TGT': 0.622754},
        'Gln': {'CAA': 0.603655, 'CAG': 0.396345},
        'Glu': {'GAA': 0.594953, 'GAG': 0.405047},
        'Gly': {'GGA': 0.364371, 'GGC': 0.155831, 'GGG': 0.12228, 'GGT': 0.357518},
        'His': {'CAC': 0.382548, 'CAT': 0.617452},
        'Ile': {'ATA': 0.236221, 'ATC': 0.297661, 'ATT': 0.466118},
        'Leu': {'CTA': 0.122977, 'CTC': 0.0842546, 'CTG': 0.152862, 'CTT': 0.170539, 'TTA': 0.178334,
                'TTG': 0.291033},
        'Lys': {'AAA': 0.519264, 'AAG': 0.480736},
        'Met': {'ATG': 1},
        'Phe': {'TTC': 0.42216, 'TTT': 0.57784},
        'Pro': {'CCA': 0.378486, 'CCC': 0.180179, 'CCG': 0.102394, 'CCT': 0.338942},
        'Ser': {'AGC': 0.104434, 'AGT': 0.158704, 'TCA': 0.208354, 'TCC': 0.174084, 'TCG': 0.0945539,
                'TCT': 0.259871},
        'Thr': {'ACA': 0.277392, 'ACC': 0.237645, 'ACG': 0.122946, 'ACT': 0.362017},
        'Trp': {'TGG': 1},
        'Tyr': {'TAC': 0.48627, 'TAT': 0.51373},
        'Val': {'GTA': 0.178012, 'GTC': 0.217067, 'GTG': 0.216627, 'GTT': 0.388294},
        'END': {'TAA': 0.399841, 'TAG': 0.339428, 'TGA': 0.260731},
    }
    # SeqUtils.seq1('END') returns 'X' (undefined amino acid), not '*' (stop) -
    # without this special case, the stop-codon frequencies silently ended up
    # filed under the wrong key and were never used for stop-codon scoring.
    return {('*' if x == 'END' else SeqUtils.seq1(x)): cub_full[x] for x in cub_full}


def _load_bundled_csv_cub(data_filename):
    with resources.files("eso.data").joinpath(data_filename).open("r") as handle:
        df = pd.read_csv(handle)

    codon_usage_table = df.groupby('aa').apply(
        lambda x: x.set_index('codon')['freq_within_aa'].to_dict()
    ).to_dict()

    if '*' not in codon_usage_table:
        codon_usage_table['*'] = {'TAA': 0.33, 'TAG': 0.33, 'TGA': 0.34}

    return codon_usage_table


def cub_human_antibody_heavy_chain():
    """Codon usage table for human antibody heavy chain, from the iGEM 2025 dataset."""
    return _load_bundled_csv_cub("human-antibody-heavy-chain-codon-frequencies.csv")


def cub_human_antibody_light_chain():
    """Codon usage table for human antibody light chain, from the iGEM 2025 dataset."""
    return _load_bundled_csv_cub("human-antibody-light-chain-codon-frequencies.csv")


CODON_USAGE_TABLES = {
    'C1': cub_c1,
    'kompas': cub_kompas,
    'human_antibody_heavy_chain': cub_human_antibody_heavy_chain,
    'human_antibody_light_chain': cub_human_antibody_light_chain,
}
