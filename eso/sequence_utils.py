"""Small sequence/region helpers shared across detection and optimization."""

import pandas as pd

COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def reverse_complement_seq(seq):
    return ''.join(COMPLEMENT[x] for x in seq[::-1])


def add_backward_sites(df):
    """Duplicate each row of a {sequence, start, end} dataframe with its reverse complement."""
    df_forward = df.copy()
    df_backward = df.copy()
    df_backward.loc[:, 'sequence'] = df_backward.sequence.apply(reverse_complement_seq)
    return pd.concat([df_forward, df_backward], ignore_index=False)


def shorten_sequences(df):
    df_short = df.copy()
    df_short.loc[:, 'sequence'] = df_short.sequence.apply(lambda x: x[:-1])
    df_short.loc[:, 'end'] = df_short.end.apply(lambda x: x - 1)
    return df_short


def parse_region(region_string):
    """Parse a region string like "start_1-end_1,start_2-end_2,..." (1-indexed, inclusive)
    into a list of 0-indexed (start, end) tuples. Returns () for '' or 'None', 'error' if malformed.
    """
    if region_string in ('', 'None'):
        return ()

    region_list = region_string.split(',')
    try:
        return [
            (int(region.strip().split('-')[0]) - 1, int(region.strip().split('-')[1]))
            for region in region_list
        ]
    except (ValueError, IndexError):
        return 'error'
