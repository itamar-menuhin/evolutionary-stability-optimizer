"""Alternate recombination/slippage detectors, developed independently inside the
STABLES project (Staubility_Code_shimshi.py) using a different algorithmic
approach than eso.detection.recombination/slippage:

- Recombination: ngram-counting (CountVectorizer, 16-mers) + exact re-match,
  vs. the Levenshtein-neighbor generation approach in eso.detection.recombination.
- Slippage: per-frameshift scan for back-to-back repeats, vs. the
  find-all-occurrences approach in eso.detection.slippage.

This module also carried a methylation-motif scorer (site_motif_grader/
calc_max_site) at one point; removed after comparison against
eso.detection.methylation found it scored candidates by raw sequence
probability with no background correction, which measurably disagrees with
(not just runs slower than) eso.detection.methylation's background-corrected
approach whenever a motif's background isn't uniform - see
docs/detector-comparisons.md for the full writeup, kept as a historical
record even though the code itself is gone.

Kept as a separate module rather than merged into the primary detectors -
these are two independently-arrived-at implementations of the same EFM
Calculator concept and haven't yet been reconciled into one canonical API.
"""

import re

import numpy as np
import pandas as pd
from sklearn.feature_extraction.text import CountVectorizer

from eso.detection._overlap import collapse_overlapping_intervals


def genome_cutter(start, end, seq):
    return seq[start:end]


def find_recombination_sites(seq, num_sites=np.inf):
    """ngram-counting variant of recombination (RMD) hotspot detection."""
    vectorizer = CountVectorizer(analyzer='char_wb', ngram_range=(16, 16))
    counter = vectorizer.fit_transform([seq]).toarray()

    sites_recombination = list(np.where(counter > 1)[1])
    empty_columns = [
        'start_1', 'end_1', 'sequence', 'start_2', 'end_2',
        'location_delta', 'site_length', 'log10_prob_recombination_ecoli', 'sequence_number',
    ]
    if not sites_recombination:
        return pd.DataFrame(columns=empty_columns)

    all_sites = vectorizer.get_feature_names_out()

    suspect_recombination = []
    for site in sites_recombination:
        curr_seq = all_sites[site]
        list_regions = [match.span() for match in re.finditer(curr_seq.upper(), seq)]
        suspect_recombination.extend(list_regions)

    suspect_recombination = sorted(suspect_recombination)
    df_recombination = pd.DataFrame(suspect_recombination, columns=['start', 'end'])

    # merge back-to-back 16-mer matches (from repeats longer than 16nt) into one site
    df_recombination.loc[:, 'start_delta'] = df_recombination['start'] - df_recombination['start'].shift()
    df_recombination.loc[:, 'end_delta'] = df_recombination['end'].shift(-1) - df_recombination['end']
    df_recombination = df_recombination[(df_recombination.start_delta != 1.0) | (df_recombination.end_delta != 1.0)]

    df_recombination.loc[(df_recombination.end_delta == 1.0), 'end'] = None
    # `df.loc[:, 'end'] = ...` would silently keep the column's existing float64
    # dtype (from the None assignment above) even after .astype(int) below;
    # only whole-column reassignment (df['end'] = ...) actually changes dtype.
    df_recombination['end'] = df_recombination.loc[:, 'end'].bfill().astype(int) - 1
    df_recombination = df_recombination[df_recombination.start_delta != 1.0][['start', 'end']]

    df_recombination.loc[:, 'sequence'] = df_recombination.apply(
        lambda x: genome_cutter(x['start'], x['end'], seq), axis=1)

    df_recombination = df_recombination.merge(df_recombination, on='sequence', suffixes=('_1', '_2'))
    df_recombination = df_recombination[df_recombination.end_1 < df_recombination.start_2]

    df_recombination.loc[:, 'location_delta'] = df_recombination.start_2 - df_recombination.end_1
    df_recombination.loc[:, 'site_length'] = df_recombination.end_1 - df_recombination.start_1

    # a=5.8 per Oliveira et al. 2008, Table 3, recA+ row (see
    # eso.detection.recombination.calc_recombination_score's docstring for
    # the full story, including why the EFM Calculator's own 8.8 is a
    # mix-up with a different row/parameter of that same table).
    a, b, c, alpha = 5.8, 1465.6, 0, 29
    base = a + df_recombination['location_delta']
    exponent = -1 * alpha / df_recombination['site_length']
    scale = df_recombination['site_length'] / (
        1 + b * df_recombination['site_length'] + c * df_recombination['location_delta'])
    df_recombination.loc[:, 'log10_prob_recombination_ecoli'] = np.log10((base ** exponent) * scale)

    df_recombination = df_recombination.sort_values('log10_prob_recombination_ecoli', ascending=False)

    if num_sites < np.inf:
        df_recombination = df_recombination.head(int(num_sites))

    for col in df_recombination:
        if df_recombination[col].isnull().all():
            del df_recombination[col]

    return df_recombination


def find_slippage_sites_length_l(seq, length):
    """Per-frameshift scan for a repeated base unit of size `length`
    (repeated 3+ times, or 12+ times for length=1 - see
    eso.detection.slippage's module docstring for why 12: below that, an
    independently-run length-2 detection of the same physical region always
    outscores it, so its result is always discarded downstream).
    """
    from textwrap import wrap

    slippage_sites = []
    homopolymer_window = 12

    for frameshift in range(length):
        curr_seq = seq[frameshift:]
        curr_seq_split = wrap(curr_seq, length)

        end_of_range = len(curr_seq_split) - 2

        for ii in range(end_of_range):
            # `length == 1` must be checked FIRST so Python's `and`
            # short-circuits before touching curr_seq_split[ii + 1/2] when
            # length > 1 (needed for is_followed2 below); is_followed1 itself
            # uses slicing, which never raises IndexError even past the list
            # end, so it needs no equivalent bounds bug to worry about - a
            # slice shorter than `homopolymer_window` just fails the length
            # check instead of crashing.
            is_followed2 = (
                length > 1
                and curr_seq_split[ii] == curr_seq_split[ii + 1]
                and curr_seq_split[ii] == curr_seq_split[ii + 2]
            )
            window = curr_seq_split[ii:ii + homopolymer_window]
            is_followed1 = length == 1 and len(window) == homopolymer_window and len(set(window)) == 1

            if is_followed2:
                slippage_sites.append((frameshift + ii * length, frameshift + length * (ii + 3)))
            if is_followed1:
                slippage_sites.append((ii, ii + homopolymer_window))

    if not slippage_sites:
        return pd.DataFrame(columns=['start', 'end', 'sequence', 'length_base_unit'])

    df_slippage = pd.DataFrame(sorted(slippage_sites), columns=['start', 'end'])

    df_slippage.loc[:, 'start_delta'] = df_slippage['start'] - df_slippage['start'].shift()
    df_slippage.loc[:, 'end_delta'] = df_slippage['end'].shift(-1) - df_slippage['end']
    df_slippage = df_slippage[(df_slippage.start_delta != 1.0) | (df_slippage.end_delta != 1.0)]
    df_slippage.loc[(df_slippage.end_delta == 1.0), 'end'] = None
    # see the identical fix in find_recombination_sites above: whole-column
    # reassignment is required for .astype(int) to actually stick here
    df_slippage['end'] = df_slippage.loc[:, 'end'].bfill().astype(int)
    df_slippage = df_slippage[df_slippage.start_delta != 1.0][['start', 'end']]

    df_slippage.loc[:, 'end'] = (
        ((df_slippage['end'] - df_slippage['start']) / length).astype(int) * length
    ) + df_slippage['start']

    df_slippage.loc[:, 'sequence'] = df_slippage.apply(
        lambda x: genome_cutter(x['start'], x['end'], seq), axis=1)
    df_slippage.loc[:, 'length_base_unit'] = length

    return df_slippage


def find_slippage_sites(seq, num_sites=np.inf):
    """ngram/frameshift-scan variant of slippage (SSR) hotspot detection."""
    slippage_sites_list = [find_slippage_sites_length_l(seq, length) for length in range(1, 16)]
    df_slippage = pd.concat(slippage_sites_list, ignore_index=True)[['start', 'end', 'length_base_unit', 'sequence']]

    df_slippage.loc[:, 'num_base_units'] = (
        df_slippage.sequence.apply(len) / df_slippage.length_base_unit
    ).astype(int)

    df_slippage.loc[:, 'log10_prob_slippage_ecoli'] = -4.749 + 0.063 * df_slippage['num_base_units']
    df_slippage.loc[df_slippage.length_base_unit == 1, 'log10_prob_slippage_ecoli'] = (
        -12.9 + 0.729 * df_slippage['num_base_units']
    )

    # this filter was applied by the STABLES caller (STABLES_full_code's
    # optimization pipeline), not inside find_slippage_sites itself - lost
    # when this function was extracted in isolation. Restored here so the
    # function is self-contained and consistent with its sibling
    # find_recombination_sites above, which does filter internally.
    df_slippage = df_slippage[df_slippage.log10_prob_slippage_ecoli > -9]

    # see eso.detection.slippage.find_slippage_sites for why exact-start dedup
    # isn't enough - phase-shifted/differently-lengthed detections of the same
    # physical repeat need an overlap-based collapse instead.
    df_slippage = collapse_overlapping_intervals(df_slippage, score_col='log10_prob_slippage_ecoli')
    df_slippage = df_slippage.sort_values(['log10_prob_slippage_ecoli', 'length_base_unit'], ascending=[False, False])

    if num_sites < np.inf:
        df_slippage = df_slippage.head(int(num_sites))

    return df_slippage


def suspect_site_extractor(seq, num_sites=np.inf, extension=''):
    """Run both recombination and slippage detection (this module's ngram/frameshift variants)."""
    return {
        'df_recombination' + extension: find_recombination_sites(seq, num_sites),
        'df_slippage' + extension: find_slippage_sites(seq, num_sites),
    }
