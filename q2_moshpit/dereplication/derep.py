# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os.path
import shutil
from typing import List, Dict, Tuple

import pandas as pd
import skbio
from q2_types.feature_data import DNAFASTAFormat

from q2_types_genomics.feature_data import MAGSequencesDirFmt
from q2_types_genomics.per_sample_data import MultiMAGSequencesDirFmt


def find_similar_bins(distance_matrix: pd.DataFrame, threshold: float):
    sample_names = list(distance_matrix.index)
    similar_bins = {}
    for i, bin1 in enumerate(sample_names):
        for j, bin2 in enumerate(sample_names[i + 1 :]):
            distance = distance_matrix.loc[bin1, bin2]
            if distance <= threshold:
                if bin1 not in similar_bins:
                    similar_bins[bin1] = []
                similar_bins[bin1].append(bin2)
    return similar_bins


def get_bin_lengths(mags: MultiMAGSequencesDirFmt) -> pd.Series:
    bin_lengths = {}
    for path, seq in mags.sequences.iter_views(DNAFASTAFormat):
        tot = 0
        for _seq in skbio.io.read(str(seq.path), format="fasta"):
            tot += len(_seq)
        bin_lengths[str(path)] = tot
    ser = pd.Series(bin_lengths)
    ser.index = ser.index.map(lambda x: x.replace(".fasta", "").replace("/", "_"))
    return ser


def remap_bins(
    all_dereplicated: List[int], all_duplicates: List[List[str]], df: pd.DataFrame
) -> Dict[str, Tuple[str, str]]:
    final_bins = {}
    for i, dupl_bins in enumerate(all_duplicates):
        for bin in dupl_bins:
            final_bins[bin] = all_dereplicated[i]
    for bin in df.index:
        if bin not in final_bins:
            final_bins[bin] = bin

    # find unique bins and assign them a new ID
    unique_bins = list(set(final_bins.values()))
    new_ids = [f"mag{i}" for i in range(1, len(unique_bins) + 1)]

    # assign new IDs to old bins
    for old_id, new_id in zip(unique_bins, new_ids):
        for bin, new_bin in final_bins.items():
            if new_bin == old_id:
                final_bins[bin] = (old_id, new_id)

    return final_bins


def reassign_bins_to_samples(
    final_bins: Dict[str, Tuple[str, str]]
) -> Dict[str, Dict[str, int]]:
    all_samples = set(
        [sample_bin.rsplit("_", maxsplit=1)[0] for sample_bin in final_bins.keys()]
    )
    all_derep_mags = set([mag_id for _, (_, mag_id) in final_bins.items()])
    samples_to_bins = {key: {mag: 0 for mag in all_derep_mags} for key in all_samples}

    for samp_bin1, (samp_bin_2, mag_id) in final_bins.items():
        sample_id = samp_bin1.rsplit("_", maxsplit=1)[0]
        samples_to_bins[sample_id][mag_id] += 1

    return samples_to_bins


def write_unique_bins(
    all_bins: MultiMAGSequencesDirFmt, bins_remapped: Dict[str, Tuple[str, str]]
) -> MAGSequencesDirFmt:
    derep_bins = MAGSequencesDirFmt()
    manifest = all_bins.manifest.view(pd.DataFrame)
    for sample_bin, (old_bin_id, new_bin_id) in bins_remapped.items():
        sample, bin = old_bin_id.rsplit("_", maxsplit=1)
        dst_bin = os.path.join(str(derep_bins), f"{new_bin_id}.fasta")
        if not os.path.isfile(dst_bin):
            src_bin = manifest.loc[(sample, bin), "filename"]
            shutil.copy(src_bin, dst_bin)
    return derep_bins


def dereplicate_mags(
    mags: MultiMAGSequencesDirFmt, distance_matrix: skbio.DistanceMatrix
) -> (MAGSequencesDirFmt, pd.DataFrame):
    distances = distance_matrix.to_data_frame()
    duplicate_mags = find_similar_bins(distances, 0.2)
    all_duplicates = [[key, *value] for key, value in duplicate_mags.items()]
    mag_lengths = get_bin_lengths(mags)
    all_dereplicated = [mag_lengths[ids].idxmax() for ids in all_duplicates]

    final_bins = remap_bins(all_dereplicated, all_duplicates, distances)

    # generate dereplicated bins
    derep_bin_seqs = write_unique_bins(mags, final_bins)

    # generate a presence-absence table
    derep_bins_per_sample = reassign_bins_to_samples(final_bins)
    presence_absence = pd.DataFrame.from_records(derep_bins_per_sample).T
    presence_absence.index.name = "sample-id"
    presence_absence.sort_index(inplace=True)
    presence_absence = presence_absence[sorted(presence_absence.columns)]

    return derep_bin_seqs, presence_absence
