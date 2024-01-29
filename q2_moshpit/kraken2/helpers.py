# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
from pathlib import Path

import pandas as pd
from q2_types.sample_data import SampleData
from qiime2.core.type import Properties

from q2_types_genomics.kraken2 import (
    Kraken2ReportDirectoryFormat,
    Kraken2OutputDirectoryFormat, Kraken2Reports
)


def collate_kraken2_reports(kraken2_reports: Kraken2ReportDirectoryFormat) \
        -> Kraken2ReportDirectoryFormat:
    collated_kraken2_reports = Kraken2ReportDirectoryFormat()
    _collate_kraken_tsvs(kraken2_reports, 'report', collated_kraken2_reports)
    return collated_kraken2_reports


def collate_kraken2_outputs(kraken2_outputs: Kraken2OutputDirectoryFormat) \
        -> Kraken2OutputDirectoryFormat:
    collated_kraken2_outputs = Kraken2OutputDirectoryFormat()
    _collate_kraken_tsvs(kraken2_outputs, 'output', collated_kraken2_outputs)
    return collated_kraken2_outputs


def _read_and_collect_df(fp, sample_id, results, mag_id=None):
    """Reads reports into DataFrames and returns a
        DF collection as a dictionary.

    There are two use cases:
    1. Reports are stored per-sample, in which case the
        mag_id needs to be provided (as it is included
        in the filename). -> returns a dictionary:
        {sample_id: {mag_id: df}}
    2. Reports are stored in a single directory, in which case
        the mag_id is not provided. -> returns a dictionary:
        {sample_id: df}
    """
    df = pd.read_csv(fp, sep='\t', header=None)
    if mag_id is not None:
        if mag_id in results[sample_id]:
            results[sample_id][mag_id].append(
                df, ignore_index=True)
        else:
            results[sample_id][mag_id] = df
    else:
        if sample_id in results:
            results[sample_id].append(df, ignore_index=True)
        else:
            results[sample_id] = df


def _handle_directory(filepath, results):
    """Handles reports/outputs stored per-directory."""
    sample_id = os.path.basename(filepath)
    if sample_id not in results:
        results[sample_id] = {}
    for mp in Path(filepath).iterdir():
        mag_id = os.path.basename(mp).split('.')[0]
        _read_and_collect_df(mp, sample_id, results, mag_id)


def _collate_kraken_tsvs(kraken2_results, kraken_type, output):
    collated_sample_outputs = {}

    for kraken2_result in kraken2_results:
        for fp in kraken2_result.path.iterdir():
            # SampleData[MAGs] were used to generate the reports
            if os.path.isdir(fp):
                _handle_directory(fp, collated_sample_outputs)
            else:
                sample_id = os.path.basename(fp).split('.')[0]
                _read_and_collect_df(fp, sample_id, collated_sample_outputs)

    for sample_id, val in collated_sample_outputs.items():
        # SampleData[MAGs] were used to generate the reports
        if isinstance(val, dict):
            for mag_id, df in val.items():
                os.makedirs(output.path / sample_id, exist_ok=True)
                df.to_csv(
                    output.path / sample_id / f'{mag_id}.{kraken_type}.txt',
                    sep='\t', index=False, header=False
                )
        else:
            val.to_csv(
                output.path / f'{sample_id}.{kraken_type}.txt', sep='\t',
                index=False, header=False
            )
