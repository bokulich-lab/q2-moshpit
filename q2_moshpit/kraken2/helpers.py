# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os

import pandas as pd

from q2_types_genomics.kraken2 import (
    Kraken2ReportDirectoryFormat,
    Kraken2OutputDirectoryFormat
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


def _collate_kraken_tsvs(kraken2_results, kraken_type, output):
    collated_sample_outputs = {}

    for kraken2_result in kraken2_results:
        for fp in kraken2_result.path.iterdir():
            # fp basename will be "sample_id.output/report.txt/tsv" we just
            # want the id
            sample_id = os.path.basename(fp).split('.')[0]
            df = pd.read_csv(fp, sep='\t', header=None)

            if sample_id in collated_sample_outputs:
                collated_sample_outputs[sample_id].append(df)
            else:
                collated_sample_outputs[sample_id] = df

    for sample_id, df in collated_sample_outputs.items():
        df.to_csv(output.path / f'{sample_id}.{kraken_type}.txt', sep='\t')
