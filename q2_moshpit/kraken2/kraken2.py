# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import subprocess
import tempfile
from copy import deepcopy
from functools import reduce

import pandas as pd
from q2_moshpit._utils import _process_common_input_params, run_command
from q2_types_genomics.per_sample_data import MultiMAGSequencesDirFmt
from q2_sapienns._metaphlan import metaphlan_taxon

from q2_moshpit.kraken2.utils import _process_kraken2_arg


def _parse_kraken2_report(report_fp: str, sample_name: str, bin_name: str) -> pd.DataFrame:
    with open(report_fp, 'r') as r:
        lines = r.readlines()

    all_taxonomies = {}
    for line in lines:
        line = line.split('\t')
        taxonomy, count = line[0], int(line[1].strip())
        all_taxonomies[taxonomy] = {f'{sample_name}/{bin_name}': count}

    return pd.DataFrame.from_dict(all_taxonomies, orient='index')


def _process_kraken2_reports(reports: dict) -> pd.DataFrame:
    # get abundances per sample
    abundances = []
    for _sample, bins in reports.items():
        for _bin, report_fp in bins.items():
            abundances.append(_parse_kraken2_report(report_fp, _sample, _bin))

    # combine all samples
    df_merged = reduce(
        lambda left, right: pd.merge(
            left, right, left_index=True, right_index=True, how='outer'
        ), abundances
    )

    # find all levels from all samples (to assign NCBI tax IDs later)
    all_levels = []
    for taxonomy in df_merged.index.tolist():
        for level in taxonomy.split('|'):
            all_levels.append(level) if level not in all_levels else False
    levels_with_ids = {
        level: _id for _id, level in enumerate(sorted(all_levels))
    }

    # assign NCBI IDs (fake for now)
    all_taxonomies = []
    for taxonomy in df_merged.index:
        _id = list(map(lambda x: str(levels_with_ids[x]), taxonomy.split('|')))
        all_taxonomies.append('|'.join(_id))

    df_merged['NCBI_tax_id'] = all_taxonomies
    df_merged.index.name = 'feature-id'

    return df_merged


def _classify_kraken(manifest, common_args) -> (pd.DataFrame, pd.DataFrame):
    base_cmd = ["kraken2", *common_args]
    kraken2_reports = {}

    with tempfile.TemporaryDirectory() as tmp_dir:
        try:
            for (_sample, _bin), fn in manifest.itertuples():
                if _sample not in kraken2_reports:
                    kraken2_reports[_sample] = {}
                cmd = deepcopy(base_cmd)
                bin_dir = os.path.join(tmp_dir, _sample, _bin)
                os.makedirs(bin_dir, exist_ok=True)
                report_fp = os.path.join(bin_dir, 'report.txt')
                output_fp = os.path.join(bin_dir, 'output.txt')
                cmd.extend([
                    '--use-mpa-style', '--report', report_fp,
                    '--output', output_fp, fn
                ])
                run_command(cmd=cmd, verbose=True)
                kraken2_reports[_sample].update({_bin: report_fp})
        except subprocess.CalledProcessError as e:
            raise Exception(
                "An error was encountered while running Kraken 2, "
                f"(return code {e.returncode}), please inspect "
                "stdout and stderr to learn more."
            )
        results_df = _process_kraken2_reports(kraken2_reports)

        # TODO: make the level configurable?
        (table, taxonomy) = metaphlan_taxon(
            stratified_table=results_df, level=7
        )
        table.fillna(0, inplace=True)

    return table, taxonomy


def classify_kraken(
        seqs: MultiMAGSequencesDirFmt, db: str, threads: int = 1,
        confidence: float = 0.0, minimum_base_quality: int = 0,
        memory_mapping: bool = False, minimum_hit_groups: int = 2,
        quick: bool = False
) -> (pd.DataFrame, pd.DataFrame):
    kwargs = {
        k: v for k, v in locals().items()
        if k not in ["seqs"]
    }
    common_args = _process_common_input_params(
        processing_func=_process_kraken2_arg, params=kwargs
    )
    manifest: pd.DataFrame = seqs.manifest.view(pd.DataFrame)

    return _classify_kraken(manifest, common_args)
