# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import re
import subprocess
import tempfile
import time
from copy import deepcopy
from functools import reduce

import numpy as np
import pandas as pd
import requests

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


def _fetch_tax_id(taxon: str) -> int:
    taxon = taxon.replace('[', '').replace(']', '')
    params = {
        'db': 'taxonomy',
        'term': f'{taxon}[Scientific Name]',
        'retmode': 'json',
        'retmax': 100,
        'retstart': 0
    }

    response = requests.get(
        url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi',
        params=params,
    )
    result = response.json()['esearchresult']['idlist']
    time.sleep(0.4)
    if len(result) == 1:
        return result[0]
    else:
        print('More than one result was found for %s', taxon)
        return min([int(x) for x in result])


def _process_kraken2_reports(reports: dict) -> pd.DataFrame:
    # get abundances per sample
    abundances = []
    for _sample, bins in reports.items():
        sample_abundances = []
        for _bin, reports in bins.items():
            bin_abundances = _parse_kraken2_report(
                reports['report'], _sample, _bin
            )
            sample_abundances.append(bin_abundances)
        # merge bins into one per sample
        bins_merged = reduce(
            lambda left, right: pd.merge(
                left, right, left_index=True, right_index=True, how='outer'
            ), sample_abundances
        )
        bins_merged = bins_merged.fillna(0).sum(axis=1)
        bins_merged.name = _sample
        abundances.append(bins_merged)

    # combine all samples
    df_merged = reduce(
        lambda left, right: pd.merge(
            left, right, left_index=True, right_index=True, how='outer'
        ), abundances
    )

    # find all levels from all samples
    all_levels = []
    for taxonomy in df_merged.index.tolist():
        for level in taxonomy.split('|'):
            all_levels.append(level) if level not in all_levels else False

    # fetch IDs from NCBI
    # levels_with_ids = {
    #     level: _fetch_tax_id(level.split('__')[-1]) for level in all_levels
    # }
    levels_with_ids = {
        level: _id for _id, level in enumerate(sorted(all_levels))
    }

    # assign NCBI IDs
    all_taxonomies = []
    for taxonomy in df_merged.index:
        _id = list(map(lambda x: str(levels_with_ids[x]), taxonomy.split('|')))
        all_taxonomies.append('|'.join(_id))

    df_merged['NCBI_tax_id'] = all_taxonomies
    df_merged.index.name = 'feature-id'

    return df_merged


def _parse_kraken2_output(
        report_fp: str, sample_name: str, bin_name: str
) -> pd.DataFrame:
    df = pd.read_csv(report_fp, sep='\t', header=None)
    df.columns = ['status', 'id', 'tax_id', 'length', 'lca_mapping']
    df['sample'] = sample_name
    df['split_tax_id'] = df['tax_id'].apply(
        lambda x: re.findall(r'(.{1,})\(taxid (\d{1,})\)', x)[0]
    )
    df[['taxon', 'tax_id']] = df['split_tax_id'].tolist()
    df.drop('split_tax_id', axis=1, inplace=True)
    df.set_index('id', drop=True, inplace=True)
    return df


def _process_kraken2_outputs(reports: dict) -> pd.DataFrame:
    contig_info_all = []
    for _sample, bins in reports.items():
        for _bin, _reports in bins.items():
            contig_info = _parse_kraken2_output(
                _reports['output'], _sample, _bin
            )
            contig_info_all.append(contig_info)
    return pd.concat(contig_info_all, axis=0)


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
                    '--use-mpa-style', '--report', report_fp, '--use-names',
                    '--output', output_fp, fn
                ])
                run_command(cmd=cmd, verbose=True)
                kraken2_reports[_sample].update(
                    {_bin: {'report': report_fp, 'output': output_fp}}
                )
        except subprocess.CalledProcessError as e:
            raise Exception(
                "An error was encountered while running Kraken 2, "
                f"(return code {e.returncode}), please inspect "
                "stdout and stderr to learn more."
            )
        results_df = _process_kraken2_reports(kraken2_reports)
        #contigs_to_taxons = _process_kraken2_outputs(kraken2_reports)

        #grouped = contigs_to_taxons.groupby(['sample', 'taxon']).sum()

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
