# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json
import os
import tempfile
from copy import deepcopy

import pandas as pd
import pkg_resources
import q2templates

from q2_moshpit._utils import _process_common_input_params, run_command
from q2_moshpit.checkm.utils import (
    _process_checkm_arg, SampleStatistics, _extract_checkm_stats
)

from q2_types_genomics.per_sample_data._format import \
    (MultiMAGSequencesDirFmt)

TEMPLATES = pkg_resources.resource_filename('q2_moshpit', 'assets')


def _evaluate_bins(
        results_dir: str, bins: MultiMAGSequencesDirFmt,
        db_path: str, common_args: list
) -> dict:
    base_cmd = ['checkm', 'lineage_wf', *common_args]
    reports = {}

    manifest: pd.DataFrame = bins.manifest.view(pd.DataFrame)
    manifest['sample_dir'] = manifest.filename.apply(
        lambda x: os.path.dirname(x)
    )
    sample_dirs = manifest['sample_dir'].unique()
    for sample_dir in sample_dirs:
        sample = os.path.split(sample_dir)[-1]
        sample_results = os.path.join(results_dir, sample)
        report_path = os.path.join(sample_results, 'report.txt')
        reports[sample] = report_path

        cmd = deepcopy(base_cmd)
        cmd.extend(
            ['-x', 'fasta', '-f', report_path, sample_dir, sample_results]
        )
        run_command(cmd, env={**os.environ, 'CHECKM_DATA_PATH': db_path})
    return reports


def _parse_checkm_reports(reports: dict) -> pd.DataFrame:
    dfs = [
        _parse_single_checkm_report(_id, fp) for _id, fp in reports.items()
    ]
    results_df = pd.concat(dfs)
    results_df.reset_index(drop=True, inplace=True)
    return results_df


def _parse_single_checkm_report(
        sample_id: str, report_fp: str
) -> pd.DataFrame:
    sample_stats = SampleStatistics(sample_id=sample_id, bins=[])

    # read the raw CheckM report
    with open(report_fp, 'r') as fh:
        lines = fh.readlines()
        lines = [line.strip() for line in lines
                 if not line.startswith('---')][1:]

    # parse the report
    for line in lines:
        bin_stats = _extract_checkm_stats(line)
        sample_stats.bins.append(bin_stats)

    # convert report to DataFrame
    df = sample_stats.to_df()

    return df


def evaluate_bins(
    output_dir: str, bins: MultiMAGSequencesDirFmt,
    db_path: str, reduced_tree: bool = None, unique: int = None,
    multi: int = None, force_domain: bool = None, no_refinement: bool = None,
    individual_markers: bool = None, skip_adj_correction: bool = None,
    skip_pseudogene_correction: bool = None, aai_strain: float = None,
    ignore_thresholds: bool = None, e_value: float = None,
    length: float = None, threads: int = None, pplacer_threads: int = None
):

    kwargs = {k: v for k, v in locals().items()
              if k not in ['output_dir', 'bins', 'db_path']}
    common_args = _process_common_input_params(
        processing_func=_process_checkm_arg, params=kwargs
    )

    # TODO: check that CheckM's database is available (or fetch?)

    with tempfile.TemporaryDirectory() as tmp:
        results_dir = os.path.join(tmp, 'results')

        # run quast
        reports = _evaluate_bins(results_dir, bins, db_path, common_args)
        checkm_results = _parse_checkm_reports(reports)
        print(checkm_results)

        # fix/remove some URLs
        # _fix_html_reports(results_dir)
        #
        # copy_tree(os.path.join(TEMPLATES, 'quast'), output_dir)
        # copy_tree(results_dir, os.path.join(output_dir, 'quast_data'))
        #

        context = {
            'tabs': [
                {
                    'title': 'QC report',
                    'url': 'index.html'
                },
                # {
                #     'title': 'Contig browser',
                #     'url': 'q2_icarus.html'
                #
                # },
            ],
            'samples': json.dumps(reports.keys())
        }

        index = os.path.join(TEMPLATES, 'checkm', 'index.html')
        # icarus = os.path.join(TEMPLATES, 'quast', 'q2_icarus.html')

        templates = [index]
        q2templates.render(templates, output_dir, context=context)

    # TODO: remove later
    return common_args

    # return _bin_contigs_metabat(
    #     contigs=contigs, alignment_maps=alignment_maps,
    #     common_args=common_args
    # )
