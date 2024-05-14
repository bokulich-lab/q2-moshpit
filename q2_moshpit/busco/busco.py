# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json

import os
import tempfile
from copy import deepcopy
from shutil import copytree
from typing import List, Dict

import pandas as pd
import q2templates

from q2_moshpit.busco.plots_detailed import _draw_detailed_plots
from q2_moshpit.busco.plots_summary import _draw_marker_summary_histograms, \
    _draw_selectable_summary_histograms

from q2_moshpit.busco.utils import (
    _parse_busco_params, _collect_summaries, _rename_columns,
    _parse_df_columns, _partition_dataframe, _calculate_summary_stats,
    _get_feature_table, _cleanup_bootstrap, _get_mag_lengths,
    _validate_lineage_dataset_input
)
from q2_moshpit._utils import _process_common_input_params, run_command
from q2_types.per_sample_sequences._format import MultiMAGSequencesDirFmt
from q2_moshpit.busco.types import BuscoDatabaseDirFmt


def _run_busco(
    output_dir: str, mags: MultiMAGSequencesDirFmt, params: List[str]
) -> Dict[str, str]:
    """Evaluates bins for all samples using BUSCO.

    Args:
        output_dir (str): Location where the final results should be stored.
        mags (MultiMAGSequencesDirFmt): The mags to be analyzed.
        params (List[str]): List of parsed arguments to pass to BUSCO.

    Returns:
        dict: Dictionary where keys are sample IDs and values are the paths
            to the `batch_summary.txt` generated by BUSCO, e.g.
            `tmp/busco_output/<sample_id>/batch_summary.txt`.
    """
    base_cmd = ["busco", *params]

    manifest: pd.DataFrame = mags.manifest.view(pd.DataFrame)
    manifest["sample_dir"] = manifest.filename.apply(
        lambda x: os.path.dirname(x)
    )

    sample_dirs = manifest["sample_dir"].unique()

    path_to_run_summaries = {}

    # For every unique sample dir run busco
    for sample_dir in sample_dirs:
        sample = os.path.split(sample_dir)[-1]

        cmd = deepcopy(base_cmd)
        cmd.extend([
            "--in",
            sample_dir,
            "--out_path",
            output_dir,
            "-o",
            sample
        ])
        run_command(cmd)

        path_to_run_summary = os.path.join(
            output_dir, sample, "batch_summary.txt"
        )
        if os.path.isfile(path_to_run_summary):
            path_to_run_summaries[sample] = path_to_run_summary
        else:
            raise FileNotFoundError(
                f"BUSCO batch summary file {path_to_run_summary} not found."
            )

    return path_to_run_summaries


def _busco_helper(bins, common_args):
    with tempfile.TemporaryDirectory() as tmp:
        path_to_run_summaries = _run_busco(
            output_dir=os.path.join(tmp, "busco_output"),
            mags=bins,
            params=common_args,
        )

        all_summaries = _collect_summaries(
            run_summaries_fp_map=path_to_run_summaries,
        )
    all_summaries = _rename_columns(all_summaries)

    lengths = _get_mag_lengths(bins)
    all_summaries = all_summaries.join(lengths, on="mag_id")

    return all_summaries


def _evaluate_busco(
    bins: MultiMAGSequencesDirFmt,
    busco_db: BuscoDatabaseDirFmt = None,
    mode: str = "genome",
    lineage_dataset: str = None,
    augustus: bool = False,
    augustus_parameters: str = None,
    augustus_species: str = None,
    auto_lineage: bool = False,
    auto_lineage_euk: bool = False,
    auto_lineage_prok: bool = False,
    cpu: int = 1,
    config: str = None,
    contig_break: int = 10,
    evalue: float = 1e-03,
    force: bool = False,
    limit: int = 3,
    long: bool = False,
    metaeuk_parameters: str = None,
    metaeuk_rerun_parameters: str = None,
    miniprot: bool = False,
    scaffold_composition: bool = False,
) -> pd.DataFrame:
    kwargs = {
        k: v for k, v in locals().items() if k not in ["bins", "busco_db"]
    }

    # Add busco_db to kwargs
    if busco_db is not None:
        kwargs["offline"] = True
        kwargs["download_path"] = f"{str(busco_db)}/busco_downloads"

    if lineage_dataset is not None:
        _validate_lineage_dataset_input(
            lineage_dataset, auto_lineage, auto_lineage_euk, auto_lineage_prok,
            busco_db, kwargs  # this may be modifies inside this function
        )

    # Filter out all kwargs that are None, False or 0.0
    common_args = _process_common_input_params(
        processing_func=_parse_busco_params, params=kwargs
    )

    return _busco_helper(bins, common_args)


def _visualize_busco(output_dir: str, busco_results: pd.DataFrame) -> None:
    busco_results.to_csv(
        os.path.join(output_dir, "busco_results.csv"),
        index=False
    )
    busco_results = _parse_df_columns(busco_results)
    dfs = _partition_dataframe(busco_results, max_rows=100)

    context = {}
    counter_left = 1
    for i, df in enumerate(dfs):
        sample_count = df['sample_id'].nunique()
        counter_right = counter_left + sample_count - 1
        sample_counter = {"from": counter_left, "to": counter_right}
        counter_left += sample_count
        subcontext = _draw_detailed_plots(
            df,
            width=600,
            height=30,
            title_font_size=20,
            label_font_size=17,
            spacing=20
        )
        context.update(
            {f"sample{i}": {
                "subcontext": subcontext,
                "sample_counter": sample_counter,
                "sample_ids": df['sample_id'].unique().tolist(),
            }}
        )

    marker_summary_spec = _draw_marker_summary_histograms(busco_results)
    selectable_summary_spec = _draw_selectable_summary_histograms(
        busco_results
    )

    vega_json = json.dumps(context)
    vega_json_summary = json.dumps(marker_summary_spec)
    vega_json_summary_selectable = json.dumps(selectable_summary_spec)

    # Copy BUSCO results from tmp dir to output_dir
    TEMPLATES = os.path.join(
        os.path.dirname(os.path.dirname(__file__)), "assets"
    )
    index = os.path.join(TEMPLATES, "busco", "index.html")
    details = os.path.join(TEMPLATES, "busco", "detailed_view.html")
    table = os.path.join(TEMPLATES, "busco", "table.html")
    copytree(
        src=os.path.join(TEMPLATES, "busco"),
        dst=output_dir,
        dirs_exist_ok=True
    )

    table_json = _get_feature_table(busco_results)
    stats_json = _calculate_summary_stats(busco_results)

    # Render
    tabbed_context = {
        "tabs": [
            {"title": "QC overview", "url": "index.html"},
            {"title": "Sample details", "url": "detailed_view.html"},
            {"title": "Feature details", "url": "table.html"}
        ],
        "vega_json": vega_json,
        "vega_summary_json": vega_json_summary,
        "vega_summary_selectable_json": vega_json_summary_selectable,
        "table": table_json,
        "summary_stats_json": stats_json,
        "page_size": 100
    }
    templates = [index, details, table]
    q2templates.render(templates, output_dir, context=tabbed_context)

    # Final cleanup, needed until we fully migrate to Bootstrap 5
    _cleanup_bootstrap(output_dir)


def evaluate_busco(
    ctx,
    bins,
    busco_db=None,
    mode="genome",
    lineage_dataset=None,
    augustus=False,
    augustus_parameters=None,
    augustus_species=None,
    auto_lineage=False,
    auto_lineage_euk=False,
    auto_lineage_prok=False,
    cpu=1,
    config=None,
    contig_break=10,
    evalue=1e-03,
    force=False,
    limit=3,
    long=False,
    metaeuk_parameters=None,
    metaeuk_rerun_parameters=None,
    miniprot=False,
    scaffold_composition=False,
    num_partitions=None
):

    kwargs = {
        k: v for k, v in locals().items()
        if k not in ["bins", "ctx", "num_partitions"]
    }

    _evaluate_busco = ctx.get_action("moshpit", "_evaluate_busco")
    partition_mags = ctx.get_action("moshpit", "partition_sample_data_mags")
    collate_busco_results = ctx.get_action("moshpit", "collate_busco_results")
    _visualize_busco = ctx.get_action("moshpit", "_visualize_busco")

    (partitioned_mags, ) = partition_mags(bins, num_partitions)
    results = []
    for mag in partitioned_mags.values():
        (busco_result, ) = _evaluate_busco(mag, **kwargs)
        results.append(busco_result)

    collated_results, = collate_busco_results(results)
    visualization, = _visualize_busco(collated_results)

    return collated_results, visualization
