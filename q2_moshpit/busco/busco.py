# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import tempfile

import pandas as pd

from q2_moshpit.busco.utils import (
    _parse_busco_params, _render_html, _run_busco,
    _collect_summaries, _draw_busco_plots,
    _zip_busco_plots,
)
from q2_moshpit._utils import _process_common_input_params
from q2_types.per_sample_sequences._format import MultiMAGSequencesDirFmt


def _rename_columns(df):
    cols = {
        "Input_file": "input_file", "Dataset": "dataset",
        "Complete": "complete", "Single": "single",
        "Duplicated": "duplicated", "Fragmented": "fragmented",
        "Missing": "missing", "n_markers": "n_markers",
        "Scaffold N50": "scaffold_n50", "Contigs N50": "contigs_n50",
        "Percent gaps": "percent_gaps", "Number of scaffolds": "scaffolds",
        "sample_id": "sample_id"
    }

    cols_reshuffled = [
        "mag_id", "sample_id", "input_file", "dataset", "complete",
        "single", "duplicated", "fragmented", "missing", "n_markers",
        "scaffold_n50", "contigs_n50", "percent_gaps", "scaffolds",
    ]

    df = df.rename(columns=cols, inplace=False)
    df["mag_id"] = df["input_file"].str.split(".", expand=True)[0]
    return df[cols_reshuffled]


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

    return all_summaries


def _evaluate_busco(
    bins: MultiMAGSequencesDirFmt,
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
        k: v for k, v in locals().items() if k not in ["bins",]
    }

    # Filter out all kwargs that are None, False or 0.0
    common_args = _process_common_input_params(
        processing_func=_parse_busco_params, params=kwargs
    )

    return _busco_helper(bins, common_args)


def _visualize_busco(output_dir: str, busco_results: pd.DataFrame) -> None:
    busco_results.to_csv(
        os.path.join(output_dir, "all_batch_summaries.csv"),
        index=False
    )
    with tempfile.TemporaryDirectory() as tmp:
        # Draw BUSCO plots for all samples
        # Result NOT included in final output
        paths_to_plots = _draw_busco_plots(
            data=busco_results,
            plots_dir=os.path.join(tmp, "plots")
        )

        _zip_busco_plots(
            paths_to_plots=paths_to_plots,
            zip_path=os.path.join(output_dir, "busco_plots.zip")
        )

        _render_html(output_dir, busco_results)


def evaluate_busco(
    ctx,
    bins,
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
    print(f"Partitioned MAGs: {partitioned_mags}")
    for mag in partitioned_mags.values():
        print(f"Processing MAG: {mag}")
        (busco_result, ) = _evaluate_busco(mag, **kwargs)
        print(f"Got BUSCO result: {busco_result}")
        results.append(busco_result)

    collated_results = collate_busco_results(results)
    visualization = _visualize_busco(collated_results)

    return collated_results, visualization
