# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import os
import json
import tempfile
import q2templates
import pandas as pd
from shutil import copytree
from .utils import (
    _parse_busco_params,
    _draw_busco_plots_for_render,
    _run_busco,
    _draw_busco_plots,
    _zip_busco_plots,
)
from .._utils import _process_common_input_params
from typing import List
from q2_types_genomics.per_sample_data._format import MultiMAGSequencesDirFmt


def busco(
    output_dir: str,
    bins: MultiMAGSequencesDirFmt,
    mode: str = None,
    lineage_dataset: str = None,
    augustus: bool = None,
    augustus_parameters: str = None,
    augustus_species: str = None,
    auto_lineage: bool = None,
    auto_lineage_euk: bool = None,
    auto_lineage_prok: bool = None,
    cpu: int = None,
    config: str = None,
    contig_break: int = None,
    datasets_version: str = None,
    download: List[str] = None,
    download_base_url: str = None,
    download_path: str = None,
    evalue: float = None,
    force: bool = None,
    limit: int = None,
    help: bool = None,
    list_datasets: bool = None,
    long: bool = None,
    metaeuk_parameters: str = None,
    metaeuk_rerun_parameters: str = None,
    miniprot: bool = None,
    offline: bool = None,
    quiet: bool = None,
    restart: bool = None,
    scaffold_composition: bool = None,
    tar: bool = None,
    update_data: bool = None,
    version: bool = None,
) -> None:
    """
    qiime2 vizualization for the BUSCO assessment tool
    <https://busco.ezlab.org/>.

    Args:
        see all possible inputs by running `qiime moshpit plot_busco`

    Output:
        plots.zip: zip file containing all of the busco plots
        busco_output: all busco outputfiles
        qiime_html: html for rendering the output plots
    """

    # Create dictionary with local varaibles
    # (kwards passed to the function or their defaults) excluding
    # "output_dir" and "bins"
    kwargs = {
        k: v for k, v in locals().items() if k not in ["output_dir", "bins"]
    }

    # Filter out all kwards that are None, False or 0.0
    common_args = _process_common_input_params(
        processing_func=_parse_busco_params, params=kwargs
    )

    # Creates output directory with path 'tmp'
    with tempfile.TemporaryDirectory() as tmp:
        # Run busco for every sample. Returns dictionary to report files.
        results_dir = os.path.join(tmp, "busco_output")
        path_to_run_summeries = _run_busco(
            output_dir=results_dir, mags=bins, params=common_args
        )

        # Collect result for each sample and save to file.
        all_summeries_list = []
        for sample_id, path_to_summary in path_to_run_summeries.items():
            df = pd.read_csv(filepath_or_buffer=path_to_summary, sep="\t")
            df["sample_id"] = sample_id
            all_summeries_list.append(df)

        # Concatenate
        all_summeries_df = pd.concat(all_summeries_list)

        # Save to file
        all_summeries_path = os.path.join(
            output_dir, "all_batch_summeries.csv"
            )
        all_summeries_df.to_csv(all_summeries_path, index=False)

        # Draw BUSCO plots for all samples
        plots_dir = os.path.join(tmp, "plots")
        paths_to_plots = _draw_busco_plots(
            path_to_run_summeries=path_to_run_summeries, output_dir=plots_dir
        )

        # Zip graphs for user download
        zip_name = os.path.join(output_dir, "busco_plots.zip")
        _zip_busco_plots(paths_to_plots=paths_to_plots, zip_path=zip_name)

        # Render qiime html report
        # Prepare context for jinja2 template
        context = {
            "tabs": [{"title": "BUSCO Results", "url": "index.html"}],
            "samples": json.dumps(list(path_to_run_summeries.keys())),
            "vega_plots_overview": _draw_busco_plots_for_render(
                all_summeries_df,
                width=600,
                height=18,
                titleFontSize=20,
                labelFontSize=17,
            ),
        }

        # Copy BUSCO results from tmp dir to output_dir
        moshpit_path = os.path.dirname(  # Path to parent dir, q2_moshpit
            os.path.dirname(  # Path to parent dir, e.g. busco
                os.path.abspath(__file__)  # Path to this file
            )
        )
        TEMPLATES = os.path.join(moshpit_path, "assets")
        index = os.path.join(TEMPLATES, "busco", "index.html")
        copytree(
            src=os.path.join(TEMPLATES, "busco"),
            dst=output_dir,
            dirs_exist_ok=True
        )
        copytree(
            src=plots_dir,
            dst=os.path.join(output_dir, "plots"),
            dirs_exist_ok=True
        )

        # Render
        templates = [index]
        q2templates.render(templates, output_dir, context=context)

        # Remove unwanted files
        # until Bootstrap 3 is replaced with v5, remove the v3 scripts as
        # the HTML files are adjusted to work with v5
        os.remove(
            os.path.join(
                output_dir, "q2templateassets", "css", "bootstrap.min.css"
                )
        )
        os.remove(
            os.path.join(
                output_dir, "q2templateassets", "js", "bootstrap.min.js"
                )
        )
