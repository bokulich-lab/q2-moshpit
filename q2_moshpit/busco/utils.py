# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json
import os
import pandas as pd
from typing import List

arguments_with_hyphens = {
    "auto_lineage": "auto-lineage",
    "auto_lineage_euk": "auto-lineage-euk",
    "auto_lineage_prok": "auto-lineage-prok",
    "list_datasets": "list-datasets",
    "update_data": "update-data",
}

MARKER_COLS = ["single", "duplicated", "fragmented", "missing", "complete"]


def _parse_busco_params(arg_key, arg_val) -> List[str]:
    """Creates a list with argument and its value to be consumed by BUSCO.
    Argument names will be converted to command line parameters by
    appending a '--' prefix and in some cases replacing "_" for "-"
    (only for e.g. `arguments_with_hyphens`)

    Args:
        arg_key (str): Argument name.
        arg_val: Argument value.
    Returns:
        [converted_arg, arg_value]: List containing a prepared command line
            parameter and, optionally, its value.
    """

    # If the key is in arguments_with_hyphens, modify key
    if arg_key in arguments_with_hyphens.keys():
        arg_key = arguments_with_hyphens[arg_key]

    if isinstance(arg_val, bool):
        return [f"--{arg_key}"]
    else:
        return [f"--{arg_key}", str(arg_val)]


def _partition_dataframe(df, max_rows):
    groups = [group for _, group in df.groupby('sample_id')]
    partitions = []
    temp = []
    total_rows = 0

    for group in groups:
        if total_rows + len(group) > max_rows:
            partitions.append(pd.concat(temp))
            temp = [group]
            total_rows = len(group)
        else:
            temp.append(group)
            total_rows += len(group)

    if temp:
        partitions.append(pd.concat(temp))

    return partitions


def _collect_summaries(run_summaries_fp_map: dict) -> pd.DataFrame:
    """
    Reads-in the sample-wise summaries and concatenates them in one
    pd.DataFrame, which is saved to file.

    Args:
        run_summaries_fp_map (dict): dict where key is sample id
            and value is path "tmp/sample_id/batch_summary.txt"

    Returns:
        all_summaries (pd.DataFrame): DataFrame composed of the individual
            run summaries.
    """

    all_summaries = []
    for sample_id, path_to_summary in run_summaries_fp_map.items():
        df = pd.read_csv(filepath_or_buffer=path_to_summary, sep="\t")
        df["sample_id"] = sample_id
        all_summaries.append(df)

    return pd.concat(all_summaries, ignore_index=True)


def _dump_spec(context, output_dir, vega_spec_fn):
    # vega_out_fp = os.path.join(output_dir, vega_spec_fn)
    # with open(vega_out_fp, 'w') as json_file:
    vega_json = json.dumps(context)
        # json_file.write(vega_json)
    return vega_json


def _get_feature_table(busco_results: pd.DataFrame):
    df = busco_results.reset_index(inplace=False, drop=False)
    for col in MARKER_COLS:
        df[col] = df[col] * df["n_markers"] / 100
    new_cols = {
        "mag_id": "MAG", "sample_id": "Sample", "dataset": "Dataset",
        "single": "Single", "duplicated": "Duplicated",
        "fragmented": "Fragmented", "missing": "Missing",
        "complete": "Complete", "n_markers": "Total markers",
        "contigs_n50": "N50 contigs", "percent_gaps": "Percent gaps",
        "scaffolds": "Contigs"
    }
    df = df[list(new_cols.keys())].rename(columns=new_cols, inplace=False)
    return df.to_json(orient='split')


def _parse_df_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds several columns required for generation of downloadable 
    BUSCO plots.

    Args:
        df (pd.DataFrame): Unformatted DataFrame

    Returns:
        df (pd.DataFrame): Formatted DataFrame
    """
    df = df.reset_index(drop=False, inplace=False)
    df = df.rename(columns={"id": "mag_id"}, inplace=False)

    # fix data types
    df["percent_gaps"] = df["percent_gaps"].str.split(
        '%', expand=True
    )[0].map(float)
    for col in MARKER_COLS:
        df[col] = df[col].map(float)
    df["n_markers"] = df["n_markers"].map(int)

    return df


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


def _cleanup_bootstrap(output_dir):
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


def _calculate_summary_stats(df: pd.DataFrame) -> json:
    stats = pd.DataFrame({
        "min": df[MARKER_COLS].min(),
        "median": df[MARKER_COLS].median(),
        "mean": df[MARKER_COLS].mean(),
        "max": df[MARKER_COLS].max(),
        "count": df[MARKER_COLS].count()
    })
    return stats.T.to_json(orient='table')
