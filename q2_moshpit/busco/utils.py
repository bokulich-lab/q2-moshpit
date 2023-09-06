import json
import os
import pandas as pd
import altair as alt
import seaborn as sns
import matplotlib.pyplot as plt
from zipfile import ZipFile
from .._utils import run_command
from copy import deepcopy
from typing import List, Dict
from q2_types_genomics.per_sample_data._format import MultiMAGSequencesDirFmt

arguments_with_hyphens = {
    "auto_lineage": "auto-lineage",
    "auto_lineage_euk": "auto-lineage-euk",
    "auto_lineage_prok": "auto-lineage-prok",
    "list_datasets": "list-datasets",
    "update_data": "update-data",
}


def _parse_busco_params(arg_key, arg_val) -> List[str]:
    """Creates a list with argument and its value to be consumed by MetaBAT 2.

    Argument names will be converted to command line parameters by
    appending a '--' prefix and concatenating words separated by a '_',
    e.g.: 'some_parameter_x' -> '--someParameterX'.

    Args:
        arg_key (str): Argument name.
        arg_val: Argument value.

    Returns:
        [converted_arg, arg_value]: List containing a prepared command line
            parameter and, optionally, its value.
    """

    # If the key is one of
    if arg_key in arguments_with_hyphens.keys():
        arg_key = arguments_with_hyphens[arg_key]

    if isinstance(arg_val, bool) and arg_val:
        return [f"--{arg_key}"]
    else:
        return [f"--{arg_key}", str(arg_val)]


def _draw_busco_plots_for_render(
    df: pd.core.frame.DataFrame,
    width: int = None,
    height: int = None,
    labelFontSize: int = None,
    titleFontSize: int = None,
) -> str:
    """
    Draws a hroizontal normalized bar plot for every sample for which BUSCO was
    run. Each barplot shows the BUSCO results for each of the MAGs in the
    sample. The plots for all samples are drwan in one composite plot which
    is then returned as a dictionary for rendering (but casted to a string).

    Args:
        df (pd.core.frame.DataFrame): tabular batch summary for all samples
        width (int): width of the plot
        height (int): height of the plot
        labelFontSize (int): size of the labels in plot
        titleFontSize (int): size of titles in plot

    Output:
        Output plot in dictionary from casted to a string.
    """

    # Clean column names
    df.columns = df.columns.str.replace(" ", "_")
    df.columns = df.columns.str.replace("[^a-zA-Z0-9_]", "", regex=True)
    df.columns = df.columns.str.lower()

    # Rename column "input_file"
    df["mag_id"] = df["input_file"].str.split(".", expand=True)[0]

    # Get number of samples
    n_samples = len(df["mag_id"].unique())

    # Pivot long
    df2 = pd.melt(
        df,
        id_vars=["sample_id", "mag_id", "dataset", "n_markers"],
        value_vars=["single", "duplicated", "fragmented", "missing"],
        value_name="BUSCO_percentage",
        var_name="category",
    )

    # Specify order
    mapping = {"single": 1, "duplicated": 2, "fragmented": 3, "missing": 4}
    df2["order"] = df2["category"].map(mapping)

    # Estimate fraction of sequences in each BUSCO category
    df2["fracc_markers"] = (
        "~"
        + round(
            df2["BUSCO_percentage"] * df2["n_markers"] / 100
        ).map(int).map(str)
        + "/124"
    )

    # Plot
    domain = ["single", "duplicated", "fragmented", "missing"]
    range_ = ["#1E90FF", "#87CEFA", "#FFA500", "#FF7F50"]

    output_plot = (
        alt.Chart(df2)
        .mark_bar()
        .encode(
            x=alt.X(
                "sum(BUSCO_percentage)",
                stack="normalize",
                title="BUSCO fracc."
            ),
            y=alt.Y("mag_id", axis=alt.Axis(title="MAG ID")),
            color=alt.Color(
                "category",
                scale=alt.Scale(domain=domain, range=range_),
                legend=alt.Legend(title="BUSCO Category"),
            ),
            order=alt.Order("order", sort="ascending"),
            tooltip=[
                alt.Tooltip("sample_id", title="Sample ID"),
                alt.Tooltip("mag_id", title="MAG ID"),
                alt.Tooltip("dataset", title="Lineage dataset"),
                alt.Tooltip(
                    "fracc_markers",
                    title="Aprox. number of markers in this category"
                ),
                alt.Tooltip("BUSCO_percentage", title="Percentage [%]"),
            ],
            opacity=alt.value(0.85),
        )
        .properties(width=width, height=height * n_samples)
        .facet(row=alt.Row("sample_id", title="Sample ID"))
        .resolve_scale(y="independent")
    )

    # Resize text and plot
    output_plot = output_plot.configure_axis(
        labelFontSize=labelFontSize, titleFontSize=titleFontSize
    )
    output_plot = output_plot.configure_legend(
        labelFontSize=labelFontSize, titleFontSize=titleFontSize
    )
    output_plot = output_plot.configure_header(
        labelFontSize=labelFontSize, titleFontSize=titleFontSize
    )

    # Return
    return json.dumps(output_plot.to_dict())


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
            `tmp/busco_output/<sampl_id>/batch_summary.txt`.
    """

    # Define base command
    base_cmd = ["busco", *params]

    # Creates pandas df "manifest" from bins
    manifest: pd.DataFrame = mags.manifest.view(pd.DataFrame)

    # Make a new column in manifest with the directories of files
    # listed in column "filename"
    manifest["sample_dir"] = manifest.filename.apply(
        lambda x: os.path.dirname(x)
    )

    # numpy.ndarray with unique dirs
    sample_dirs = manifest["sample_dir"].unique()

    # Initialize dictionary with paths to run summeries
    path_to_run_summeries = {}

    # For every unique sample dir run busco
    for sample_dir in sample_dirs:
        # Get sample ide from tip dirname
        sample = os.path.split(sample_dir)[-1]

        # Deep copy base comand extend it with the sample specific
        # info and run it
        cmd = deepcopy(base_cmd)
        cmd.extend([
            "--in",
            sample_dir,
            "--out_path",
            output_dir,
            "-o",
            sample
        ])
        run_command(cmd, env={**os.environ})

        # Check for output
        path_to_run_summerie = os.path.join(
            output_dir, sample, "batch_summary.txt"
        )
        if os.path.isfile(path_to_run_summerie):
            path_to_run_summeries[sample] = path_to_run_summerie
        else:
            raise FileNotFoundError(
                f"BUSCO batch summary file {path_to_run_summerie} not found."
            )

    # Return a dict where key is sample name and value is path
    # "tmp/sample/batch_summary.txt"
    return path_to_run_summeries


def _draw_busco_plots(
        path_to_run_summeries: dict, output_dir: str
        ) -> Dict[str, str]:
    """
    Generates plots for all `batch_summary.txt` (one for every sample)
    and saves them to `output_dir`.

    Args:
        output_dir (str): Location where the final results should be stored.
        mags (MultiMAGSequencesDirFmt): The mags to be analyzed.
        params (List[str]): List of parsed arguments to pass to BUSCO

    Returns:
        dict: Dictionary where keys are sample IDs and values are the paths
            to the generated plots, e.g.
            `tmp/plots/<sampl_id>/plot_batch_summary.svg`.
    """

    # Initialize output dictionary
    paths_to_plots = {}

    # For every sample make a plot
    for sample_id, path_to_summary in path_to_run_summeries.items():
        # Read in text file as dataframe
        df = pd.read_csv(filepath_or_buffer=path_to_summary, sep="\t")

        # Compute cumulative percentages
        df["single"] = df["Single"]
        df["duplicated"] = df["Single"] + df["Duplicated"]
        df["fragmented"] = df["duplicated"] + df["Fragmented"]
        df["missing"] = df["fragmented"] + df["Missing"]

        # Get sample id without extension (.fasta)
        df["input_file"] = df["Input_file"].str.split(".", expand=True)[0]

        # Set the style
        sns.set(style="whitegrid")
        sns.set_palette("colorblind")

        # Create a horizontal stacked barplot
        plt.figure(figsize=(10, 6))
        sns.barplot(
            data=df, y="input_file", x="missing", color="r", label="Missing"
        )
        sns.barplot(
            data=df,
            y="input_file",
            x="fragmented",
            color="tab:orange",
            label="Fragmented",
        )
        sns.barplot(
            data=df,
            y="input_file",
            x="duplicated",
            color="tab:cyan",
            label="Duplicated",
        )
        sns.barplot(
            data=df,
            y="input_file",
            x="single",
            color="tab:blue",
            label="Single",
        )

        # Customize the plot
        plt.xlabel("%BUSCO")
        plt.ylabel("MAG ID's")
        plt.legend(loc="lower right")
        plt.title(f"Sample ID: {sample_id}")

        # Save figure to file
        output_name = os.path.join(
            output_dir, sample_id, "plot_batch_summary.svg"
        )
        os.makedirs(os.path.dirname(output_name), exist_ok=True)
        plt.savefig(output_name, format="svg")

        # Save path to dictionary
        paths_to_plots[sample_id] = output_name

    # Return paths to all generated plots
    return paths_to_plots


def _zip_busco_plots(paths_to_plots: dict, zip_path: str) -> None:
    """
    Creates a single zip archive containing all plots produced by BUSCO,
    one for each sample.

    Args:
        paths_to_plots: Dictionary mapping sample to plot path.
        zip_path (str): The path to the zip archive.
    """

    # Get shortest comon path between files
    common_path = os.path.commonpath(paths_to_plots.values())

    # Write to zipfile
    with ZipFile(zip_path, "w") as zf:
        for _, path_to_plot in paths_to_plots.items():
            arcname = os.path.relpath(path_to_plot, common_path)
            zf.write(path_to_plot, arcname=arcname)