import os
import q2templates
from shutil import copytree
import pandas as pd
import altair as alt
import seaborn as sns
import matplotlib.pyplot as plt
from zipfile import ZipFile
from .._utils import run_command
from copy import deepcopy
from typing import List, Dict
from q2_types_genomics.per_sample_data._format import MultiMAGSequencesDirFmt


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

    if isinstance(arg_val, bool) and arg_val:
        return [f"--{arg_key.replace('_', '-')}"]
    else:
        return [f"--{arg_key.replace('_', '-')}", str(arg_val)]


def _draw_busco_plots_for_render(
    df: pd.DataFrame,
    width: int = None,
    height: int = None,
    labelFontSize: int = None,
    titleFontSize: int = None,
) -> str:
    """
    Draws a horizontal normalized bar plot for every sample for which BUSCO was
    run. Each barplot shows the BUSCO results for each of the MAGs in the
    sample. The plots for all samples are drawn in one composite plot which
    is then returned as a dictionary for rendering (but casted to a string).

    Args:
        df (pd.DataFrame): tabular batch summary for all samples
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

    # Cast into percent_gaps col to float
    df["percent_gaps"] = df["percent_gaps"].str.split(
        '%', expand=True
    )[0].map(float)

    # Get number of samples
    n_samples = len(df["mag_id"].unique())

    # Format data for plotting
    busco_plot_data = pd.melt(
        df,
        id_vars=["sample_id", "mag_id", "dataset", "n_markers"],
        value_vars=["single", "duplicated", "fragmented", "missing"],
        value_name="BUSCO_percentage",
        var_name="category",
    )

    secondary_plot_data = df[[
        "sample_id",
        "mag_id",
        'scaffold_n50',
        'contigs_n50',
        'percent_gaps',
        'number_of_scaffolds',
    ]]

    # Specify order
    mapping = {"single": 1, "duplicated": 2, "fragmented": 3, "missing": 4}
    busco_plot_data["order"] = busco_plot_data["category"].map(mapping)

    # Estimate fraction of sequences in each BUSCO category
    busco_plot_data["fracc_markers"] = (
        "~"
        + round(
            busco_plot_data["BUSCO_percentage"] *
            busco_plot_data["n_markers"] / 100
        ).map(int).map(str)
        + "/" + busco_plot_data["n_markers"].map(str)
    )

    # Plot
    domain = ["single", "duplicated", "fragmented", "missing"]
    range_ = ["#1E90FF", "#87CEFA", "#FFA500", "#FF7F50"]

    busco_plot = (
        alt.Chart(busco_plot_data)
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
                legend=alt.Legend(title="BUSCO Category", orient="top"),
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

    # Secondary plot
    # Drop down menu
    dropdown = alt.binding_select(
        options=[
            'scaffold_n50',
            'contigs_n50',
            'percent_gaps',
            'number_of_scaffolds',
        ],
        name="Assambly Statistics: "
    )

    xcol_param = alt.param(
        value='scaffold_n50',
        bind=dropdown
    )

    secondary_plot = alt.Chart(secondary_plot_data).mark_bar().encode(
        x=alt.X('x:Q').title('Assembly Statistic'),
        y=alt.Y('mag_id', axis=None),
        tooltip=[alt.Tooltip('x:Q', title="value")],
        opacity=alt.value(0.85)
    ).transform_calculate(
        x=f'datum[{xcol_param.name}]'
    ).add_params(
        xcol_param
    ).properties(
        width=width,
        height=height * n_samples
    ).facet(
        row=alt.Row(
            "sample_id",
            title=None,
            header=alt.Header(labelFontSize=0)
        )
    ).resolve_scale(
        y="independent"
    )

    # concatenate plots horizontally
    output_plot = alt.hconcat(
        busco_plot, secondary_plot, spacing=3
    ).configure_axis(
        labelFontSize=labelFontSize, titleFontSize=titleFontSize
    ).configure_legend(
        labelFontSize=labelFontSize, titleFontSize=titleFontSize
    ).configure_header(
        labelFontSize=labelFontSize, titleFontSize=titleFontSize
    )

    # Return
    return output_plot.to_json()


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

    # Initialize dictionary with paths to run summaries
    path_to_run_summaries = {}

    # For every unique sample dir run busco
    for sample_dir in sample_dirs:
        # Get sample ide from tip dirname
        sample = os.path.split(sample_dir)[-1]

        # Deep copy base command extend it with the sample specific
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
        run_command(cmd)

        # Check for output
        path_to_run_summary = os.path.join(
            output_dir, sample, "batch_summary.txt"
        )
        if os.path.isfile(path_to_run_summary):
            path_to_run_summaries[sample] = path_to_run_summary
        else:
            raise FileNotFoundError(
                f"BUSCO batch summary file {path_to_run_summary} not found."
            )

    # Return a dict where key is sample id and value is path
    # "tmp/sample_id/batch_summary.txt"
    return path_to_run_summaries


def _draw_busco_plots(
        path_to_run_summaries: dict, plots_dir: str
        ) -> Dict[str, str]:
    """
    Generates plots for all `batch_summary.txt` (one for every sample)
    and saves them to `plots_dir`.

    Args:
        plots_dir (str): Path where the results should be stored.
        dict: Dictionary where keys are sample IDs and values are the paths
            to the `batch_summary.txt` generated by BUSCO, e.g.
            `tmp/busco_output/<sample_id>/batch_summary.txt`.

    Returns:
        dict: Dictionary where keys are sample IDs and values are the paths
            to the generated plots, e.g.
            `tmp/plots/<sampl_id>/plot_batch_summary.svg`.
    """

    # Initialize output dictionary
    paths_to_plots = {}

    # For every sample make a plot
    for sample_id, path_to_summary in path_to_run_summaries.items():
        # Read in text file as data frame
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
            plots_dir, sample_id, "plot_batch_summary.svg"
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

    # Get shortest common path between files
    common_path = os.path.commonpath(paths_to_plots.values())

    # Write to zipfile
    with ZipFile(zip_path, "w") as zf:
        for _, path_to_plot in paths_to_plots.items():
            arcname = os.path.relpath(path_to_plot, common_path)
            zf.write(path_to_plot, arcname=arcname)


def _collect_summaries_and_save(
        all_summaries_path: str,
        path_to_run_summaries: dict
        ) -> pd.DataFrame:
    """
    Reads-in the sample wise summaries and concatenates them it one
    pd.DataFrame, which is saved to file.

    Args:
        all_summaries_path (str): Directory path where to write the 
            pd.DataFrame 
        path_to_run_summaries (dict): dict where key is sample id
        and value is path "tmp/sample_id/batch_summary.txt"

    Returns:
        all_summaries_df (pd.DataFrame): Data frame composed of the individual
        run summaries.
    """

    all_summaries_list = []
    for sample_id, path_to_summary in path_to_run_summaries.items():
        df = pd.read_csv(filepath_or_buffer=path_to_summary, sep="\t")
        df["sample_id"] = sample_id
        all_summaries_list.append(df)

    # Concatenate
    all_summaries_df = pd.concat(all_summaries_list, ignore_index=True)

    # Save to file
    all_summaries_df.to_csv(all_summaries_path, index=False)

    return all_summaries_df


def _render_html(
        output_dir: str,
        all_summaries_df: pd.DataFrame,
        ):
    """
    Renders an qiime2 html file with the plots summarizing the BUSCO output.

    Args:
        output_dir (str): Directory path where to write the pd.DataFrame
        all_summaries_df (pd.DataFrame): Data frame composed of the individual
            run summaries.
    """
    # Prepare context for jinja2 template
    context = {
        "tabs": [{"title": "BUSCO Results", "url": "index.html"}],
        "vega_plots_overview": _draw_busco_plots_for_render(
            all_summaries_df,
            width=600,
            height=18,
            titleFontSize=20,
            labelFontSize=17,
        ),
    }

    # Copy BUSCO results from tmp dir to output_dir
    moshpit_path = os.path.dirname(  # Path to parent dir, q2_moshpit
        os.path.dirname(__file__)
    )
    TEMPLATES = os.path.join(moshpit_path, "assets")
    index = os.path.join(TEMPLATES, "busco", "index.html")
    copytree(
        src=os.path.join(TEMPLATES, "busco"),
        dst=output_dir,
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
