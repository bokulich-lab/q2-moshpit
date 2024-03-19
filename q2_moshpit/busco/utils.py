import json
import os
import q2templates
from shutil import copytree
import pandas as pd
import altair as alt
from typing import List

arguments_with_hyphens = {
    "auto_lineage": "auto-lineage",
    "auto_lineage_euk": "auto-lineage-euk",
    "auto_lineage_prok": "auto-lineage-prok",
    "list_datasets": "list-datasets",
    "update_data": "update-data",
}


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


def _draw_busco_plots(
    df: pd.DataFrame,
    width: int = None,
    height: int = None,
    labelFontSize: int = None,
    titleFontSize: int = None,
    spacing: int = None
) -> str:
    """
    Draws a horizontal normalized bar plot for every sample for which BUSCO was
    run. Each barplot shows the BUSCO results for each of the MAGs in the
    sample. The plots for all samples are drawn in one composite plot which
    is then returned as a dictionary for rendering (but casted to a string).

    Args:
        df (pd.DataFrame): tabular batch summary for all samples
        width (int): width of the plot
        height (int): height of each bar in the plot
        labelFontSize (int): size of the labels in plot
        titleFontSize (int): size of titles in plot

    Output:
        Output plot in dictionary from casted to a string.
    """
    # Format data for plotting
    df = df.reset_index(drop=False, inplace=False)
    df = df.rename(columns={"id": "mag_id"}, inplace=False)

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
        'scaffolds',
    ]]

    # Specify order
    mapping = {"single": 1, "duplicated": 2, "fragmented": 3, "missing": 4}
    busco_plot_data["order"] = busco_plot_data["category"].map(mapping)

    # Estimate fraction of sequences in each BUSCO category
    busco_plot_data["frac_markers"] = (
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

    # Get the first 10 sample ids
    # if len(df['sample_id'].unique()) <= 10:
    #     default_regex = ""
    # else:
    #     default_regex = df['sample_id'].unique()[0:10]
    #     default_regex = '$|^'.join(default_regex)
    #     default_regex = '^' + default_regex + "$"

    # Make BUSCO bar plots (the plots on the left)
    busco_plot = (
        alt.Chart(busco_plot_data)
        .mark_bar()
        .encode(
            x=alt.X(
                "sum(BUSCO_percentage)",
                stack="normalize",
                title="BUSCO fraction"
            ),
            y=alt.Y("mag_id", axis=alt.Axis(titleFontSize=0)),
            color=alt.Color(
                "category",
                scale=alt.Scale(domain=domain, range=range_),
                legend=alt.Legend(title="BUSCO category", orient="top"),
            ),
            order=alt.Order("order", sort="ascending"),
            tooltip=[
                alt.Tooltip("sample_id", title="Sample ID"),
                alt.Tooltip("mag_id", title="MAG ID"),
                alt.Tooltip("dataset", title="Lineage dataset"),
                alt.Tooltip(
                    "frac_markers",
                    title="Approx. number of markers in this category"
                ),
                alt.Tooltip("BUSCO_percentage", title="% BUSCOs"),
            ],
            opacity=alt.value(0.85),
        )
        .properties(
            width=width,
            height={"step": height}
        )
        .facet(
            row=alt.Row(
                "sample_id",
                title="Sample ID / MAG ID"
            ),
            spacing=spacing
        )
        .resolve_scale(y="independent")
    )

    # Secondary plot
    # Drop down menu
    dropdown = alt.binding_select(
        options=[
            'scaffold_n50',
            'contigs_n50',
            'percent_gaps',
            'scaffolds',
        ],
        name="Assembly statistics: "
    )

    xcol_param = alt.param(
        value='scaffold_n50',
        bind=dropdown
    )

    secondary_plot = alt.Chart(secondary_plot_data).mark_bar().encode(
        x=alt.X('x:Q').title('Assembly statistic'),
        y=alt.Y('mag_id', axis=None),
        tooltip=[alt.Tooltip('x:Q', title="value")],
        opacity=alt.value(0.85)
    ).transform_calculate(
        x=f'datum[{xcol_param.name}]'
    ).add_params(
        xcol_param
    ).properties(
        width=width,
        height={"step": height}
    ).facet(
        row=alt.Row(
            "sample_id",
            title=None,
            header=alt.Header(labelFontSize=0),
        ),
        spacing=spacing
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
    return output_plot.to_dict()


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

    # Concatenate
    all_summaries = pd.concat(all_summaries, ignore_index=True)

    return all_summaries


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
    dfs = _partition_dataframe(all_summaries_df, max_rows=100)

    context = {}
    counter_left = 1
    for i, df in enumerate(dfs):
        sample_count = df['sample_id'].nunique()
        counter_right = counter_left + sample_count - 1
        sample_counter = {"from": counter_left, "to": counter_right}
        counter_left += sample_count
        df = _parse_df_columns(df)
        subcontext = _draw_busco_plots(
            df,
            width=600,
            height=30,
            titleFontSize=20,
            labelFontSize=17,
            spacing=20
        )
        context.update(
            {f"sample{i}": {
                "subcontext": subcontext,
                "sample_counter": sample_counter,
                "sample_ids": df['sample_id'].unique().tolist(),
            }}
        )

    vega_out_fp = os.path.join(output_dir, "vega.json")
    with open(vega_out_fp, 'w') as json_file:
        vega_json = json.dumps(context)
        json_file.write(vega_json)

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
    q2templates.render(index, output_dir, context={"vega_json": vega_json})

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


def _parse_df_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds several columns required for generation of downloadable 
    BUSCO plots.

    Args:
        df (pd.DataFrame): Unformatted DataFrame

    Returns:
        df (pd.DataFrame): Formatted DataFrame
    """
    # Cast into percent_gaps col to float
    df["percent_gaps"] = df["percent_gaps"].str.split(
        '%', expand=True
    )[0].map(float)

    for col in ["single", "duplicated", "fragmented", "missing"]:
        df[col] = df[col].map(float)

    df["n_markers"] = df["n_markers"].map(int)

    # Make new columns for downloadable plots
    # (only used in _draw_busco_plots)
    df["single_"] = df["single"]
    df["duplicated_"] = df["single_"] + df["duplicated"]
    df["fragmented_"] = df["duplicated_"] + df["fragmented"]
    df["missing_"] = df["fragmented_"] + df['missing']

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
