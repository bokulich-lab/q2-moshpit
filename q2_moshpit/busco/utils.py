import pandas as pd
import altair as alt
import json
from types import List

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
    df: pd.core.frame.DataFrame, width: int = None, height: int = None
) -> dict:
    """
    Draws a hroizontal normalized bar plot for every sample for which BUSCO was
    run. Each barplot shows the BUSCO results for each of the MAGs in the sample.
    The plots for all samples are drwan in one composite plot which is then returned as
    a dictionary for rendering.

    Args:
        df (pd.core.frame.DataFrame): tabular batch summary for all samples
        width (int): width of the plot
        height (int): height of the plot

    Output:
        Output plot in dictionary from.
    """

    # Clean column names
    df.columns = df.columns.str.replace(" ", "_")
    df.columns = df.columns.str.replace("[^a-zA-Z0-9_]", "", regex=True)
    df.columns = df.columns.str.lower()

    # Rename column "input_file"
    df["mag_id"] = df["input_file"].str.split(".", expand=True)[0]

    # Pivot long
    df2 = pd.melt(
        df2,
        id_vars=["sample_id", "mag_id"],
        value_vars=["single", "duplicated", "fragmented", "missing"],
        value_name="BUSCO_percentage",
        var_name="category",
    )

    # Specify order
    mapping = {"single": 1, "duplicated": 2, "fragmented": 3, "missing": 4}
    df2["order"] = df2["category"].map(mapping)

    # Plot
    domain = ["single", "duplicated", "fragmented", "missing"]
    range_ = ["#1E90FF", "#87CEFA", "#FFA500", "#FF7F50"]

    output_plot = (
        alt.Chart(df2)
        .mark_bar()
        .encode(
            x=alt.X("sum(BUSCO_percentage)", stack="normalize", title="BUSCO fracc."),
            y=alt.Y("mag_id", axis=alt.Axis(title="MAG ID")),
            # color=alt.Color('category').scale(domain=domain, range=range_),
            color=alt.Color(
                "category",
                scale=alt.Scale(domain=domain, range=range_),
                legend=alt.Legend(title="BUSCO Category"),
            ),
            order=alt.Order(
                # Sort the segments of the bars by this field
                "order",
                sort="ascending",
            ),
            row=alt.Row("sample_id", title="Sample ID"),
        )
        .configure_axis(labelFontSize=12, titleFontSize=15)
    )

    # Return
    if width and height:
        return json.dumps(output_plot.properties(width=width, height=height).to_dict())
    else:
        return json.dumps(output_plot.to_dict())
