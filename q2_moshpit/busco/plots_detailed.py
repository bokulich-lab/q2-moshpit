# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import altair as alt
import pandas as pd

alt.data_transformers.disable_max_rows()


def _draw_detailed_plots(
    df: pd.DataFrame,
    width: int = None,
    height: int = None,
    labelFontSize: int = None,
    titleFontSize: int = None,
    spacing: int = None
) -> dict:
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

    return output_plot.to_dict()
