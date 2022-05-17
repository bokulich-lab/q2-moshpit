# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import altair as alt
import pandas as pd


def _draw_final_plot(df: pd.DataFrame) -> dict:
    # rename columns for better plot labels
    col_names = {
        'count0': '0', 'count1': '1', 'count2': '2',
        'count3': '3', 'count4': '4', 'count5_or_more': '5+'
    }
    df = df.rename(columns=col_names, inplace=False)

    # convert genome size to Mbp
    df['genome_size'] = df['genome_size'] / 10**6

    # prepare required selectors
    sample_ids = df['sample_id'].unique()
    sample_dropdown = alt.binding_select(
        options=sample_ids,
        name='Sample ID  '
    )
    sample_selection = alt.selection_single(
        fields=['sample_id'],
        bind=sample_dropdown,
        init={'sample_id': sample_ids[0]}
    )

    bin_selection = alt.selection_interval()

    base = alt.Chart(df).transform_fold(list(col_names.values()))

    # prep and concatenate all plots
    final_plot = _concatenate_plots(
        completeness_plot=_prep_scatter_plot(
            base, 'completeness', 'contamination', 'Completeness [%]',
            'Contamination [%]', sample_selection, bin_selection,
            interactive=False,
            more_tooltips=[
                alt.Tooltip('marker_lineage:N', title='Marker lineage')
            ],
            reverse_y=True
        ),
        gc_plot=_prep_scatter_plot(
            base, 'gc', 'coding_density', 'GC content', 'Coding density',
            sample_selection, bin_selection, interactive=False
        ),
        marker_plot=_prep_markers_plot(base, sample_selection, bin_selection),
        contig_plots=_prep_contig_plots(base, sample_selection, bin_selection),
        genes_plot=_prep_scatter_plot(
            base, 'genome_size', 'predicted_genes', 'Genome size [Mbp]',
            'Predicted genes', sample_selection, bin_selection,
            interactive=False
        ),
        contig_count_plot=_prep_scatter_plot(
            base, 'contigs', 'genome_size', 'Contigs', 'Genome size [Mbp]',
            sample_selection, bin_selection, interactive=False
        )
    )

    return final_plot.to_dict()


def _concatenate_plots(
        completeness_plot, gc_plot, marker_plot, contig_plots,
        genes_plot, contig_count_plot
):
    plot = alt.vconcat(
        alt.hconcat(
            completeness_plot, gc_plot, spacing=40
        ),
        alt.hconcat(
            genes_plot, contig_count_plot, spacing=40
        ),
        marker_plot,
        *contig_plots.values(),
        spacing=40
    ).resolve_scale(
        color='independent'
    ).configure_axis(
        labelFontSize=12, titleFontSize=15
    ).configure_legend(
        labelFontSize=12, titleFontSize=14
    )
    return plot


def _prep_contig_plots(base_plot, sample_selection, bin_selection):
    contig_plots = {}
    contig_cols = {
        'longest_contig': 'Longest contig length [bp]',
        'n50_contigs': 'N50 contigs [bp]',
        'mean_contig_length': 'Mean contig length [bp]',
        'ambiguous_bases': 'Count of ambiguous bases'
    }
    for y, title in contig_cols.items():
        contig_plots[y] = \
            base_plot.mark_bar().encode(
                x=alt.X('bin_id:N', title='Bin'),
                y=alt.Y(f'{y}:Q', title=title),
            ).transform_filter(
                sample_selection
            ).transform_filter(
                bin_selection
            ).properties(
                width=880, height=250
            )
    return contig_plots


def _prep_scatter_plot(
        base_plot, x_col: str, y_col: str, x_title: str, y_title: str,
        sample_selection, bin_selection, more_tooltips=None,
        color_map='tableau10', interactive=True, reverse_x=False,
        reverse_y=False
):
    more_tooltips = [] if more_tooltips is None else more_tooltips
    plot = base_plot.mark_point(
        size=80, filled=True, opacity=0.7,
        strokeWidth=0.2, strokeOpacity=0.8, stroke='lightgrey'
    ).encode(
        x=alt.X(
            f'{x_col}:Q', title=x_title, scale=alt.Scale(reverse=reverse_x)
        ),
        y=alt.Y(
            f'{y_col}:Q', title=y_title, scale=alt.Scale(reverse=reverse_y)
        ),
        color=alt.condition(
            bin_selection,
            alt.Color(
                'bin_id:N',
                scale=alt.Scale(scheme=color_map),
                legend=alt.Legend(title='Bin ID')
            ),
            alt.value('lightgray')
        ),
        tooltip=[
            alt.Tooltip('sample_id:N', title='Sample ID'),
            alt.Tooltip('bin_id:N', title='Bin ID'),
            alt.Tooltip(f'{x_col}:Q', title=x_title, format='.2'),
            alt.Tooltip(f'{y_col}:Q', title=y_title, format='.2'),
            *more_tooltips
        ],
        opacity=alt.value(0.75)
    ).add_selection(
        sample_selection, bin_selection
    ).transform_filter(
        sample_selection
    ).properties(
        width=400, height=400
    )
    return plot.interactive() if interactive else plot


def _prep_markers_plot(base_plot, sample_selection, bin_selection):
    plot_markers = base_plot.mark_bar().encode(
        x=alt.X('bin_id', title='Bin'),
        y=alt.Y('sum(value):Q', title='Marker count'),
        color=alt.Color(
            'key:N',
            title='Markers',
            scale=alt.Scale(scheme='blues'),
            sort='ascending'
        )
    ).transform_filter(
        sample_selection
    ).transform_filter(
        bin_selection
    ).properties(
        width=880, height=250
    )
    return plot_markers
