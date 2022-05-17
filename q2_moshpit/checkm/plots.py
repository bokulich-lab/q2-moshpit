# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import altair as alt
import pandas as pd


def _draw_detailed_plots(df: pd.DataFrame) -> dict:
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
    final_plot = _concatenate_detailed_plots(
        completeness_plot=_prep_scatter_plot(
            base, 'completeness', 'contamination', 'Completeness [%]',
            'Contamination [%]', primary_selection=bin_selection,
            primary_filter=sample_selection, interactive=False,
            more_tooltips=[
                alt.Tooltip('marker_lineage:N', title='Marker lineage')
            ],
            reverse_y=True
        ),
        gc_plot=_prep_scatter_plot(
            base, 'gc', 'coding_density', 'GC content', 'Coding density',
            primary_selection=bin_selection, primary_filter=sample_selection,
            interactive=False
        ),
        marker_plot=_prep_bar_plot(
            base, sample_selection, x_col='bin_id', y_col='sum(value):Q',
            x_title='Bin ID', y_title='Marker count', color_shorthand='key:N',
            color_title='Markers', color_map='blues', sort='ascending',
            bin_selection=bin_selection, width=880, height=250
        ),
        contig_plots=_prep_contig_plots(base, sample_selection, bin_selection),
        genes_plot=_prep_scatter_plot(
            base, 'genome_size', 'predicted_genes', 'Genome size [Mbp]',
            'Predicted genes', primary_selection=bin_selection,
            primary_filter=sample_selection, interactive=False
        ),
        contig_count_plot=_prep_scatter_plot(
            base, 'contigs', 'genome_size', 'Contigs', 'Genome size [Mbp]',
            primary_selection=bin_selection, primary_filter=sample_selection,
            interactive=False
        )
    )

    return final_plot.to_dict()


def _draw_overview_plots(df: pd.DataFrame) -> dict:
    # rename columns for better plot labels
    col_names = {
        'count0': '0', 'count1': '1', 'count2': '2',
        'count3': '3', 'count4': '4', 'count5_or_more': '5+'
    }
    df = df.rename(columns=col_names, inplace=False)

    # convert genome size to Mbp
    df['genome_size'] = df['genome_size'] / 10**6

    # prepare required selectors
    sample_selection = alt.selection_interval()

    base = alt.Chart(df).transform_fold(list(col_names.values()))

    # prep and concatenate all plots
    final_plot = _concatenate_overview_plots(
        completeness_samples_plot=_prep_scatter_plot(
            base, 'completeness', 'contamination', 'Completeness [%]',
            'Contamination [%]', primary_selection=sample_selection,
            selection_col='sample_id:N', selection_title='Sample ID',
            primary_filter=None, interactive=False, reverse_y=True
        ),
        completeness_contigs_plot=_prep_scatter_plot(
            base, 'completeness', 'contamination', 'Completeness [%]',
            'Contamination [%]', primary_selection=sample_selection,
            selection_col='sum(contigs):Q', selection_title='Contig count',
            color_map='blues', primary_filter=None, interactive=False,
            reverse_y=True, more_tooltips=[
                alt.Tooltip('sum(contigs):Q', title='Contigs')
            ]
        ),
        completeness_summary_plot=_prep_bar_plot(
            base, sample_selection, x_col='sample_id', y_col='count(bin_id):Q',
            x_title='Sample ID', y_title='Bin count',
            color_shorthand='qc_category:N', color_title='Genome completeness',
            color_map='blues', sort='ascending', bin_selection=None,
            width=900, height=350
        )
    )

    return final_plot.to_dict()


def _concatenate_detailed_plots(
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


def _concatenate_overview_plots(
        completeness_samples_plot, completeness_contigs_plot,
        # coverage_plot,
        completeness_summary_plot
):
    plot = alt.vconcat(
        alt.hconcat(
            completeness_samples_plot, completeness_contigs_plot, spacing=40
        ).resolve_scale(color='independent'),
        # alt.hconcat(
        #     coverage_plot, completeness_summary_plot, spacing=40
        # ),
        completeness_summary_plot,
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
        primary_selection, selection_col: str = 'bin_id:N',
        selection_title: str = 'Bin ID', primary_filter=None,
        more_tooltips=None, color_map='tableau10', interactive=True,
        reverse_x=False, reverse_y=False
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
            primary_selection,
            alt.Color(
                selection_col,
                scale=alt.Scale(scheme=color_map),
                legend=alt.Legend(title=selection_title)
            ),
            alt.value('gainsboro')
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
        primary_selection
    )

    if primary_filter:
        plot = plot.add_selection(
            primary_filter
        ).transform_filter(
            primary_filter
        )

    plot = plot.properties(width=400, height=400)
    return plot.interactive() if interactive else plot


def _prep_bar_plot(
        base_plot, sample_selection, x_col: str, y_col: str, x_title: str,
        y_title: str, color_shorthand: str, color_title: str, color_map: str,
        sort: str = 'ascending', bin_selection=None, width=880, height=250
):
    plot_markers = base_plot.mark_bar().encode(
        x=alt.X(x_col, title=x_title),
        y=alt.Y(y_col, title=y_title),
        color=alt.Color(
            color_shorthand,
            title=color_title,
            scale=alt.Scale(scheme=color_map),
            sort=sort
        )
    ).transform_filter(
        sample_selection
    )

    if bin_selection:
        plot_markers = plot_markers.transform_filter(
            bin_selection
        )

    return plot_markers.properties(width=width, height=height)
