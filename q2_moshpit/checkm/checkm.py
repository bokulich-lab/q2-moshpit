# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import json
import os
import tempfile
from copy import deepcopy
from distutils.dir_util import copy_tree
from zipfile import ZipFile

import altair as alt
import pandas as pd
import pkg_resources
import q2templates

from q2_moshpit._utils import _process_common_input_params, run_command
from q2_moshpit.checkm.utils import _process_checkm_arg, _get_plots_per_sample

from q2_types_genomics.per_sample_data._format import \
    (MultiMAGSequencesDirFmt)

TEMPLATES = pkg_resources.resource_filename('q2_moshpit', 'assets')


def _evaluate_bins(
        results_dir: str, bins: MultiMAGSequencesDirFmt,
        db_path: str, common_args: list
) -> dict:
    base_cmd = ['checkm', 'lineage_wf', *common_args]
    stats_fps = {}

    manifest: pd.DataFrame = bins.manifest.view(pd.DataFrame)
    manifest['sample_dir'] = manifest.filename.apply(
        lambda x: os.path.dirname(x)
    )
    sample_dirs = manifest['sample_dir'].unique()
    for sample_dir in sample_dirs:
        sample = os.path.split(sample_dir)[-1]
        sample_results = os.path.join(results_dir, sample)

        cmd = deepcopy(base_cmd)
        cmd.extend(['-x', 'fasta', sample_dir, sample_results])
        run_command(cmd, env={**os.environ, 'CHECKM_DATA_PATH': db_path})

        stats_fp = os.path.join(sample_results, 'storage', 'bin_stats_ext.tsv')
        if os.path.isfile(stats_fp):
            stats_fps[sample] = stats_fp
        else:
            raise FileNotFoundError(
                f'CheckM stats file {stats_fp} could not be found.'
            )

    return stats_fps


def _draw_checkm_plots(
        results_dir: str, bins: MultiMAGSequencesDirFmt,
        db_path: str, plot_type: str = 'gc'
) -> dict:
    base_cmd = [
        'checkm', f'{plot_type}_plot', '-x', 'fasta', '--image_type', 'svg',
        '--font_size', '10'
    ]
    # TODO: the numbers should probably be configurable
    dist_values = [] if plot_type == 'nx' else ['50', '75', '90']
    # width = '6.5' if plot_type == 'nx' else '6.5'
    # height = '6.5' if plot_type == 'nx' else '3.5'
    # base_cmd.extend(['--width', width, '--height', height])
    plots = {}

    manifest: pd.DataFrame = bins.manifest.view(pd.DataFrame)
    manifest['sample_dir'] = manifest.filename.apply(
        lambda x: os.path.dirname(x)
    )
    sample_dirs = manifest['sample_dir'].unique()
    for sample_bins in sample_dirs:
        sample = os.path.split(sample_bins)[-1]
        sample_plots = os.path.join(results_dir, 'plots', plot_type, sample)
        checkm_files = os.path.join(results_dir, sample)
        plots[sample] = sample_plots

        cmd = deepcopy(base_cmd)
        cmd.append(checkm_files) if plot_type == 'coding' else False
        cmd.extend([sample_bins, sample_plots, *dist_values])
        run_command(cmd, env={**os.environ, 'CHECKM_DATA_PATH': db_path})
    return plots


def _parse_checkm_reports(reports: dict) -> pd.DataFrame:
    dfs = [
        _parse_single_checkm_report(_id, fp) for _id, fp in reports.items()
    ]
    results_df = pd.concat(dfs)
    results_df.reset_index(drop=True, inplace=True)
    return results_df


def _parse_single_checkm_report(
        sample_id: str, report_fp: str
) -> pd.DataFrame:
    # read the raw CheckM report
    with open(report_fp, 'r') as fh:
        stats = {
            k: json.loads(v.replace('\'', '"')) for [k, v] in
            [line.split('\t') for line in fh.readlines()]
        }

    # convert report to DataFrame
    df = pd.DataFrame.from_dict(stats, orient='index')
    df.reset_index(drop=False, inplace=True)
    col_names = {
        'index': 'bin_id',
        'marker lineage': 'marker_lineage', '# genomes': 'genomes',
        '# markers': 'markers', '# marker sets': 'marker_sets',
        '0': 'count0', '1': 'count1', '2': 'count2', '3': 'count3',
        '4': 'count4', '5+': 'count5_or_more', 'Completeness': 'completeness',
        'Contamination': 'contamination', 'GC': 'gc', 'GC std': 'gc_std',
        'Genome size': 'genome_size', '# ambiguous bases': 'ambiguous_bases',
        '# scaffolds': 'scaffolds', '# contigs': 'contigs',
        'Longest scaffold': 'longest_scaffold',
        'Longest contig': 'longest_contig',
        'N50 (scaffolds)': 'n50_scaffolds', 'N50 (contigs)': 'n50_contigs',
        'Mean scaffold length': 'mean_scaffold_length',
        'Mean contig length': 'mean_contig_length',
        'Coding density': 'coding_density',
        'Translation table': 'translation_table',
        '# predicted genes': 'predicted_genes',
        'GCN0': 'gcn0', 'GCN1': 'gcn1', 'GCN2': 'gcn2', 'GCN3': 'gcn3',
        'GCN4': 'gcn4', 'GCN5+': 'gcn5_or_more'
    }
    df.rename(columns=col_names, inplace=True)
    df['sample_id'] = sample_id

    # reorder columns
    df = df[['sample_id', *col_names.values()]]

    return df


def _draw_plots(df: pd.DataFrame) -> dict:
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

    opacity_selection = alt.selection_single(
        fields=['opacity'],
        bind=alt.binding_range(
            min=0.1, max=1, name='Opacity: other samples  '
        ),
        init={'opacity': 0.15}
    )

    base = alt.Chart(df).transform_calculate(
        opacity=opacity_selection.opacity
    ).transform_fold(
        list(col_names.values())
    )

    # prep and concatenate all plots
    final_plot = _concatenate_plots(
        completeness_plot=_prep_completeness_plot(
            base, opacity_selection, sample_selection
        ),
        gc_plot=_prep_gc_plot(base, opacity_selection, sample_selection),
        marker_plot=_prep_markers_plot(base, sample_selection),
        contig_plots=_prep_contig_plots(base, sample_selection),
        genes_plot=_prep_genes_plot(base, opacity_selection, sample_selection),
        contig_count_plot=_prep_contig_count_plot(
            base, opacity_selection, sample_selection
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


def _prep_contig_plots(base_plot, sample_selection):
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
            ).add_selection(
                sample_selection
            ).transform_filter(
                sample_selection
            ).properties(
                width=880, height=250
            )
    return contig_plots


def _prep_gc_plot(base_plot, opacity_selection, sample_selection):
    plot_gc = base_plot.mark_point(
        size=80, filled=True, opacity=0.7,
        strokeWidth=0.2, strokeOpacity=0.8, stroke='lightgrey'
    ).encode(
        x=alt.X('gc:Q', title='GC content'),
        y=alt.Y('coding_density:Q', title='Coding density'),
        color=alt.Color(
            'sample_id:N',
            scale=alt.Scale(scheme='tableau10'),
            legend=alt.Legend(title='Sample ID')
        ),
        tooltip=[
            alt.Tooltip('sample_id:N', title='Sample ID'),
            alt.Tooltip('bin_id:N', title='Bin ID'),
            alt.Tooltip('coding_density:Q', title='Coding density'),
            alt.Tooltip('gc:Q', title='GC content'),
        ],
        opacity=alt.condition(
            sample_selection, alt.value(1),
            alt.Opacity('opacity:Q', scale=None)
        )
    ).add_selection(
        sample_selection, opacity_selection
    ).properties(
        width=400, height=400
    ).interactive()
    return plot_gc


def _prep_genes_plot(base_plot, opacity_selection, sample_selection):
    plot_genes = base_plot.mark_point(
        size=80, filled=True, opacity=0.7,
        strokeWidth=0.2, strokeOpacity=0.8, stroke='lightgrey'
    ).encode(
        x=alt.X('genome_size:Q', title='Genome size [Mbp]'),
        y=alt.Y('predicted_genes:Q', title='Count of predicted genes'),
        color=alt.Color(
            'sample_id:N',
            scale=alt.Scale(scheme='tableau10'),
            legend=alt.Legend(title='Sample ID')
        ),
        tooltip=[
            alt.Tooltip('sample_id:N', title='Sample ID'),
            alt.Tooltip('bin_id:N', title='Bin ID'),
            alt.Tooltip('genome_size:Q', title='Genome size [Mbp]'),
            alt.Tooltip('predicted_genes:Q', title='Predicted genes'),
        ],
        opacity=alt.condition(
            sample_selection, alt.value(1),
            alt.Opacity('opacity:Q', scale=None)
        )
    ).add_selection(
        sample_selection, opacity_selection
    ).properties(
        width=400, height=400
    ).interactive()
    return plot_genes


def _prep_contig_count_plot(base_plot, opacity_selection, sample_selection):
    plot_contigs = base_plot.mark_point(
        size=80, filled=True, opacity=0.7,
        strokeWidth=0.2, strokeOpacity=0.8, stroke='lightgrey'
    ).encode(
        x=alt.X('contigs:Q', title='Count of contigs'),
        y=alt.Y('genome_size:Q', title='Genome size [Mbp]'),
        color=alt.Color(
            'sample_id:N',
            scale=alt.Scale(scheme='tableau10'),
            legend=alt.Legend(title='Sample ID')
        ),
        tooltip=[
            alt.Tooltip('sample_id:N', title='Sample ID'),
            alt.Tooltip('bin_id:N', title='Bin ID'),
            alt.Tooltip('genome_size:Q', title='Genome size [Mbp]'),
            alt.Tooltip('contigs:Q', title='Contigs'),
        ],
        opacity=alt.condition(
            sample_selection, alt.value(1),
            alt.Opacity('opacity:Q', scale=None)
        )
    ).add_selection(
        sample_selection, opacity_selection
    ).properties(
        width=400, height=400
    ).interactive()
    return plot_contigs


def _prep_markers_plot(base_plot, sample_selection):
    plot_markers = base_plot.mark_bar().encode(
        x=alt.X('bin_id', title='Bin'),
        y=alt.Y('sum(value):Q', title='Marker count'),
        color=alt.Color(
            'key:N',
            title='Markers',
            scale=alt.Scale(scheme='blues'),
            sort='ascending'
        )
    ).add_selection(
        sample_selection
    ).transform_filter(
        sample_selection
    ).properties(
        width=880, height=250
    )
    return plot_markers


def _prep_completeness_plot(base_plot, opacity_selection, sample_selection):
    plot_completeness = base_plot.mark_point(
        size=80, filled=True, opacity=0.7,
        strokeWidth=0.2, strokeOpacity=0.8, stroke='lightgrey'
    ).encode(
        x=alt.X('completeness:Q', title='Completeness [%]'),
        y=alt.Y('contamination:Q', title='Contamination [%]'),
        color=alt.Color(
            'sample_id:N',
            scale=alt.Scale(scheme='viridis'),
            legend=alt.Legend(title='Sample ID')
        ),
        tooltip=[
            alt.Tooltip('sample_id:N', title='Sample ID'),
            alt.Tooltip('bin_id:N', title='Bin ID'),
            alt.Tooltip('marker_lineage:N', title='Marker lineage'),
            alt.Tooltip('completeness:Q', title='% completeness'),
            alt.Tooltip('contamination:Q', title='% contamination')
        ],
        opacity=alt.condition(
            predicate=sample_selection,
            if_true=alt.value(1),
            if_false=alt.Opacity('opacity:Q', scale=None)
        )
    ).add_selection(
        sample_selection, opacity_selection
    ).properties(
        width=400, height=400
    ).interactive()
    return plot_completeness


def evaluate_bins(
    output_dir: str, bins: MultiMAGSequencesDirFmt,
    db_path: str, reduced_tree: bool = None, unique: int = None,
    multi: int = None, force_domain: bool = None, no_refinement: bool = None,
    individual_markers: bool = None, skip_adj_correction: bool = None,
    skip_pseudogene_correction: bool = None, aai_strain: float = None,
    ignore_thresholds: bool = None, e_value: float = None,
    length: float = None, threads: int = None, pplacer_threads: int = None
):

    kwargs = {k: v for k, v in locals().items()
              if k not in ['output_dir', 'bins', 'db_path']}
    common_args = _process_common_input_params(
        processing_func=_process_checkm_arg, params=kwargs
    )

    # TODO: check that CheckM's database is available (or fetch?)

    with tempfile.TemporaryDirectory() as tmp:
        results_dir = os.path.join(tmp, 'results')

        # run CheckM's lineage_wf pipeline and draw all the QC plots
        reports = _evaluate_bins(results_dir, bins, db_path, common_args)
        all_plots = {}
        for plot_type in ['gc', 'nx', 'coding']:
            plot_dirs = _draw_checkm_plots(
                results_dir, bins, db_path, plot_type=plot_type
            )
            all_plots[f'plots_{plot_type}'] = plot_dirs

        # TODO: move this to a new function and test
        plots_per_sample = _get_plots_per_sample(all_plots)

        # TODO: instead of displaying the ugly plots directly in the viz
        #  only provide a link to the zipped plots for download - zip them here
        zip_path = os.path.join(TEMPLATES, 'checkm', 'checkm_plots.zip')
        with ZipFile(zip_path, 'w') as zf:
            for sample_id in plots_per_sample.keys():
                plot_fps = [
                    glob.glob(os.path.join(v, '*.svg'))
                    for k, v in plots_per_sample[sample_id].items()
                ]
                plot_fps = [x for sublist in plot_fps for x in sublist]
                common_path = os.path.commonpath(plot_fps)
                for plot_fp in plot_fps:
                    arcname = os.path.relpath(plot_fp, common_path)
                    zf.write(plot_fp, arcname=arcname)

        checkm_results = _parse_checkm_reports(reports)
        checkm_results.to_csv(
            os.path.join(TEMPLATES, 'checkm', 'results.tsv'),
            sep='\t', index=False,
        )

        context = {
            'samples': json.dumps(list(reports.keys())),
            'vega_plots': json.dumps(
                _draw_plots(checkm_results)
            ),
            **{k: json.dumps(v) for k, v in all_plots.items()}
        }

        index = os.path.join(TEMPLATES, 'checkm', 'index.html')
        gc_plots = os.path.join(TEMPLATES, 'checkm', 'gc_plots.html')

        copy_tree(os.path.join(TEMPLATES, 'checkm'), output_dir)
        copy_tree(
            os.path.join(results_dir, 'plots'),
            os.path.join(output_dir, 'plots')
        )

        templates = [index, gc_plots]
        q2templates.render(templates, output_dir, context=context)
