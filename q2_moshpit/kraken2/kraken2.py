# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import subprocess
import tempfile
from copy import deepcopy
from dataclasses import dataclass, field, make_dataclass
from functools import reduce
from typing import Union

import pandas as pd

from q2_moshpit._utils import _process_common_input_params, run_command
from q2_types_genomics.per_sample_data import MultiMAGSequencesDirFmt

from q2_moshpit.kraken2.utils import _process_kraken2_arg


RANKS = {x: y for x, y in zip('dpcofgs', range(7))}
LEVELS = {x: y for x, y in zip(range(7), 'dpcofgs')}


@dataclass
class Node:
    rank: str
    name: str
    count: int
    level: int
    parent_taxonomy: str = None
    taxonomy: str = None
    children: Union[list, None] = field(default_factory=list)

    def __post_init__(self):
        self.taxonomy = f'{self.rank}__{self.name}'


# prepare taxonomic data classes
LEVEL_CLASSES = []
for i in range(7):
    base_cls = (Node,) if i == 0 else (LEVEL_CLASSES[i-1],)
    LEVEL_CLASSES.append(make_dataclass(
        cls_name=f'L{i}',
        fields=[],
        bases=base_cls
    ))


def _find_parent(all_levels: dict, parent_rank: str, parent_name: str):
    parent_candidates = [
        x for x in all_levels[parent_rank] if x.name == parent_name
    ]
    return parent_candidates[0]


def _parse_kraken2_report(
        report_fp: str, sample_name: str, bin_name: str
) -> pd.DataFrame:
    with open(report_fp, 'r') as r:
        lines = r.readlines()

    all_taxonomies = {}
    levels = {x: [] for x in 'dpcofgs'}
    for line in lines:
        line = line.split('\t')
        taxonomy, count = line[0], int(line[1].strip())
        all_taxonomies[taxonomy] = {f'{sample_name}/{bin_name}': count}

        line = [tuple(x.split('__')) for x in line[0].split('|')]
        for i, (rank, name) in enumerate(line):
            existing_ranks = [node.name for node in levels[rank]]

            # does it exist already?
            if name not in existing_ranks:
                current_level = RANKS[rank]
                current_node = LEVEL_CLASSES[current_level](
                    rank=rank, name=name, count=count, level=current_level
                )
                # add the node to the list of ranks
                levels[rank].append(current_node)
                if rank != 'd':
                    (parent_rank, parent_name) = line[i - 1]
                    parent = _find_parent(levels, parent_rank, parent_name)
                    # add the node to the list of children of a higher rank
                    parent.children.append(current_node)
                    # add the parent to the current node
                    current_node.parent_taxonomy = parent

    for _l in range(5, -1, -1):
        current_rank = LEVELS[_l]
        for node in levels[current_rank]:
            child_rank, child_level = LEVELS[_l + 1], _l + 1
            children_counts = sum([x.count for x in node.children])
            if children_counts < node.count:
                unclassified_node = LEVEL_CLASSES[child_level](
                    rank=child_rank, name='unclassified',
                    count=node.count - children_counts, level=child_level,
                    children=None, parent_taxonomy=node)
                node.children.append(unclassified_node)
                levels[child_rank].append(unclassified_node)

    # expand taxonomies
    for _l in range(7):
        current_rank = LEVELS[_l]
        for node in levels[current_rank]:
            if _l > 0:
                node.taxonomy = f'{node.parent_taxonomy.taxonomy};' \
                                f'{node.taxonomy}'

    # find childless nodes
    features = []
    for _, taxonomies in levels.items():
        childless = [x for x in taxonomies if not x.children]
        features.extend(childless)

    abundances = {x.taxonomy: x.count for x in features}
    abundances = pd.DataFrame.from_dict(
        abundances, orient='index', columns=[f'{sample_name}/{bin_name}']
    )
    return abundances


def _process_kraken2_reports(reports: dict) -> pd.DataFrame:
    # get abundances per sample
    abundances = []
    for _sample, bins in reports.items():
        sample_abundances = []
        for _bin, reports in bins.items():
            bin_abundances = _parse_kraken2_report(
                reports['report'], _sample, _bin
            )
            sample_abundances.append(bin_abundances)
        # merge bins into one per sample
        bins_merged = reduce(
            lambda left, right: pd.merge(
                left, right, left_index=True, right_index=True, how='outer'
            ), sample_abundances
        )
        bins_merged = bins_merged.fillna(0).sum(axis=1)
        bins_merged.name = _sample
        abundances.append(bins_merged)

    # combine all samples
    df_merged = reduce(
        lambda left, right: pd.merge(
            left, right, left_index=True, right_index=True, how='outer'
        ), abundances
    )

    df_merged.index.name = 'feature-id'

    return df_merged


def _fill_missing_levels(
        taxon: str, sep: str =';',
        fill_with: str = 'unclassified',
        prepend: bool =True
) -> str:
    ranks = 'dpcofgs'
    taxon = taxon.split(sep)
    for i in range(len(taxon), 7):
        level = f'{ranks[i]}__{fill_with}' if prepend else fill_with
        taxon.append(level)
    return sep.join(taxon)


def _classify_kraken(manifest, common_args) -> (pd.DataFrame, pd.DataFrame):
    base_cmd = ["kraken2", *common_args]
    kraken2_reports = {}

    with tempfile.TemporaryDirectory() as tmp_dir:
        try:
            for (_sample, _bin), fn in manifest.itertuples():
                if _sample not in kraken2_reports:
                    kraken2_reports[_sample] = {}
                cmd = deepcopy(base_cmd)
                bin_dir = os.path.join(tmp_dir, _sample, _bin)
                os.makedirs(bin_dir, exist_ok=True)
                report_fp = os.path.join(bin_dir, 'report.txt')
                output_fp = os.path.join(bin_dir, 'output.txt')
                cmd.extend([
                    '--use-mpa-style', '--report', report_fp, '--use-names',
                    '--output', output_fp, fn
                ])
                run_command(cmd=cmd, verbose=True)
                kraken2_reports[_sample].update(
                    {_bin: {'report': report_fp, 'output': output_fp}}
                )
        except subprocess.CalledProcessError as e:
            raise Exception(
                "An error was encountered while running Kraken 2, "
                f"(return code {e.returncode}), please inspect "
                "stdout and stderr to learn more."
            )
        results_df = _process_kraken2_reports(kraken2_reports)

        # fill missing levels
        results_df.index = results_df.index.map(_fill_missing_levels)
        results_df.fillna(0, inplace=True)

        # create "fake" taxonomy for now
        # should we use the lowest level's NCBI ID here?
        taxonomy = pd.DataFrame(
            results_df.index.tolist(),
            index=results_df.index,
            columns=['Taxon']
        )
        taxonomy.index.name = 'Feature ID'
        taxonomy['Taxon'] = taxonomy['Taxon'].apply(
            lambda x: x.replace('|', ';')
        )

    return results_df.T, taxonomy


def classify_kraken(
        seqs: MultiMAGSequencesDirFmt, db: str, threads: int = 1,
        confidence: float = 0.0, minimum_base_quality: int = 0,
        memory_mapping: bool = False, minimum_hit_groups: int = 2,
        quick: bool = False
) -> (pd.DataFrame, pd.DataFrame):
    kwargs = {
        k: v for k, v in locals().items()
        if k not in ["seqs"]
    }
    common_args = _process_common_input_params(
        processing_func=_process_kraken2_arg, params=kwargs
    )
    manifest: pd.DataFrame = seqs.manifest.view(pd.DataFrame)

    return _classify_kraken(manifest, common_args)
