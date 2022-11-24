# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import shutil
import subprocess
import tempfile
from copy import deepcopy
from dataclasses import dataclass, field, make_dataclass
from functools import reduce
from typing import Optional

import pandas as pd

from q2_moshpit._utils import _process_common_input_params, run_command
from q2_types_genomics.kraken2 import Kraken2ReportDirectoryFormat
from q2_types_genomics.per_sample_data import MultiMAGSequencesDirFmt

from q2_moshpit.kraken2.utils import _process_kraken2_arg


COLUMNS = [
        'perc_frags_covered', 'no_frags_covered',
        'no_frags_assigned', 'rank', 'ncbi_tax_id', 'name'
    ]

RANK_CODES = 'URDPCOFGS'
RANKS = {x: y for x, y in zip(RANK_CODES, range(9))}
LEVELS = {x: y for x, y in zip(range(9), RANK_CODES)}


@dataclass
class Node:
    rank: str
    name: str
    count: int
    level: int
    parent_taxonomy: str = None
    taxonomy: str = None
    children: Optional[list] = field(default_factory=list)

    def __post_init__(self):
        self.taxonomy = f'{self.rank}__{self.name}'


# prepare taxonomic data classes
LEVEL_CLASSES = []
for i in range(9):
    base_cls = (Node,) if i == 0 else (LEVEL_CLASSES[i-1],)
    LEVEL_CLASSES.append(make_dataclass(
        cls_name=f'L{i}',
        fields=[],
        bases=base_cls
    ))


def _find_parent(all_levels: dict, df: pd.DataFrame, i: int):
    # reverse the df up to the given level for easier look-up of a parent
    df_rev = df[i::-1]

    # crawl through the "tree" and find the taxon with level smaller by 1
    parent_rank, parent_name, parent_level = None, None, None
    for j, row in df_rev.iterrows():
        parent_rank = df_rev.loc[j - 1]['rank'][0]
        parent_name = df_rev.loc[j - 1]['name']
        parent_level = df_rev.loc[j - 1]['level']
        if parent_level <= df.iloc[i]['level'] - 1:
            break
    if not parent_rank:
        raise ValueError('Could not find parent for %s.', df.iloc[i]['name'])

    parent_candidates = [
        x for x in all_levels[parent_rank] if x.name == parent_name
    ]
    return parent_candidates[0]


def _assign_level(rank: str) -> float:
    """Assign a numeric level value based on the provided rank.

    Args:
        rank (str): taxon's rank provided in one-letter code

    Returns:
        level (float): taxon's taxonomic level. If provided rank was not one
            of the 10 standard ones, returned level will be the one of the
            closest ancestor + 0.5.
    """
    if len(rank) == 1:
        return RANKS[rank]
    else:
        return RANKS[rank[0]] + 0.5


def _parse_kraken2_report(
        report_fp: str, sample_name: str, bin_name: str
) -> pd.DataFrame:
    df = pd.read_csv(report_fp, sep='\t', header=None)
    df.columns = COLUMNS
    df['level'] = df['rank'].apply(_assign_level)
    df['name'] = df['name'].apply(lambda x: x.strip())

    levels = {x: [] for x in RANK_CODES}
    for i, row in df.iterrows():
        existing_ranks = [node.name for node in levels[row['rank'][0]]]

        # create the node if it does not yet exist in our register
        name = row['name'].strip()
        if name not in existing_ranks:
            current_node = LEVEL_CLASSES[int(row.level)](
                rank=row['rank'], name=name,
                count=row['no_frags_covered'], level=row['level']
            )
            # add the node to the list of ranks
            levels[row['rank'][0]].append(current_node)
            if row['rank'][0] not in ['U', 'R']:
                parent = _find_parent(levels, df, i)
                # add the node to the list of children of a higher rank
                parent.children.append(current_node)
                # add the parent to the current node
                current_node.parent_taxonomy = parent

    # add "unclassified" nodes wherever necessary (if current level's
    # count is bigger than the sum of its children's counts)
    for _l in range(7, -1, -1):
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

    # expand taxonomies up to the domain rank
    for _l in range(9):
        current_rank = LEVELS[_l]
        for node in levels[current_rank]:
            if _l > 2:
                node.taxonomy = f'{node.parent_taxonomy.taxonomy};' \
                                f'{node.taxonomy}'

    # find childless nodes - these are the ones we will want to report;
    # only consider the ranks defined by RANK_CODES (ignore the intermediate
    # ones like G1, R1, etc.)
    features = []
    for _, taxonomies in levels.items():
        childless = [
            x for x in taxonomies
            if not x.children and x.parent_taxonomy.rank in RANK_CODES
        ]
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
        taxon: str, sep: str = ';',
        fill_with: str = 'unclassified',
        prepend: bool = True
) -> str:
    """Fill missing levels in the provided taxonomy string, down to level 7.

    Args:
        taxon (str): Full taxon name to be filled.
        sep (str): Character used to separate taxonomic levels in the taxon.
        fill_with (str): String used to fill the missing levels,
            down to level 7.
        prepend (bool): If True will prepend every level with one-letter
            abbreviation of the corresponding rank.

    Returns:
        str:  Full taxonomy filled down to level 7.
    """
    ranks = RANK_CODES[2:]
    taxon = taxon.split(sep)
    for i in range(len(taxon), 7):
        level = f'{ranks[i]}__{fill_with}' if prepend else fill_with
        taxon.append(level)
    return sep.join(taxon)


def _classify_kraken(manifest, common_args) -> (
        pd.DataFrame, pd.DataFrame, Kraken2ReportDirectoryFormat
):
    base_cmd = ["kraken2", *common_args]
    kraken2_reports = {}
    kraken2_reports_dir = Kraken2ReportDirectoryFormat()

    with tempfile.TemporaryDirectory() as tmp_dir:
        try:
            for (_sample, _bin), fn in manifest.itertuples():
                if _sample not in kraken2_reports:
                    kraken2_reports[_sample] = {}
                cmd = deepcopy(base_cmd)
                sample_dir = os.path.join(kraken2_reports_dir.path, _sample)
                os.makedirs(sample_dir, exist_ok=True)
                report_fp = os.path.join(sample_dir, f'{_bin}.report.txt')
                output_fp = os.path.join(sample_dir, f'{_bin}.output.txt')
                cmd.extend([
                    '--report', report_fp, '--use-names',
                    # '--output', output_fp,
                    fn
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

        # kraken2_reports_dir = _copy_reports_to_dirfmt(kraken2_reports)
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

    return results_df.T, taxonomy, kraken2_reports_dir


def classify_kraken(
        seqs: MultiMAGSequencesDirFmt, db: str, threads: int = 1,
        confidence: float = 0.0, minimum_base_quality: int = 0,
        memory_mapping: bool = False, minimum_hit_groups: int = 2,
        quick: bool = False
) -> (pd.DataFrame, pd.DataFrame, Kraken2ReportDirectoryFormat):
    kwargs = {
        k: v for k, v in locals().items()
        if k not in ["seqs"]
    }
    common_args = _process_common_input_params(
        processing_func=_process_kraken2_arg, params=kwargs
    )
    manifest: pd.DataFrame = seqs.manifest.view(pd.DataFrame)

    return _classify_kraken(manifest, common_args)
