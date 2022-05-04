# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re
from dataclasses import dataclass
from typing import List, Type, Union

import pandas as pd


# TODO: make it frozen to prevent modification
@dataclass
class MarkerCounts:
    count0: int
    count1: int
    count2: int
    count3: int
    count4: int
    count5_or_more: int

    def to_series(self) -> pd.Series:
        index = get_attrs(self)
        return pd.Series([getattr(self, k) for k in index], index=index)


@dataclass
class BinStatistics:
    """Class for storing MAG statistics calculated by CheckM"""
    bin_id: str = None
    marker_lineage: str = None
    marker_lineage_uid: str = None
    genomes: int = None
    markers: int = None
    marker_sets: int = None
    marker_counts: MarkerCounts = None
    completeness: float = None
    contamination: float = None
    strain_heterogeneity: float = None

    # TODO: add __post_init__ to validate number ranges

    def to_series(self) -> pd.Series:
        index = get_attrs(self, ('marker_counts',))
        main_series = pd.Series([getattr(self, k) for k in index], index=index)
        marker_counts = self.marker_counts.to_series()
        return pd.concat([main_series, marker_counts])

    def validate_marker_counts(self):
        attrs = get_attrs(self.marker_counts)
        counts = [getattr(self.marker_counts, a) for a in attrs]
        if sum(counts) != self.markers:
            raise ValueError(
                f'Sum of marker counts ({sum(counts)}) is different from '
                f'the expected total marker count ({self.markers}).'
            )


@dataclass
class SampleStatistics:
    """Class for storing CheckM statistics of all the bins
        belonging to one sample"""
    sample_id: str
    bins: List[BinStatistics] = None

    def to_df(self) -> pd.DataFrame:
        all_series = [x.to_series() for x in self.bins]
        df = pd.concat(all_series, axis=1).T
        df['sample_id'] = self.sample_id

        # reorder columns
        cols_int = [
            'genomes', 'markers', 'marker_sets',
            'count0', 'count1', 'count2', 'count3', 'count4', 'count5_or_more'
        ]
        cols_float = ['completeness', 'contamination', 'strain_heterogeneity']
        cols_all = [
            'sample_id', 'bin_id', 'marker_lineage', 'marker_lineage_uid',
            *cols_int, *cols_float
        ]
        df = df[cols_all]

        # fix dtypes
        for col in cols_int:
            df[col] = df[col].astype(int)
        for col in cols_float:
            df[col] = df[col].astype(float)

        return df


def _process_checkm_arg(arg_key, arg_val):
    """Creates a list with argument and its value to be consumed by CheckM.

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
        return [f'--{arg_key}']
    else:
        return [f'--{arg_key}', str(arg_val)]


def _extract_pattern(
        pattern: str, text: str, desired_type: Type[Union[int, float]] = None
) -> list:
    result = re.findall(pattern, text)
    if desired_type:
        result = [desired_type(x.rstrip()) for x in result]
    return result


def _extract_checkm_stats(text: str) -> BinStatistics:
    words = _extract_pattern(r'(\S+)\s+(\D__\w+)\s\((UID\d+)\)', text)[0]
    ints = _extract_pattern(r'\s+(\d{1,})', text, int)[:10]
    floats = _extract_pattern(r'\d{1,3}\.\d{2}\s*', text, float)[:10]

    bin_stats = BinStatistics(
        bin_id=words[0],
        marker_lineage=words[1],
        marker_lineage_uid=words[2],
        genomes=ints[0],
        markers=ints[1],
        marker_sets=ints[2],
        marker_counts=MarkerCounts(
            count0=ints[3], count1=ints[4], count2=ints[5],
            count3=ints[6], count4=ints[7], count5_or_more=ints[8]
        ),
        completeness=floats[0],
        contamination=floats[1],
        strain_heterogeneity=floats[2]
    )
    bin_stats.validate_marker_counts()
    return bin_stats


def get_attrs(obj, excluded=()):
    return [k for k, v in vars(obj).items()
            if k not in excluded and not k.startswith('__')]
