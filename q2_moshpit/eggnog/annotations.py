# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import numpy as np
import pandas as pd

from q2_types.feature_data_mag import OrthologAnnotationDirFmt


def _extract_cog(data: pd.DataFrame) -> pd.DataFrame:
    data = data.set_index("seed_ortholog", inplace=False)
    data = data["COG_category"].fillna("")
    data = data.apply(lambda x: pd.Series(list(x)) if x is not np.nan else x)
    data = data.stack().reset_index(level=1, drop=True)
    data = data.reset_index(name='COG_category')
    data = data.pivot_table(
        index='seed_ortholog', columns='COG_category',
        aggfunc='size', fill_value=0
    )
    data = data.drop("-", axis=1, inplace=False)
    return data


def _merge(
        cogs: pd.DataFrame, frequencies: pd.DataFrame, mag_id: str
) -> pd.DataFrame:
    ids = cogs.index.to_list()
    result = frequencies.loc[ids, mag_id].dot(cogs)
    return result


def extract_annotations(
        ortholog_frequencies: pd.DataFrame,
        ortholog_annotations: OrthologAnnotationDirFmt,
        annotation: str
) -> pd.DataFrame:
    if annotation == "cog":
        extract_method = _extract_cog
    else:
        raise NotImplementedError(
            f"Annotation method {annotation} not supported"
        )

    cogs = []
    for _id, fp in ortholog_annotations.annotation_dict().items():
        annot_df = pd.read_csv(
            fp, sep="\t", skiprows=4, index_col=0
        )
        cog_df = extract_method(annot_df)
        cog_freqs = _merge(cog_df, ortholog_frequencies.T, _id)
        cogs.append(cog_freqs)

    result = pd.concat(cogs, axis=1).T
    result.index.name = "id"

    return result
