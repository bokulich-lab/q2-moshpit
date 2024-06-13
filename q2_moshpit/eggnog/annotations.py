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
        index='COG_category', columns='seed_ortholog',
        aggfunc='size', fill_value=0
    )
    data = data.drop("-", axis=0, inplace=False)
    return data


def _ensure_dim(table1: pd.DataFrame, table2: pd.DataFrame):
    # check that tables are compatible for dot product
    if table1.shape[1] != table2.shape[0]:
        raise ValueError(
            f"Tables do not have compatible dimensions for dot product: "
            f"{table1.shape[1]} != {table2.shape[0]}"
        )


def _merge(
        cogs: pd.DataFrame, frequencies: pd.DataFrame, mag_id: str
) -> pd.DataFrame:
    frequencies_filtered = frequencies.loc[mag_id, cogs.columns]
    _ensure_dim(cogs, frequencies_filtered)
    result = cogs.dot(frequencies_filtered)
    result.name = mag_id
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

    annotations = []
    for _id, fp in ortholog_annotations.annotation_dict().items():
        # we need to skip the first 4 rows as they contain comments
        annot_df = pd.read_csv(
            fp, sep="\t", skiprows=4, index_col=0
        )
        annot_df = extract_method(annot_df)
        annot_freqs = _merge(annot_df, ortholog_frequencies, _id)
        annotations.append(annot_freqs)

    result = pd.concat(annotations, axis=1).T
    result.index.name = "id"
    return result


def collapse_tables(
        mags_pa: pd.DataFrame,
        annotation_frequency: pd.DataFrame
) -> pd.DataFrame:
    _ensure_dim(mags_pa, annotation_frequency)
    result = mags_pa.dot(annotation_frequency)
    return result
