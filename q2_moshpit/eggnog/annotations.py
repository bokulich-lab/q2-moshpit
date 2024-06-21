# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pandas as pd

from q2_types.feature_data_mag import OrthologAnnotationDirFmt


def _extract_generic(
        data: pd.DataFrame, column: str, transform_func: callable
) -> pd.Series:
    data = data.set_index("seed_ortholog", inplace=False)[column]
    data = data.apply(transform_func if data.notna().any() else lambda x: x)
    data = data.stack().reset_index(level=1, drop=True)
    data = data.value_counts()
    for char in ("", "-"):
        if char in data.index:
            data = data.drop(char, axis=0, inplace=False)
    return data


def _extract_cog(data: pd.DataFrame) -> pd.Series:
    return _extract_generic(
        data, "COG_category", lambda x: pd.Series(list(x))
    )


def _extract_kegg_ko(data: pd.DataFrame) -> pd.Series:
    return _extract_generic(
        data, "KEGG_ko", lambda x: pd.Series([i[3:] for i in x.split(",")])
    )


def _extract_kegg_pathway(data: pd.DataFrame) -> pd.Series:
    return _extract_generic(
        data, "KEGG_Pathway", lambda x: pd.Series(
            [i for i in x.split(",") if i.startswith("map")]
        )
    )


def _extract_kegg_module(data: pd.DataFrame) -> pd.Series:
    return _extract_generic(
        data, "KEGG_Module", lambda x: pd.Series(x.split(","))
    )


def _extract_kegg_reaction(data: pd.DataFrame) -> pd.Series:
    return _extract_generic(
        data, "KEGG_Reaction", lambda x: pd.Series(x.split(","))
    )


def _extract_brite(data: pd.DataFrame) -> pd.Series:
    return _extract_generic(
        data, "BRITE", lambda x: pd.Series(x.split(","))
    )


def _extract_caz(data: pd.DataFrame) -> pd.Series:
    return _extract_generic(
        data, "CAZy", lambda x: pd.Series(x.split(","))
    )


def _filter(
        data: pd.DataFrame, max_evalue: float, min_score: float
) -> pd.DataFrame:
    data = data[(data["evalue"] <= max_evalue) & (data["score"] >= min_score)]
    if len(data) == 0:
        raise ValueError(
            "E-value/score filtering resulted in an empty table - "
            "please adjust your thresholds and try again."
        )
    return data


def extract_annotations(
        ortholog_annotations: OrthologAnnotationDirFmt,
        annotation: str,
        max_evalue: float = 1.0,
        min_score: float = 0.0
) -> pd.DataFrame:
    extract_method = globals().get(f"_extract_{annotation}")
    if not extract_method:
        raise NotImplementedError(
            f"Annotation method {annotation} not supported"
        )

    annotations = []
    for _id, fp in ortholog_annotations.annotation_dict().items():
        annot_df = pd.read_csv(
            fp, sep="\t", skiprows=4, index_col=0
        )  # skip the first 4 rows as they contain comments
        annot_df = annot_df.iloc[:-3, :]  # remove the last 3 comment rows
        annot_df = _filter(annot_df, max_evalue, min_score)
        annot_df = extract_method(annot_df)
        annot_df.name = _id
        annotations.append(annot_df)

    result = pd.concat(annotations, axis=1).fillna(0).T
    result.index.name = "id"
    return result
