# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import base64
import glob
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Union

import pandas as pd
from q2_types.per_sample_sequences import (
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)

from q2_moshpit._utils import run_command, _process_common_input_params
from q2_types_genomics.kaiju import KaijuDBDirectoryFormat

DEFAULT_PREFIXES = ["d__", "p__", "c__", "o__", "f__", "g__", "s__", "ssp__"]


def _get_sample_paths(df_index, df_row, paired):
    if paired:
        sample_name, fps = df_index, df_row.tolist()
    else:
        sample_name, fps = df_index, [df_row["forward"], ""]
    return sample_name, fps


def _rename_taxon(_id: str, id_to_taxon: dict) -> str:
    _id = id_to_taxon[_id].rstrip(";").split(";")
    x = zip(_id, DEFAULT_PREFIXES)
    return ";".join([f"{prefix}{taxon}" for taxon, prefix in x])


def _clean_terminal_ranks(x: str) -> str:
    taxa = x.split(";")
    while taxa[-1].endswith('__'):
        taxa = taxa[:-1]
    return ";".join(taxa)


def _encode_unclassified_ids(table: pd.DataFrame, text: str) -> pd.DataFrame:
    taxon = table.loc[table["taxon_name"].str.startswith(text), "taxon_name"].iloc[0]
    encoded = base64.b64encode(taxon.encode()).decode()
    table.loc[table["taxon_name"].str.startswith(text), "taxon_id"] = encoded[:8]
    return table


def _fix_id_types(table: pd.DataFrame) -> pd.DataFrame:
    table["taxon_id"].fillna(0, inplace=True)
    table['taxon_id'] = table['taxon_id'].astype(int)
    table = _encode_unclassified_ids(table, "cannot be assigned")
    table = _encode_unclassified_ids(table, "unclassified")
    return table


def _construct_feature_table(table_fp: str) -> (pd.DataFrame, pd.DataFrame):
    table = pd.read_csv(table_fp, sep="\t")

    # clean up taxon IDs
    table = _fix_id_types(table)

    # extract sample name from the file path
    table["sample"] = table["file"].map(lambda x: Path(x).stem)

    # clean up all the NAs
    table["taxon_name"] = table["taxon_name"].str.replace("NA", "Unspecified")

    # rename taxon IDs to taxon names and clean up
    taxa = table.set_index('taxon_id')['taxon_name'].to_dict()
    table["taxon_name"] = table["taxon_id"].map(lambda x: _rename_taxon(x, taxa))
    # table["taxon_name"] = table["taxon_name"].map(lambda x: _clean_terminal_ranks(x))

    # create taxonomy table
    taxonomy = table[["taxon_id", "taxon_name"]].drop_duplicates()
    taxonomy.set_index("taxon_id", inplace=True, drop=True)
    taxonomy.index = taxonomy.index.astype(str)
    taxonomy.index.name = "Feature ID"
    taxonomy.columns = ["Taxon"]

    # convert to sample x feature format
    table = table.groupby(["taxon_id", "sample"], as_index=False)["reads"].sum()
    table = table.pivot(index="sample", columns="taxon_id", values="reads")

    # convert column names to strings
    table.columns = table.columns.astype(str)

    return table, taxonomy


def _process_kaiju_reports(tmpdir, all_args):
    table_args = [
        "-r", all_args["r"],
        "-t", os.path.join(str(all_args["db"].path), "nodes.dmp"),
        "-n", os.path.join(str(all_args["db"].path), "names.dmp"),
        "-l", "superkingdom,phylum,class,order,family,genus,species",
    ]
    if all_args["exp"]:
        table_args.append("-e")
    if all_args["u"]:
        table_args.append("-u")
    if 0 <= all_args["c"] <= 1:
        table_args.extend(["-m", str(all_args["c"])])
    else:
        table_args.extend(["-c", str(all_args["c"])])

    report_fps = sorted(glob.glob(os.path.join(tmpdir, "*.out")))
    table_fp = os.path.join(tmpdir, "results.tsv")

    cmd = [
        "kaiju2table", "-v", "-o", table_fp, *table_args, *report_fps
    ]
    try:
        run_command(cmd=cmd, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running Kaiju2Table, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )

    return _construct_feature_table(table_fp)


def _classify_kaiju(manifest: pd.DataFrame, all_args: dict) -> (pd.DataFrame, pd.DataFrame):
    kaiju_args = [
        "-z", str(all_args["z"]),
        "-a", str(all_args["a"]),
        "-e", str(all_args["e"]),
        "-m", str(all_args["m"]),
        "-s", str(all_args["s"]),
        "-E", str(all_args["evalue"]),
        "-x" if all_args["x"] else "-X",
        "-t", os.path.join(str(all_args["db"].path), "nodes.dmp"),
        "-f", glob.glob(os.path.join(str(all_args["db"].path), "kaiju_*.fmi"))[0],
    ]

    base_cmd = ["kaiju-multi", "-v", *kaiju_args]
    paired = "reverse" in manifest.columns

    samples_fwd, samples_rev, output_fps = [], [], []
    with tempfile.TemporaryDirectory() as tmpdir:
        for index, row in manifest.iterrows():
            sample_name, fps = _get_sample_paths(index, row, paired)
            samples_fwd.append(fps[0])
            samples_rev.append(fps[1])
            output_fps.append(f"{os.path.join(tmpdir, sample_name)}.out")

        base_cmd.extend(["-i", ",".join(samples_fwd)])
        if paired:
            base_cmd.extend(["-j", ",".join(samples_rev)])
        base_cmd.extend(["-o", ",".join(output_fps)])

        try:
            run_command(cmd=base_cmd, verbose=True)
        except subprocess.CalledProcessError as e:
            raise Exception(
                "An error was encountered while running Kaiju, "
                f"(return code {e.returncode}), please inspect "
                "stdout and stderr to learn more."
            )

        table, taxonomy = _process_kaiju_reports(tmpdir, all_args)

    return table, taxonomy


def classify_kaiju(
    seqs: Union[
        SingleLanePerSamplePairedEndFastqDirFmt,
        SingleLanePerSampleSingleEndFastqDirFmt,
    ],
    db: KaijuDBDirectoryFormat,
    z: int = 1,
    a: str = "greedy",
    e: int = 3,
    m: int = 11,
    s: int = 65,
    evalue: float = 0.01,
    x: bool = True,
    r: str = "species",
    c: float = 0.0,
    exp: bool = False,
    u: bool = False,
) -> (pd.DataFrame, pd.DataFrame):
    manifest: pd.DataFrame = seqs.manifest.view(pd.DataFrame)
    return _classify_kaiju(manifest, dict(locals().items()))
