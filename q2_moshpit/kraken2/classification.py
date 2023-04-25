# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import re
import subprocess
import tempfile
from copy import deepcopy
from typing import Union, List

import pandas as pd
from q2_types.per_sample_sequences import (
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt
)

from q2_moshpit._utils import run_command, _process_common_input_params
from q2_moshpit.kraken2.utils import _process_kraken2_arg
from q2_types_genomics.kraken2 import (
    Kraken2ReportDirectoryFormat,
    Kraken2OutputDirectoryFormat,
    Kraken2DBDirectoryFormat, BrackenDBDirectoryFormat
)
from q2_types_genomics.per_sample_data import MultiMAGSequencesDirFmt


def _get_seq_paths(df_index, df_row, df_columns):
    if "filename" in df_columns:
        _sample, _bin, fn = df_index[0], df_index[1], [df_row["filename"]]
    elif "reverse" in df_columns:
        _sample, _bin, fn = df_index, df_index, df_row.tolist()
    else:
        _sample, _bin, fn = df_index, df_index, [df_row["forward"]]
    return _sample, _bin, fn


def _construct_output_paths(
    _sample, _bin, kraken2_outputs_dir, kraken2_reports_dir
):
    sample_dir_report = os.path.join(kraken2_reports_dir.path, _sample)
    sample_dir_output = os.path.join(kraken2_outputs_dir.path, _sample)
    for s in [sample_dir_report, sample_dir_output]:
        os.makedirs(s, exist_ok=True)
    report_fp = os.path.join(sample_dir_report, f"{_bin}.report.txt")
    output_fp = os.path.join(sample_dir_output, f"{_bin}.output.txt")
    return output_fp, report_fp


def _classify_kraken(
    manifest, common_args
) -> (Kraken2ReportDirectoryFormat, Kraken2OutputDirectoryFormat):
    base_cmd = ["kraken2", *common_args]
    base_cmd.append("--paired") if "reverse" in manifest.columns else False

    kraken2_reports_dir = Kraken2ReportDirectoryFormat()
    kraken2_outputs_dir = Kraken2OutputDirectoryFormat()

    try:
        for index, row in manifest.iterrows():
            _sample, _bin, fn = _get_seq_paths(index, row, manifest.columns)
            output_fp, report_fp = _construct_output_paths(
                _sample, _bin, kraken2_outputs_dir, kraken2_reports_dir
            )
            cmd = deepcopy(base_cmd)
            cmd.extend(
                ["--report", report_fp, "--output", output_fp,
                 "--use-names", *fn]
            )
            run_command(cmd=cmd, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running Kraken 2, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )

    return kraken2_reports_dir, kraken2_outputs_dir


def _run_bracken_one_sample(
        bracken_db: str, kraken2_report_dir: str, tmp_dir: str,
        threshold: int, read_len: int, level: str
) -> pd.DataFrame:
    sample_id = os.path.basename(kraken2_report_dir)
    kraken2_report_fp = os.path.join(
        kraken2_report_dir, f"{sample_id}.report.txt"
    )
    bracken_output_fp = os.path.join(
        tmp_dir, f"{sample_id}.bracken.output.txt"
    )
    bracken_report_fp = os.path.join(
        tmp_dir, f"{sample_id}.bracken.report.txt"
    )
    cmd = [
        "bracken",
        "-d", bracken_db,
        "-i", kraken2_report_fp,
        "-o", bracken_output_fp,
        "-w", bracken_report_fp,
        "-t", str(threshold),
        "-r", str(read_len),
        "-l", level,
    ]
    try:
        run_command(cmd=cmd, verbose=True)
    except subprocess.CalledProcessError as e:
        # TODO: what should be the behaviour when no reads could be classified?
        raise Exception(
            "An error was encountered while running Bracken, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )
    bracken_table = pd.read_csv(bracken_output_fp, sep="\t", index_col=0)
    bracken_table["sample_id"] = sample_id

    return bracken_table


def _estimate_bracken(
        kraken_reports: Kraken2ReportDirectoryFormat,
        bracken_db: BrackenDBDirectoryFormat,
        threshold: int,
        read_len: int,
        level: str
) -> pd.DataFrame:
    bracken_tables = []
    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            for report_dir in kraken_reports.path.iterdir():
                bracken_table = _run_bracken_one_sample(
                    bracken_db=str(bracken_db),
                    kraken2_report_dir=report_dir,
                    tmp_dir=tmpdir, threshold=threshold,
                    read_len=read_len, level=level
                )
                bracken_tables.append(bracken_table)
        except subprocess.CalledProcessError as e:
            raise Exception(
                "An error was encountered while running Bracken, "
                f"(return code {e.returncode}), please inspect "
                "stdout and stderr to learn more."
            )

    bracken_table = pd.concat(bracken_tables).reset_index()
    bracken_table = bracken_table.pivot(
        index="sample_id", columns="name", values="new_est_reads"
    )

    return bracken_table


def _get_available_lens(bracken_db: BrackenDBDirectoryFormat) -> List[int]:
    pattern = r'.+database(\d{2,})mers\.kmer_distrib$'
    lengths = []
    for db in bracken_db.path.iterdir():
        lengths.extend(re.findall(pattern, str(db)))
    return sorted([int(x) for x in lengths])


def _is_reads(seqs):
    return isinstance(
            seqs, (SingleLanePerSamplePairedEndFastqDirFmt,
                   SingleLanePerSampleSingleEndFastqDirFmt)
    )


def classify_kraken(
    seqs: Union[
        SingleLanePerSamplePairedEndFastqDirFmt,
        SingleLanePerSampleSingleEndFastqDirFmt,
        MultiMAGSequencesDirFmt,
    ],
    kraken2_db: Kraken2DBDirectoryFormat,
    bracken_db: BrackenDBDirectoryFormat,
    threads: int = 1,
    confidence: float = 0.0,
    minimum_base_quality: int = 0,
    memory_mapping: bool = False,
    minimum_hit_groups: int = 2,
    quick: bool = False,
    threshold: int = 0,
    read_len: int = 100,
    level: str = 'S'
) -> (
        Kraken2ReportDirectoryFormat,
        Kraken2OutputDirectoryFormat,
        pd.DataFrame
):
    kraken2_args = [
        "threads", "confidence", "minimum_base_quality",
        "memory_mapping", "minimum_hit_groups", "quick"
    ]
    kwargs = {k: v for k, v in locals().items() if k in kraken2_args}
    common_args = _process_common_input_params(
        processing_func=_process_kraken2_arg, params=kwargs
    )
    common_args.extend(["--db", str(kraken2_db.path)])
    manifest: pd.DataFrame = seqs.manifest.view(pd.DataFrame)

    if _is_reads(seqs):
        # Check if read length is available in the Bracken DB
        available_lens = _get_available_lens(bracken_db)
        if read_len not in available_lens:
            raise ValueError(
                f"Provided read length ({read_len}) is not available in the "
                f"Bracken DB. The available values are: "
                f"{', '.join([str(x) for x in available_lens])}."
            )

    k2_reports, k2_outputs = _classify_kraken(manifest, common_args)

    # Process with Bracken
    # TODO: should we try to estimate read length from the input seqs
    #  or let user decide?
    if _is_reads(seqs):
        bracken_table = _estimate_bracken(
            kraken_reports=k2_reports, bracken_db=bracken_db,
            threshold=threshold, read_len=read_len, level=level
        )
    else:
        # TODO: for now, if MAGs were provided, output an empty table
        bracken_table = pd.DataFrame()
        bracken_table.index = [f'{x[0]}_{x[1]}' for x in manifest.index]

    return k2_reports, k2_outputs, bracken_table
