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

import pandas as pd

from q2_moshpit._utils import run_command
from q2_moshpit.kraken2.select import kraken2_to_features
from q2_types_genomics.kraken2 import (
    Kraken2ReportDirectoryFormat,
    BrackenDBDirectoryFormat,
)


def _run_bracken_one_sample(
        bracken_db: str, kraken2_report_dir: str, bracken_report_dir: str,
        tmp_dir: str, threshold: int, read_len: int, level: str
) -> pd.DataFrame:
    sample_id = os.path.basename(kraken2_report_dir)
    kraken2_report_fp = os.path.join(
        kraken2_report_dir, f"{sample_id}.report.txt"
    )
    bracken_output_fp = os.path.join(
        tmp_dir, f"{sample_id}.bracken.output.txt"
    )
    sample_dir = os.path.join(bracken_report_dir, sample_id)
    os.makedirs(sample_dir, exist_ok=True)
    bracken_report_fp = os.path.join(
        sample_dir, "bracken.report.txt"
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
    bracken_table["taxonomy_id"] = bracken_table["taxonomy_id"].astype(str)

    return bracken_table


def _estimate_bracken(
        kraken_reports: Kraken2ReportDirectoryFormat,
        bracken_db: BrackenDBDirectoryFormat,
        threshold: int,
        read_len: int,
        level: str
) -> (pd.DataFrame, Kraken2ReportDirectoryFormat):
    bracken_tables = []
    bracken_reports = Kraken2ReportDirectoryFormat()

    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            for report_dir in kraken_reports.path.iterdir():
                bracken_table = _run_bracken_one_sample(
                    bracken_db=str(bracken_db),
                    kraken2_report_dir=report_dir,
                    bracken_report_dir=str(bracken_reports),
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
        index="sample_id", columns="taxonomy_id", values="new_est_reads"
    )
    bracken_table.fillna(0, inplace=True)

    return bracken_table, bracken_reports


def _assert_read_lens_available(
        bracken_db: BrackenDBDirectoryFormat, read_len: int
):
    pattern = r'.+database(\d{2,})mers\.kmer_distrib$'
    lengths = []
    for db in bracken_db.path.iterdir():
        lengths.extend(re.findall(pattern, str(db)))
    lengths = sorted([int(x) for x in lengths])
    if read_len not in lengths:
        raise ValueError(
            f"Provided read length ({read_len}) is not available in the "
            f"Bracken DB. The available values are: "
            f"{', '.join([str(x) for x in lengths])}."
        )


def estimate_bracken(
    kraken_reports: Kraken2ReportDirectoryFormat,
    bracken_db: BrackenDBDirectoryFormat,
    threshold: int = 0,
    read_len: int = 100,
    level: str = 'S'
) -> (Kraken2ReportDirectoryFormat, pd.DataFrame, pd.DataFrame):
    _assert_read_lens_available(bracken_db, read_len)

    table, reports = _estimate_bracken(
        kraken_reports=kraken_reports, bracken_db=bracken_db,
        threshold=threshold, read_len=read_len, level=level
    )

    _, taxonomy = kraken2_to_features(
        reports=reports, coverage_threshold=0.0
    )

    return reports, taxonomy, table
