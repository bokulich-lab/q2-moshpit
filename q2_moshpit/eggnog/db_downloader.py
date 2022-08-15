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
from itertools import chain as chain
import pathlib
import shutil

from q2_types_genomics.feature_data import BinaryReferenceDBFmt


def create_reference_db(mode: str, target_taxa: str, name: str = None,
                        simulate: bool = False) -> BinaryReferenceDBFmt:

    # temp directory to download into & generate filepath
    test_dir = tempfile.TemporaryDirectory()
    test_dir_path = pathlib.Path(test_dir.name)
    os.environ["DATA_PATH"] = test_dir_path.as_posix()

    # arrange parameters
    if not name:
        name = "{}RefDB".format(mode)

    flags = ['-y']
    if simulate:
        flags.append('-s')

    taxa_type, taxa_vals = _check_taxa(target_taxa)

    cmds = list(chain(['create_dbs.py', '-m', mode, '--dbname', name,
                       taxa_type, taxa_vals, '--data_dir', test_dir_path],
                      flags))

    # filepath of download target
    download_log_fp = pathlib.Path(test_dir_path, "download_log")
    downloaded_db_fp = pathlib.Path(test_dir_path, "{}.dmnd".format(name))

    # do the actual downloading
    with open(download_log_fp, "w") as dl_fp_log:
        # Have sub-process write its own actual log file.
        subprocess.run(cmds, stdout=dl_fp_log, stderr=dl_fp_log)

    # instantiate format object
    download_db = BinaryReferenceDBFmt()

    # we need a return either way, but if not simulate want the downloaded
    # data actually written into the object

    if not simulate:
        shutil.copy(downloaded_db_fp, download_db.path / "eggnog.tsv")

    return download_db


def _check_taxa(target_taxa: list):
    if target_taxa[0].isnumeric():
        taxa_type = "--taxids"
        for taxon in target_taxa:
            if not taxon.isnumeric():
                raise ValueError("All taxa inputs must be the same type,"
                                 " either all taxids as integers or all"
                                 " string labels.")

    else:
        taxa_type = "--taxa"
        for taxon in target_taxa:
            if taxon.isnumeric():
                raise ValueError(
                    "All taxa inputs must be the same type,"
                    " either all taxids as integers or all string labels.")

    return taxa_type, ",".join(target_taxa)
