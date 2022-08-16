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

from ._utils import _check_taxa
from q2_types_genomics.feature_data import BinaryReferenceDBFmt


def create_reference_db(mode: str, target_taxa: str, name: str = None,
                        simulate: bool = False) -> BinaryReferenceDBFmt:

    # temp directory to download into & generate filepath
    test_dir = tempfile.TemporaryDirectory()
    test_dir_path = pathlib.Path(test_dir.name)
    os.environ["DATA_PATH"] = test_dir_path.as_posix()

    # arrange parameters
    if mode == 'diamond':
        file_ext = 'dmnd'
    elif mode == 'mmseqs':
        file_ext = 'mmseqs'
    else:
        raise ValueError("Please supply a valid mode")

    if not name:
        name = "{}RefDB".format(mode)

    saved_file_name = "{}.{}".format(name, file_ext)
    print(saved_file_name)

    flags = ['-y']
    if simulate:
        flags.append('-s')

    taxa_type, taxa_vals = _check_taxa(target_taxa)

    # filepath of download target
    downloaded_db_fp = pathlib.Path(test_dir_path, saved_file_name)
    print(downloaded_db_fp)

    # raise Exception("temp directory: {}\n-----\nfilename:"
    #     " {}\n-----\ndownloaded db: {}".format(test_dir_path,
    #             saved_file_name, downloaded_db_fp))

    cmds = list(chain(['create_dbs.py', '-m', mode, '--dbname', name,
                       taxa_type, taxa_vals, '--data_dir', test_dir_path],
                      flags))

    # do the actual downloading
    subprocess.run(cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # instantiate format object
    download_db = BinaryReferenceDBFmt()

    # we need a return either way, but if not simulate want the downloaded
    # data actually written into the object

    if not simulate:
        shutil.copy(downloaded_db_fp, download_db.path / "reference_database")

    return download_db
