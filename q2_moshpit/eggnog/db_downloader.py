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
from qiime2.core.type.primitive import Bool, Str


def create_reference_db(mode: Str, target_taxa: Str, name: Str = None,
                        simulate: Bool = False) -> BinaryReferenceDBFmt:

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


def _check_taxa(taxa: Str) -> (str, str):
    try:
        split_taxa = [val.strip() for val in taxa.split(",")]
    except Exception as e:
        raise e
    #    raise ValueError("Unable to parse provided `Taxa` values, target"
    #                     " values should be provided as string with each"
    #                     " target taxa value seperated by a `,`")

    # if first value is an integer(isnumeric checks this because it looks for
    # only ascii values 0x30-0x39 only, so the "." in non-integer
    # numbers will return false as it finds 0x2E), then checks to make
    # sure all subsequent values are integers as well.
    if split_taxa[0].isnumeric():
        taxa_type = "--taxids"
        try:
            taxid_vals = ",".join({val.strip() for val in split_taxa})
            return (taxa_type, taxid_vals)

        except ValueError:
            raise ValueError("All taxa inputs must be the same type, either"
                             " all taxids as integers or all string labels.")
        except Exception as e:
            raise e

    else:
        taxa_type = "--taxa"
        try:
            for each in split_taxa:
                assert not each.isnumeric()
            return (taxa_type, ",".join(split_taxa))
        except ValueError:
            raise ValueError("All taxa inputs must be the same type, either"
                             " all taxids as integers or all string labels.")
