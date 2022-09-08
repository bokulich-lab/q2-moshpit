# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import subprocess
import tempfile
from itertools import chain as chain
import pathlib
import shutil

from ._utils import _check_taxa
from q2_types_genomics.eggnog import EggnogRefDirFmt, ReferenceDB, Eggnog


def get_references(mode: str = None, target_taxa: str = None, name: str =
        None, simulate: bool = False) -> ReferenceDB[Eggnog]:
    pass

def create_reference_db(mode: str, target_taxa: str = None, name: str = None,
                        simulate: bool = False) -> ReferenceDB[Eggnog]:

    # temp directory to download into & generate filepath
    test_dir = tempfile.TemporaryDirectory()
    test_dir_path = pathlib.Path(test_dir.name)

    # arrange parameters
    allowed_modes = {'diamond': 'dmnd', 'mmseqs': 'mmseqs', }
    if mode in allowed_modes.keys():
        file_ext = allowed_modes[mode]
    else:
        raise ValueError("Please supply a valid mode from the following"
                         " choices: {}".format(allowed_modes.keys()))

    if not name:
        name = "{}RefDB".format(mode)

    cmds = ['create_dbs.py', '-m', mode, '--dbname', name,
            taxa_type, taxa_vals, '--data_dir', test_dir_path,
            ]

    saved_file_name = "{}.{}".format(name, file_ext)

    flags = ['-y']
    if simulate:
        flags.extend('-s')

    cmds.extend(flags)

    cmds.extend(list(_check_taxa(target_taxa)))

    # filepath of download target
    downloaded_db_fp = pathlib.Path(test_dir_path, saved_file_name)
    print(downloaded_db_fp)

    # raise Exception("temp directory: {}\n-----\nfilename:"
    #     " {}\n-----\ndownloaded db: {}".format(test_dir_path,
    #             saved_file_name, downloaded_db_fp))


    # do the actual downloading
    subprocess.run(cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # instantiate format object
    download_db = EggnogRefDirFmt()

    # we need a return either way, but if not simulate want the downloaded
    # data actually written into the object

    if not simulate:
        shutil.copy(downloaded_db_fp, download_db.path)

    return download_db


def download_references(target_taxa: str = None,
                        simulate: bool = False) -> ReferenceDB[Eggnog]:
    working_dir = tempfile.TemporaryDirectory()

    # setup download commands
    cmds = ["download_eggnog_data.py", "--data_dir", working_dir.name, "-P",
            "-M", "-F", "-y",
            ]
    if simulate:
        cmds.extend(["-s"])

    if target_taxa is not None:
        taxa_type, taxa_vals = _check_taxa(target_taxa)
        cmds.extend([taxa_type, taxa_vals])

    # do the actual downloading
    subprocess.run(cmds, check=True)

    downloaded_things = sorted(pathlib.Path(working_dir.name).glob("*"))
    print(downloaded_things)

    # # artifact to store this downloaded data.
    # reference_artifact = EggnogRefDirFmt()
    # # now copy the downloaded data into the artifact.

    # return reference_artifact
