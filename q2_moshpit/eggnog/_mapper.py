# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import pathlib
import tempfile

import subprocess

from q2_types.feature_data import DNAFASTAFormat
from q2_types_genomics.feature_data import (
    BinaryReferenceDBDirFmt, ArbitraryHeaderTSVDirFmt
    )


def e_mapper(main_db: BinaryReferenceDBDirFmt,
             query_sequences: DNAFASTAFormat,
             ancillary_db: BinaryReferenceDBDirFmt = None,
             itype: str = None,
             result_name: str = "egg_output",
             mode: str = None,
             ) -> ArbitraryHeaderTSVDirFmt:

    working_dir = tempfile.TemporaryDirectory()
    working_dir_path = pathlib.Path(working_dir.name)

    cmds = ['emapper.py', "--data_dir", str(main_db), '-i',
            str(query_sequences), '-o', result_name, "--output_dir",
            working_dir.name,
            ]

    if mode is not None:
        cmds.extend(["-m", mode])

    if itype is not None:
        cmds.extend(["--itype", itype])

    subprocess.run(cmds, check=True)

    # annotation_fps =
    print("*" * 13, sorted(working_dir_path.glob('*')), "*" * 13)
    assert False

    # anno_out = ArbitraryHeaderTSVDirFmt()
    # shutil.copy(src=annotation_fps[0], dst=anno_out.path)
    # return performed_annotations
    # return anno_out
