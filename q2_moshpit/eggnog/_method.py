# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import tempfile

import subprocess

from q2_types_genomics.eggnog import Ortholog, Seed, ReferenceDB, Eggnog, Diamond, BinaryReferenceDBDirFmt, OrthologDirFmt, ArbitraryHeaderTSVFmt
from q2_types.feature_data import DNAFASTAFormat, FeatureData

from q2_types.sample_data import SampleData
from q2_types_genomics.per_sample_data import Contigs

# diamond search
# inputs DNAFASTAFormat data, eggnog reference database, diamond reference database
# outputs: seed orthologs(fmt: ArbitraryHeaderTSVDirFmt, type: Ortholog[Seed])

def diamond_search(input_sequences: DNAFASTAFormat,
                   diamond_db: BinaryReferenceDBDirFmt,
                   eggnog_db: BinaryReferenceDBDirFmt)->  ArbitraryHeaderTSVFmt:
    result_name = "turkey"

    search_results = ArbitraryHeaderTSVFmt(mode='w')

    cmds = ['emapper.py', "--data_dir", str(eggnog_db), '--dmnd_db', str(diamond_db), '-i',
            str(query_sequences), '-o', result_name, "--output_dir",
            search_results.name,
            ]

    print(search_results.path)
    print(cmds)
    assert False
    subprocess.run(cmds, check=True)
    return search_results



#def e_mapper(main_db: BinaryReferenceDBDirFmt,
#             query_sequences: DNAFASTAFormat,
#             ancillary_db: BinaryReferenceDBDirFmt = None,
#             itype: str = None,
#             result_name: str = "egg_output",
#             mode: str = None,
#             ) -> ArbitraryHeaderTSVDirFmt:
#
#    working_dir = tempfile.TemporaryDirectory()
#    working_dir_path = pathlib.Path(working_dir.name)
#
#
#    if mode is not None:
#        cmds.extend(["-m", mode])
#
#    if itype is not None:
#        cmds.extend(["--itype", itype])
#
#
#    # annotation_fps =
#    print("*" * 13, sorted(working_dir_path.glob('*')), "*" * 13)
#    assert False
#
