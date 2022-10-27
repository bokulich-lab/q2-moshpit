# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import tempfile
import os
import subprocess
import shutil
import glob

from q2_types_genomics.eggnog import  EggnogRefDirFmt, DiamondRefDirFmt

from q2_types_genomics.per_sample_data import ContigSequencesDirFmt
from q2_types.feature_data import DNAFASTAFormat, FeatureData

from q2_types.sample_data import SampleData
from q2_types_genomics.per_sample_data import Contigs

# inputs DNAFASTAFormat data, eggnog reference database, diamond reference database
# outputs: seed orthologs(fmt: ArbitraryHeaderTSVDirFmt, type: Ortholog[Seed])


# diamond search
def diamond_search(input_sequences: ContigSequencesDirFmt,
                   diamond_db:DiamondRefDirFmt,
                   eggnog_db:EggnogRefDirFmt,
                   )-> OrthologDirFmt:

    # TODO: wire up where to write results....
    #tmp = tempfile.TemporaryDirectory()

    # Instantiate an directory format object to which our processed data can be written
    # temp_result = tempfile.TemporaryDirectory()

    result = OrthologDirFmt()
    diamond_db_fp = diamond_db.path / 'eggnog_proteins.dmnd'

    for relpath, obj_path in input_sequences.sequences.iter_views(DNAFASTAFormat):
        sample_label = str(relpath).rsplit(r'\.', 2)[0] #+ '_seed_ortholog'
        print(sample_label)


        _diamond_search_runner(input_path=obj_path,
                               eggnog_db=eggnog_db,
                               diamond_db=diamond_db_fp,
                               sample_label=sample_label,
                               output=str(result))

    #if not os.path.isdir(os.path.join(result.path, 'data')):
    #    os.mkdir(os.path.join(result.path, 'data'))

    #for dirpath, dirname, filenames in os.walk(temp_result.name):
    #for entry in os.scandir(temp_result.name):
    #    if re.match("{}\.emapper\..*".format(sample_label), entry.name):
    #           shutil.copy(os.path.join(os.getcwd(), entry.name), os.path.join(result.path, 'data/'))
    #       else:
    #           print("sorry")
    #temp_files = glob.glob("%s.*".format(sample_label), root_dir=temp_result.name)
    #for file in temp_files:
    #    shutil.copy(os.path.join(temp_result.name, file), result.path)
    return result


#@def eggnog_mmseq2_search(input_sequences: ContigSequencesDirFmt,
#@                         mmseq_db: MMseqsDirFmt,
#@                         eggnog_db:EggnogRefDirFmt,
#@                         )-> OrthologDirFmt:
#@
#@    # TODO: wire up where to write results....
#@
#@    # Instantiate an directory format object to which our processed data can be written
#@    result = OrthologDirFmt()
#@
#@    for relpath, obj_path in input_sequences.sequences.iter_views(DNAFASTAFormat):
#@        sample_label = str(relpath).rsplit('_', 1)[0] #+ '_seed_ortholog'
#@        _mmseqs_search_runner(input_path=obj_path,
#@                              eggnog_db=eggnog_db,
#@                              mmseqs_db=mmseq_db,
#@                              sample_label=sample_label,
#@                              output=result)
#@        return result
#@
#@
# apply emapper to each sample sequences.
# this could be made more abstract to work for any mapping that we would like to perform with the emapper.py script or made more specifict to generate a helper method for each main method we would like to use emapper.py to process our data with. To start with I am making a more specific verision for just search diamond, and leave further abstraction for a later point.
def _diamond_search_runner(input_path, eggnog_db, diamond_db, sample_label, output):
    cmds = ['emapper.py',
            '--data_dir', str(eggnog_db),
            '--dmnd_db', str(diamond_db),
            '-i', str(input_path),
            '--itype', 'metagenome',
            '-o', sample_label,
            '--output_dir', output,
            '--no_annot',
            ]
    subprocess.run(cmds, check=True)
    return None

def _mmseqs_search_runner(input_path, eggnog_db, mmseqs_db, sample_label, output):
    cmds = ['emapper.py',
            '-m', 'mmseqs',
            '--data_dir', str(eggnog_db),
            '--mmseqs_db', str(mmseqs_db),
            '-i', str(input_path),
            '--itype', 'metagenome',
            '-o', sample_label,
            '--output_dir', str(output),
            '--no_annot',
            ]
    subprocess.run(cmds, check=True)
    return None

def _eggnog_annotate_runner(input_path, eggnog_db, output):
    cmds = [
            ]
    subprocess.run(cmds, check=True)
    return None

    # no return because the emapper.py script will write the results of subprocess.run to the location specied by output_target, which is pointed at our directory format object, which can be thought of as our real place on disk to write this to.

    #search_results = ArbitraryHeaderTSVFmt(mode='w')


    #print(search_results.path)
    #print(cmds)
    #assert False
    #subprocess.run(cmds, check=True)
    #return search_results



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
#    ['emapper.py','-m', 'mmseqs','--data_dir', str(eggnog_db),'--mmseqs_db', str(mmseqs_db),'-i', str(input_path),'--itype', 'metagenome','-o', sample_label, '--output_dir', str(output), '--no_annot']
