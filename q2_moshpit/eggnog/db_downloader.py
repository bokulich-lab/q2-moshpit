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

from q2_types_genomics.feature_data import DiamondDB, ArbitraryHeaderTSVDirFmt
from q2_types.feature_data import FeatureData
from qiime2.core.type.primitive import Bool, Str

def create_reference_db(mode: Str, target_taxa: Str, name: Str = None,
        simulate: Bool = False) -> FeatureData[DiamondDB]:

    #setup working directory
    test_dir = tempfile.TemporaryDirectory()
    test_dir_path = test_dir.name
    os.environ["DATA_PATH"] = test_dir_path

    # setup defaults
    if not name:
        name = str(mode) + "_DB"

    flags = ['-y']
    if simulate:
        flags.append('-s')

    # instantiate artifact
    db_to_return = ArbitraryHeaderTSVDirFmt()
    output_path = str(db_to_return)
    print(output_path)

    # parse/check taxa inputs....
    try:
        taxa_type, taxa_vals = _check_taxa(target_taxa)
    except Exception as e:
        raise e
    print("Taxa type: {}, taxa_vals: {}".format(taxa_type, taxa_vals))
    print("flags: {}".format(flags))
    # setup call inputs for call to script
    #cmd_str = r'create_dbs.py --dbname {} {} {} -y'.format(name,
    #        taxa_type, taxa_vals)
    #cmds = cmd_str.split(" ")
    cmds = list(chain(['create_dbs.py', '--dbname', name, taxa_type, taxa_vals,
            '--data_dir', test_dir_path], flags))
    print("commands chain: {}".format(cmds))
    # actual download
    local_db = subprocess.run(cmds, capture_output=True)
            #stderr=subprocess.STDOUT)
    # we could use the following regex to get the path to the downloaded data
    # file: re.search("(?<=\-\-db)\s.*\s", <output-from our subprocess run>
    # or better yet because we know the directory and the name of the database
    # file and we could create a path from that information to the file that
    # we can use to instantiate an artifact object to to feed into our Format
    # object to creat our Artifact.


    # framework saves as correct dirfmt up return...
    return local_db

def _check_taxa(taxa: Str) -> (str, str):
    try:
        split_taxa = [val.strip() for val in taxa.split(",")]
    except:
        raise ValueError("Unable to parse provided `Taxa` values, target"
                         " values should be provided as string with each"
                         " target taxa value seperated by a `,`")

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
        #except ValueError:
        #    raise TypeError("The supplied taxa value: {} cannot be recognized"
        #                    " as a valid taxa input. Make sure you input a"
        #                    " list of taxa strings or a list of integers as"
        #                    " taxa ids")
