# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


#from os import strerror
import subprocess

from q2_types_genomics.feature_data import DiamondDB, ArbitraryHeaderTSVDirFmt
from q2_types.feature_data import FeatureData
from qiime2.core.type.primitive import Str

def create_reference_db(mode: Str, target_taxa: Str, name: Str = None) -> FeatureData[DiamondDB]:

    # setup default name if none provided
    if not name:
        name = str(mode) + "_DB"

    # generate single taxa input
    #target_taxa = [select_taxa for select_taxa in [taxa, taxids] if
    #        select_taxa is not None]
    #print(target_taxa)

    # instantiate artifact
    db_to_return = ArbitraryHeaderTSVDirFmt()
    output_path = str(db_to_return)

    # parse/check taxa inputs....
    try:
        taxa_type, taxa_vals = _check_taxa(target_taxa)
    except Exception as e:
        raise e

    # setup call inputs for call to script
    cmd_str = r'create_dbs.py --dbname {} {} {}'.format(name,
            taxa_type, taxa_vals)
    cmds = cmd_str.split(" ")


    # actual download
    local_db = subprocess.run(cmds, check=True, stderr=subprocess.PIPE)


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
            taxid_vals = str({int(val.strip()) for val in split_taxa})
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
