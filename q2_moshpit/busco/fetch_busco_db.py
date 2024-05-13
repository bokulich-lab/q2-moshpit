# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from q2_moshpit._utils import colorify, run_command
from q2_moshpit.busco.types import BuscoDatabaseDirFmt


def fetch_busco_db(virus: bool, prok: bool, euk: bool) -> BuscoDatabaseDirFmt:
    busco_db = BuscoDatabaseDirFmt(path=None, mode='w')

    # Parse kwargs
    if all([virus, prok, euk]):
        args = ["all"]
    else:
        args = []
        variable_and_flag = [
            ('virus', virus), ('prokaryota', prok), ('eukaryota', euk)
        ]
        for variable_name, flag in variable_and_flag:
            if flag:
                args.append(variable_name)

    # Download
    print(colorify("Downloading BUSCO database..."))
    run_command(cmd=["busco", "--download", *args], cwd=str(busco_db))

    # Let user know that the process is compleat but it still needs
    # some time to copy files over.
    print(colorify(
        "Download completed. \n"
        "Copying files from temporary directory to final location..."
    ))

    return busco_db
