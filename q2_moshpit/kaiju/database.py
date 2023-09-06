# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import tarfile
from urllib.parse import urljoin

from tqdm import tqdm

from q2_types_genomics.kaiju import KaijuDBDirectoryFormat


from bs4 import BeautifulSoup
import requests
import pandas as pd

CHUNK_SIZE = 8192
KAIJU_SERVER_URL = "https://kaiju.binf.ku.dk/server"
ERR_MSG = (
    "Unable to connect to the Kaiju server. Please try again later. "
    "The error was: {}"
)


def _fetch_and_extract_db(db_uri: str, db_dir: str):
    latest_db = os.path.basename(db_uri)
    db_path = os.path.join(db_dir, latest_db)
    try:
        response = requests.get(db_uri, stream=True)
        response.raise_for_status()
        total_size = int(response.headers.get("content-length", 0))
        if total_size > 0:
            progress_bar = tqdm(
                desc=f'Downloading the "{latest_db}" database',
                total=total_size,
                unit="B",
                unit_scale=True,
                unit_divisor=1024,
            )

        with open(db_path, "wb") as file:
            for chunk in response.iter_content(chunk_size=CHUNK_SIZE):
                file.write(chunk) if chunk else False
                if total_size > 0:
                    progress_bar.update(len(chunk))
            progress_bar.close() if total_size > 0 else False
    except requests.exceptions.ConnectionError as e:
        raise Exception(ERR_MSG.format(e))

    msg = "Download finished. Extracting database files..."
    print(f"{msg}", end="", flush=True)
    with tarfile.open(db_path, "r:gz") as tar:
        tar.extractall(path=db_dir)
    print(f"\r{msg} Done.", flush=True)

    os.remove(db_path)


def _find_latest_db_url(database_type, sidebox_element, url):
    # Extract the databases and dates
    df = _find_all_dbs(sidebox_element)

    # Filter databases based on target_database type
    filtered_df = df[df.index.str.contains(database_type)]

    # Find the latest database
    latest_database = filtered_df["Date"].idxmax()
    # latest_database = filtered_df.loc[latest_index, "Database"]
    download_link = sidebox_element.find("a", string=latest_database)["href"]
    download_link = urljoin(url, download_link)

    return download_link


def _find_all_dbs(sidebox_element):
    databases, dates = [], []
    for link in sidebox_element.find_all("a"):
        database = link.get_text()
        date = database.split()[-2]  # Last element is the date
        databases.append(database)
        dates.append(date)
    df = pd.DataFrame({"Database": databases, "Date": dates})
    df.set_index("Database", inplace=True)
    df.loc[:, "Date"] = pd.to_datetime(df.loc[:, "Date"])
    return df


def fetch_kaiju_db(
    database_type: str,
) -> KaijuDBDirectoryFormat:

    try:
        response = requests.get(KAIJU_SERVER_URL)
    except requests.exceptions.RequestException as e:
        raise Exception(ERR_MSG.format(e))
    soup = BeautifulSoup(response.content, "html.parser")
    sidebox_db = soup.find("div", id="sidebox_db")

    download_link = _find_latest_db_url(
        database_type, sidebox_db, KAIJU_SERVER_URL
    )

    db = KaijuDBDirectoryFormat()
    _fetch_and_extract_db(download_link, str(db.path))

    return db
