import pytest
from q2_types_genomics.per_sample_data._format import MultiMAGSequencesDirFmt
from zipfile import ZipFile
import os


@pytest.fixture(scope="session")
def MAGS_test_data(tmp_path_factory):
    """
    Fixture for BUSCO testing. It unzips a directoy with some test data
    and initializes a MultiMAGSequencesDirFmt which gets passed to the tests.

    Returns:
        test_data (MultiMAGSequencesDirFmt): object with the data for tests
    """

    # Get temporary directory
    tmp_path = tmp_path_factory.mktemp("busco_pytest_tmpdir")

    # Unzip files into temp dir
    with ZipFile("./data/test_data.zip", "r") as f:
        f.extractall(tmp_path)

    # NOTE: tmp_path_factory is a fixed name defined in pytest which is
    # passed on to your testfunction if you have it as an argument name.

    # Initialize object with path to unziped files
    test_data = MultiMAGSequencesDirFmt(
        path=os.path.join(tmp_path),
        mode="r",
    )

    # Return
    return test_data
