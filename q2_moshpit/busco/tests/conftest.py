import pytest
from q2_types_genomics.per_sample_data._format import MultiMAGSequencesDirFmt


@pytest.fixture(scope="session")
def MAGS_test_data():
    """
    Fixture for BUSCO testing. It unzips a directoy with some test data and initializes a
    MultiMAGSequencesDirFmt which gets passed to the tests.

    Returns:
        test_data (MultiMAGSequencesDirFmt): object with the data for tests
    """

    # TODO: make temp dir for test data

    # TODO: expand zip into temp dir

    # TODO: Initialize object with downloaded data
    test_data = MultiMAGSequencesDirFmt(
        path="/Users/santiago/moshpit_tutorial/0a2e152b-c5bb-4f61-ad75-d31d7a9b5ab6",
        mode="r",
    )

    # Return
    return test_data
