import pytest
from q2_types_genomics.per_sample_data._format import MultiMAGSequencesDirFmt


@pytest.fixture(scope="session")
def image_file(tmpdir_factory):
    img = compute_expensive_image()
    fn = tmpdir_factory.mktemp("data").join("img.png")
    img.save(str(fn))
    return fn


@pytest.fixture(scope="session")
def MultiMAGSequencesDirFmt_test_data():
    """
    Fixture for BUSCO testing. Downloads directory with MAGs for three samples and stores it in a temoprary directory.

    Returns:
        test_data (MultiMAGSequencesDirFmt): object with the data for tests
    """

    # TODO: make temp dir for test data

    # TODO: Download data from github into temp dir

    # Initialize object with downloaded data
    test_data = MultiMAGSequencesDirFmt(
        path="/Users/santiago/moshpit_tutorial/0a2e152b-c5bb-4f61-ad75-d31d7a9b5ab6",
        mode="r",
    )

    return test_data
