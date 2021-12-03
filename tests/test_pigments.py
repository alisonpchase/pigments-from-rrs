import pandas as pd
import pytest

from pigments_from_rrs.spectra import G1, G2, LNOT, load_data


@pytest.fixture(scope="session")
def temporary_directory(tmp_path_factory):
    return tmp_path_factory.mktemp("test-directory")


@pytest.fixture(scope="session")
def example_dataframe():
    data_dict = dict(
        first_column_name=[1, 2, 3],
        second_column_name=[4, 5, 6],
    )
    dataframe = pd.DataFrame.from_dict(data_dict)
    return dataframe


def test_constants():
    assert G1 == 0.0949
    assert G2 == 0.0794
    assert LNOT == 400


def test_load_data(example_dataframe, temporary_directory):
    # Define a filepath to write a test CSV to
    example_csv_filepath = f"{temporary_directory}/test.csv"
    # Write our test CSV to that path
    example_dataframe.to_csv(example_csv_filepath, index=False)
    # Use our load function to load the CSV
    dataframe_from_csv = load_data(path=example_csv_filepath)
    # Use pandas' built-in testing to ensure `load_data` worked
    pd.testing.assert_frame_equal(example_dataframe, dataframe_from_csv)
