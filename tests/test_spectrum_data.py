from pathlib import Path

import pytest

from simepy.extract_meta_data import extract_meta_data
from simepy.extract_meta_data import get_file_type
from simepy.extract_scans import extract_scan_data


def test_extract_scans():
    input_file = Path(__file__).parent / "data" / "BSA1.mzML"
    df = extract_scan_data(input_file)
    assert df["filename"].nunique() == 1
    assert df["filename"].iloc[0] == "BSA1.mzML"
    assert df["spectrum_id"].nunique() == 1684
    assert pytest.approx(df["mz"].mean()) == 449.397186848143


def test_extract_meta_data_mzml():
    input_file = Path(__file__).parent / "data" / "BSA1.mzML"

    run, instrument, scans, noise = extract_meta_data(
        input_file=input_file,
    )

    assert scans["spectrum_id"].nunique() == 1684
    assert pytest.approx(scans["rt"].mean()) == 2037.37997769969

    assert len(run) == 1
    assert run["spectrum_count"].iloc[0] == 1684

    assert len(noise) == 0
    assert len(instrument) == 0


def test_extract_meta_data_raw():
    pass


@pytest.mark.parametrize(
    "path,file_type",
    [
        (Path("some_file.thermorawfile.mzml"), "mzml"),
        (Path("another_file.mzml.gz"), "mzml.gz"),
        (Path("one_more_file.raw"), "raw"),
    ],
)
def test_get_file_type(path, file_type):
    assert get_file_type(path) == file_type
