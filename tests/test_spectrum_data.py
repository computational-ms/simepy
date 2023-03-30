from pathlib import Path
from spectrum_data.extract_scans import extract_scan_data
from spectrum_data.extract_meta_data import extract_meta_data

from tempfile import NamedTemporaryFile
import pandas as pd
import pytest


def test_extract_scans():
    input_file = Path(__file__).parent / "data" / "BSA1.mzML"

    with NamedTemporaryFile() as tmp_file:
        extract_scan_data(input_file, tmp_file.name)
        df = pd.read_csv(tmp_file.name)
    assert df["filename"].nunique() == 1
    assert df["filename"].iloc[0] == "BSA1.mzML"
    assert df["spectrum_id"].nunique() == 1684
    assert pytest.approx(df["mz"].mean(), 449.397186848143)


def test_extract_meta_data_mzml():
    input_file = Path(__file__).parent / "data" / "BSA1.mzML"

    with NamedTemporaryFile() as scans, NamedTemporaryFile() as run, NamedTemporaryFile() as noise, NamedTemporaryFile() as instrument:
        extract_meta_data(
            input_file=input_file,
            run_output_file=run.name,
            spec_output_file=scans.name,
            noise_output_file=noise.name,
            unit_output_file=instrument.name,
            # object_name=args.object_name,
            # lineage_root=args.lineage_root,
        )

        scans = pd.read_csv(scans.name)
        run = pd.read_csv(run.name)
        noise = pd.read_csv(noise.name)
        instrument = pd.read_csv(instrument.name)

    assert scans["spectrum_id"].nunique() == 1684
    assert pytest.approx(scans["rt"].mean(), 2037.37997769969)

    assert len(run) == 1
    assert run["spectrum_count"].iloc[0] == 1684

    assert len(noise) == 0
    assert len(instrument) == 0


def test_extract_meta_data_raw():
    pass
