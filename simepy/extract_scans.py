"""Resource to extract scan info from a mzml file."""
import logging
from pathlib import Path

import pandas as pd
import pymzml


def extract_scan_data(input_file):
    """
    Extract scan info like mass, intensity, id and filename from an input mzml file.

    Args:
        input_file (PosixPath): path to input file

    Returns:
        scans_out (pandas.DataFrame): Dataframe containing scans information
    """
    logging.info("Starting scan info extraction!")
    with pymzml.run.Reader(str(input_file)) as reader:
        scans = []
        for spectrum in reader:
            scans_df = pd.DataFrame()
            scans_df["mz"] = spectrum.get_array("m/z array")
            scans_df["intensity"] = spectrum.get_array("intensity array")
            scans_df["spectrum_id"] = spectrum.ID
            scans_df["retention_time_seconds"] = spectrum.scan_time_in_minutes() * 60.0
            scans_df["filename"] = Path(input_file).resolve().name
            scans.append(scans_df)

    scans_out = pd.concat(scans)
    logging.info("Scans info extraction completed!")
    return scans_out
