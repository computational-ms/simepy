"""Spectrum Meta Data resource."""
import logging
from datetime import datetime, timedelta
from pathlib import Path
from typing import Optional

import pandas as pd
import pymzml
from dateutil import parser as dparser
from pydantic import BaseModel, Field

try:
    import czxcalibur

    raw_readable = True

except ModuleNotFoundError:
    raw_readable = False


class SpectrumMetaData(BaseModel):
    """Spectrum Meta Data extraction class."""

    lineage_root: str = Field(description="Root files of the execution graph")
    file: str = Field(description="Name of the file containing spectra")
    spectrum_id: int = Field(description="data file spectrum identifier")
    ms_level: int = Field(description="the ms level of the spectrum")
    rt: float = Field(description="retention time of the spectrum")
    rt_unit: str = Field(description="time unit for the rt")
    preset_scan_configuration: Optional[int] = Field(
        default=None,
        description='A user-defined scan configuration that specifies the instrumental settings in which a spectrum is acquired. An instrument may cycle through a list of preset scan configurations to acquire data. This is a more generic term for the Thermo "scan event", which is defined in the Thermo Xcalibur glossary as: a mass spectrometer scan that is defined by choosing the necessary scan parameter settings. Multiple scan events can be defined for each segment of time.',
    )
    total_ion_current: Optional[float] = Field(
        default=None, description="summed spectrum intensities"
    )
    ion_injection_time: Optional[float] = Field(
        default=None, description="the time taken to accumulate ions for the spectrum"
    )
    ion_injection_time_unit: Optional[str] = Field(
        default=None, description="time unit for the ion_injection_time"
    )
    filter_string: Optional[str] = Field(
        default=None, description="Xcalibur spectrum header string"
    )
    scan_window_lower_limit: Optional[float] = Field(
        default=None,
        description="The lower m/z bound of a mass spectrometer scan window",
    )
    scan_window_upper_limit: Optional[float] = Field(
        default=None,
        description="The upper m/z bound of a mass spectrometer scan window",
    )
    precursor_mz: Optional[float] = Field(
        default=None, description="Xcalibur recalculated precursor m/z"
    )
    precursor_charge: Optional[int] = Field(
        default=None, description="charge state of the precursor ion"
    )
    # precursor_intensity: Optional[float]
    #  ^--- wrong in mzML .. dunno where in raw
    precursor_id: Optional[int] = Field(
        default=None, description="assigned id for precursor ion"
    )
    resolution: Optional[int] = Field(
        default=None, description="spectrum nominal resolution"
    )
    ms2_isolation_width: Optional[float] = Field(
        default=None, description="isolation width for precursor ion"
    )
    ms2_isolation_mz: Optional[float] = Field(
        default=None, description="Xcalibur set mass for the isolation of the precursor"
    )
    gex_id: Optional[int] = Field(
        default=None,
        description="the gex_id (number) is used to separate different cycles or 'experiments' applied during acquisition.",
    )
    faims_cv: Optional[float] = Field(
        default=None, description="the FAIMS voltage applied"
    )


class SpectrumNoise(BaseModel):
    """Individual entry data for noise data."""

    spectrum_id: int = Field(description="data file spectrum identifier")
    mz: float = Field(description="the m/z value for the noise measurement")
    noise: float = Field(description="the measured noise value at the m/z value")
    baseline: float = Field(description="the Xcalibur calculated baseline")


class InstrumentMetaUnit(BaseModel):
    """Individual entry data for the instrument unit (MS/MS scanning parameters)."""

    unit: int = Field(description="the id of the unit")
    gex_id: int = Field(
        description="the gex_id (number) is used to separate different cycles or 'experiments' applied during acquisition."
    )
    order: int = Field(description="the order in which the gex is processed")
    ms_level: int = Field(description="the ms level of the unit spectra")
    activation: str = Field(description="the type of excitation used for fragmentation")
    detector: str = Field(
        description="the section of the instrument used to measure the spectrum"
    )
    activation_time: float = Field(
        description="the activation time, used in fragmentation"
    )
    energy: float = Field(description="the energy used in fragmentation")
    isolation: float = Field(description="the instrument ioslation width setting")
    low_mz: float = Field(description="the lowest m/z in the spectrum")
    resolution: int = Field(description="the resolution setting for the spectrum")
    use: str = Field(description="the fragmentation and spectrum usage regime")
    scan_events: Optional[str] = Field(
        default=None,
        description="a list of all spectrum orders that are part of the gex",
    )
    col_energy_steps: Optional[int] = Field(
        default=None,
        description="the number of steps used to increase the collision energy",
    )
    col_energy_width: Optional[float] = Field(
        default=None, description="the size of each step in stepped collisions"
    )
    faims_cv: float = Field(description="the FAIMS voltage applied")


###----------


class RunInfoMetaData(BaseModel):
    """RunInfo Meta Data extraction class."""

    lineage_root: str = Field(description="Root files of the execution graph")
    file: str = Field(description="Name of the file containing spectra")
    run_id: str = Field(description="Unique identifier for the run")
    spectrum_count: int = Field(description="Number of spectra in run")
    file_size_mb: float = Field(description="Size of the file in MB")
    start_time: str = Field(description="Start time of the run")
    stop_time: Optional[str] = Field(default=None, description="Stop time of the run")


###----------


### mzML


def _extract_spectrum_info_from_mzml(
    file_path, object_name, lineage_root, time_format, gex_lookup=None
):
    """Extract spectrum meta data from the mzML file.

    Args:
        input_file (PosixPath): oath to input file
        object_name (str): Object Name in a generalized way
        lineage_root (str): name of the root file associated to the one, from which metadata is queried
        time_format (str): string defining the time_format to format into
        gex_lookup: (dict): containing the faims_cv voltage to gex_id lookup values

    Returns:
        Populated SpectrumMetaData object
    """
    with pymzml.run.Reader(str(file_path)) as reader:
        for spectrum in reader:
            scan_time, scan_time_unit = spectrum.scan_time
            # scan time unit should be either second or minute
            # see https://github.com/HUPO-PSI/psi-ms-CV/blob/fa10067dfe12c26840bd1a87d0a776fb15d0f20a/psi-ms.obo#L308-L317
            scan_time_unit = scan_time_unit.lower()
            if scan_time_unit == "minute":
                scan_time_in_seconds = scan_time * 60
                scan_time_unit = "second"
            elif scan_time_unit == "second":
                scan_time_in_seconds = scan_time
            else:
                message = f"Cant recognize unit {scan_time_unit}"
                logging.exception(message)
                raise Exception(message)

            data = {
                "lineage_root": str(lineage_root),
                "file": str(object_name),
                "spectrum_id": spectrum.ID,
                "ms_level": spectrum.ms_level,
                "rt": scan_time_in_seconds,
                "rt_unit": scan_time_unit,
            }

            element_lookup = [
                "preset scan configuration",
                "ion injection time",
                "total ion current",
                "filter string",
                "scan window lower limit",
                "scan window upper limit",
            ]

            for val in element_lookup:
                ele = spectrum.get_element_by_name(val)
                if ele != None:
                    data.update({val.replace(" ", "_"): ele.get("value")})
                    if val == "ion injection time":
                        data.update({"ion_injection_time_unit": ele.get("unitName")})

            if len(spectrum.selected_precursors) > 0:
                for n, precursor in enumerate(spectrum.selected_precursors):
                    data["precursor_mz"] = precursor.get("mz", "")
                    data["precursor_charge"] = precursor.get("charge", None)
                    data["precursor_intensity"] = precursor.get("i", None)
                    # ^--- this value is not represetned in the MS1 scan and I fail to understand
                    #      where this comes from, looking at two files
                    data["precursor_id"] = n
                    yield SpectrumMetaData(**data)
            else:
                yield SpectrumMetaData(**data)


def _extract_run_info_from_mzml(
    file_path, object_name, lineage_root, time_format, gex_lookup=None
):
    """Extract simple run data from the mzML file.

    Args:
        input_file (PosixPath): oath to input file
        object_name (str): Object Name in a generalized way
        lineage_root (str): name of the root file associated to the one, from which metadata is queried
        time_format (str): string defining the time_format to format into
        gex_lookup: (dict): containing the faims_cv voltage to gex_id lookup values

    Returns:
        Populated RunInfoMetaData object
    """
    with pymzml.run.Reader(str(file_path)) as reader:
        start_time = dparser.isoparse(reader.info["start_time"]).strftime(time_format)

        data = {
            "lineage_root": str(lineage_root),
            "file": str(object_name),
            "run_id": reader.info["run_id"],
            "spectrum_count": reader.info["spectrum_count"],
            "file_size_mb": file_path.stat().st_size / 1e6,
            "start_time": start_time,
        }

        # Try to estimate the stop time!
        # The last TIC, is the start time in minutes relative to the run start time.
        # Format it to seconds (*60) for easier timedelta calculation.
        # could also keep it in minutes and use timedelta with minutes=last_spec_time
        try:
            last_spec_time = round(reader["TIC"].time[-1] * 60, 2)
        except (StopIteration, KeyError):
            logging.warning(
                "Your mzml does not contain TIC information. Therefore, the end time of the ms run cannot be estimated and will be set to None."
            )
            last_spec_time = None

        if last_spec_time is not None:
            stop_time = datetime.strptime(start_time, time_format) + timedelta(
                seconds=last_spec_time
            )
            data.update({"stop_time": stop_time.strftime(time_format)})

        yield RunInfoMetaData(**data)


### RAW


def _extract_run_info_from_raw(
    input_file, object_name, lineage_root, time_format, gex_lookup=None
):
    """Extract simple run data from the Xcalibur file.

    Args:
        input_file (PosixPath): oath to input file
        object_name (str): Object Name in a generalized way
        lineage_root (str): name of the root file associated to the one, from which metadata is queried
        time_format (str): string defining the time_format to format into

    Returns:
        Populated RunInfoMetaData object
    """
    raw_io = czxcalibur.CZXRawFile(str(input_file))
    scan_range = raw_io.scan_range
    start_time = raw_io.get_creation_date()
    last_scan_time = raw_io.get_rt_from_scan_number(scan_range[1])
    end_time = start_time + timedelta(
        minutes=raw_io.get_rt_from_scan_number(scan_range[1])
    )
    data = {
        "lineage_root": str(lineage_root),
        "file": str(object_name),
        "run_id": Path(raw_io.filename).name,
        "spectrum_count": scan_range[1],
        "start_time": start_time.strftime(time_format),
        "stop_time": end_time.strftime(time_format),
        "file_size_mb": Path(input_file).stat().st_size / 1e6,
    }

    yield RunInfoMetaData(**data)


def _extract_units_info_from_raw(
    input_file, object_name, lineage_root, time_format, gex_lookup=None
):
    """Extract the  from the Xcalibur file.

    Args:
        input_file (PosixPath): oath to input file
        object_name (str): Object Name in a generalized way
        lineage_root (str): name of the root file associated to the one, from which metadata is queried
        time_format (str): string defining the time_format to format into
        gex_lookup: (dict): containing the faims_cv voltage to gex_id lookup values

    Returns:
        Populated RunInfoMetaData object
    """
    raw_io = czxcalibur.CZXRawFile(str(input_file))
    method = raw_io.get_formatted_method_data()
    events = raw_io.extract_instrument_events(method)

    mapping = {
        "acttime": "activation_time",
        "lowmz": "low_mz",
        "scanevents": "scan_events",
        "colenergysteps": "col_energy_steps",
        "colenergywidth": "col_energy_width",
    }
    unit_id = 0
    for experiment in events:
        data = events[experiment]
        if len(data["msms_event"]) == 1:
            use = "IQ"
        else:
            use = "I"

        for idx, unit in enumerate(data["msms_event"]):
            unit_id += 1
            unit["unit"] = unit_id
            unit["order"] = idx + 1
            if unit["ms_level"] == 2:
                unit["use"] = use
            else:
                unit["use"] = "Q"

            for key in mapping:
                if key in unit:
                    unit[mapping[key]] = unit[key]
                    del unit[key]
            yield InstrumentMetaUnit(**unit)


def _extract_spectrum_info_from_raw(
    input_file, object_name, lineage_root, time_format, gex_lookup=None
):
    """Extract spectrum meta data from the Xcalibur file.

    Args:
        input_file (PosixPath): oath to input file
        object_name (str): Object Name in a generalized way
        lineage_root (str): name of the root file associated to the one, from which metadata is queried
        time_format (str): string defining the time_format to format into
        gex_lookup: (dict): containing the faims_cv voltage to gex_id lookup values

    Returns:
        Populated SpectrumMetaData object
    """
    raw_io = czxcalibur.CZXRawFile(str(input_file))
    levels = {"ms": 1, "ms2": 2, "ms3": 3}
    for spectrum_id in raw_io.scan_range:
        filter = raw_io.get_parsed_filter(spectrum_id)
        extra = raw_io.get_formatted_trailer_extra_for_scan_number(spectrum_id)
        header = raw_io.get_scan_header_info_for_scan_number(spectrum_id)
        try:
            faims_cv = float(filter["cv"])
            gex_id = gex_lookup[faims_cv]
        except (TypeError, KeyError, ValueError):
            faims_cv = 0.0
            gex_id = 1

        precursor_mz = extra.get("Monoisotopic M/Z", None)
        if precursor_mz == 0:
            if levels[filter["scan"]] == 1:
                precursor_mz = None
            else:
                precursor_mz = float(filter["setmass1"])

        data = {
            "lineage_root": str(lineage_root),
            "file": str(input_file),
            "spectrum_id": spectrum_id,
            "ms_level": levels[filter["scan"]],
            "rt": header["start_time"] * 60,
            "rt_unit": "second",
            # "preset_scan_configuration": Optional[int],
            "total_ion_current": header["tic"],
            "ion_injection_time": extra["Ion Injection Time (ms)"],
            "ion_injection_time_unit": "ms",
            "filter_string": filter["filter_string"],
            "scan_window_lower_limit": header["low_mass"],
            "scan_window_upper_limit": header["high_mass"],
            "precursor_mz": precursor_mz,
            "precursor_charge": extra.get("Charge State", None),
            "precursor_intensity": None,
            "precursor_id": None,
            "resolution": extra.get("Orbitrap Resolution", None),
            "ms2_isolation_width": extra.get("MS2 Isolation Width", None),
            "ms2_isolation_mz": filter["setmass1"],
            "gex_id": gex_id,
            "faims_cv": faims_cv,
        }
        yield SpectrumMetaData(**data)


def _extract_spectrum_noise_from_raw(
    input_file, object_name, lineage_root, time_format, gex_lookup=None
):
    """Extract spectrum noise data from the Xcalibur file.

    Args:
        input_file (PosixPath): oath to input file
        object_name (str): Object Name in a generalized way
        lineage_root (str): name of the root file associated to the one, from which metadata is queried
        time_format (str): string defining the time_format to format into
        gex_lookup: (dict): containing the faims_cv voltage to gex_id lookup values

    Returns:
        Populated SpectrumNoise object
    """
    raw_io = czxcalibur.CZXRawFile(str(input_file))

    for scan in raw_io.scan_range:
        noise_list = raw_io.get_noise_data_for_scan_number(scan)
        for noise_data in noise_list:
            data = {
                "spectrum_id": scan,
                "mz": noise_data[0],
                "noise": noise_data[1],
                "baseline": noise_data[2],
            }

            yield SpectrumNoise(**data)


def extract_data(
    input_file, object_name, lineage_root, file_type, time_format, mode, gex_lookup=None
):
    """
    Return Meta data extraction function.

    Args:
        input_file (PosixPath): oath to input file
        object_name (str): Object Name in a generalized way
            (no abs path but rather data location independent view of and object.
            Think more url, urn in ursgal terms, where object_name is the urn)
        lineage_root (str): name of the root file associated to the one, from which metadata is queried
        file_type (str): file type of the input file for which the extraction function should be determined
        time_format (str): string defining the time_format to format into
        mode (str): defines if spectrum_metadata or run_metadata dict should be used to identify the correct extraction function
        gex_lookup: (dict): containing the faims_cv voltage to gex_id lookup values

    Raises:
        IOError: if file type cannot be used to identify appropriate
            metadata extraction function

    Returns:
        extraction function: function will iter over file and return
            metadata object for each spectrum
    """
    spec_md_extractors = {
        "mzml": _extract_spectrum_info_from_mzml,
        "mzml.gz": _extract_spectrum_info_from_mzml,
        "raw": _extract_spectrum_info_from_raw,
        # "mgf": _extract_spectrum_info_from_mgf, #Not implemented yet!
    }

    run_md_extractors = {
        "mzml": _extract_run_info_from_mzml,
        "mzml.gz": _extract_run_info_from_mzml,
        "raw": _extract_run_info_from_raw,
        # "mgf": _extract_run_info_from_mgf, #Not implemented yet!
    }

    spec_noise_extractors = {
        "raw": _extract_spectrum_noise_from_raw,
    }

    instrument_units_extractors = {
        "raw": _extract_units_info_from_raw,
    }

    if mode == "spectrum_metadata":
        extract_function = spec_md_extractors.get(file_type, None)
    elif mode == "run_metadata":
        extract_function = run_md_extractors.get(file_type, None)
    elif mode == "run_noise":
        extract_function = spec_noise_extractors.get(file_type, None)
    elif mode == "instrument_unit":
        extract_function = instrument_units_extractors.get(file_type, None)

    if extract_function is None:
        raise NotImplementedError(
            f"{mode} extraction from {file_type} files is not implemented!"
        )
    elif file_type == "raw" and raw_readable is False:
        logging.error("CZXcalibur not installed! Require czxcalibur to read raw.")
        raise ModuleNotFoundError
    else:
        return extract_function(input_file, object_name, lineage_root, time_format)


def extract_meta_data(
    input_file,
    time_format="%Y-%m-%d %H:%M:%S",
    object_name=None,
    lineage_root=None,
):
    """
    Extract Run and Spectrum metadata from MS data.

    Args:
        input_file (PosixPath/str): path to ms data input file (mzml|raw|mgf)
        time_format (str): string defining the time_format to format into
        object_name (str/NoneType): string to be used as file identifier
        lineage_root (str/NoneType): name of the root file associated to the one, from which metadata is queried

    Returns:
        rmd_df (pandas.DataFrame): df containing run metadata
        imu_df (pandas.DataFrame): df containing instrument unit info
        smd_df (pandas.DataFrame): df containing spectrum metadata
        sn_df (pandas.DataFrame): df containing spectrum noise info
    """
    if isinstance(input_file, str):
        input_file = Path(input_file)

    if object_name is None:
        object_name = Path(input_file).resolve().name

    if lineage_root is None:
        lineage_root = Path(input_file).resolve().name

    file_type = "".join(Path(input_file).suffixes).lower()[1:]
    logging.info(
        f"Identified file type is '{file_type}'. Proceeding with metadata extraction..."
    )
    # 1. RunInfoMetaData
    rmd_rows = []
    for meta_data_line in extract_data(
        input_file,
        object_name,
        lineage_root,
        file_type,
        time_format,
        mode="run_metadata",
    ):
        rmd_rows.append(
            pd.DataFrame(meta_data_line.dict(exclude_unset=True), index=[0])
        )
    rmd_df = pd.concat(rmd_rows, ignore_index=True)
    logging.info("Finished run metadata info extraction!")

    # 2. InstrumentUnit
    gex_lookup = {}
    if input_file.suffix.lower() == ".raw":
        imu_rows = []
        for unit_data_line in extract_data(
            input_file,
            object_name,
            lineage_root,
            file_type,
            time_format,
            mode="instrument_unit",
        ):
            imu_rows.append(
                pd.DataFrame(unit_data_line.dict(exclude_unset=True), index=[0])
            )
            gex_lookup[unit_data_line.faims_cv] = unit_data_line.gex_id
        imu_df = pd.concat(imu_rows, ignore_index=True)
        logging.info("Finished instrument unit info extraction!")

    elif "".join(input_file.suffixes).lower() in [".mzml", ".mzml.gz"]:
        imu_df = pd.DataFrame(
            [], columns=[x for x in InstrumentMetaUnit.schema()["properties"].keys()]
        )

    # 3. SpectrumMetaData
    smd_rows = []
    for meta_data_line in extract_data(
        input_file,
        object_name,
        lineage_root,
        file_type,
        time_format,
        mode="spectrum_metadata",
        gex_lookup=gex_lookup,
    ):
        smd_rows.append(
            pd.DataFrame(meta_data_line.dict(exclude_unset=True), index=[0])
        )
    smd_df = pd.concat(smd_rows, ignore_index=True)
    logging.info("Finished spectrum metadata info extraction!")

    # 4. SpectrumNoise
    if input_file.suffix.lower() == ".raw":
        sn_rows = []
        for noise_data_line in extract_data(
            input_file,
            object_name,
            lineage_root,
            file_type,
            time_format,
            mode="run_noise",
        ):
            sn_rows.append(
                pd.DataFrame(noise_data_line.dict(exclude_unset=True), index=[0])
            )
        sn_df = pd.concat(sn_rows, ignore_index=True)
        logging.info("Finished spectrum noise extraction!")

    elif "".join(input_file.suffixes).lower() in [".mzml", ".mzml.gz"]:
        sn_df = pd.DataFrame(
            [], columns=[x for x in SpectrumNoise.schema()["properties"].keys()]
        )

    return rmd_df, imu_df, smd_df, sn_df
