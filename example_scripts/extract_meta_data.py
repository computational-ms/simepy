import argparse
from pathlib import Path

from simepy.extract_meta_data import extract_meta_data


def main(
    input_file,
    run_output_file,
    spec_output_file,
    noise_output_file,
    unit_output_file,
    time_format,
    object_name,
    lineage_root,
):
    """
    Extract scan info from mzml file using extract_scan_data().

    Args:
        input_file (PosixPath): path to input file
        run_output_file (str/NoneType): path to run metadata output_file
        spec_output_file (str/NoneType): path to spectrum metadata output_file
        noise_output_file (str/NoneType): path to spectrum noise output_file
        units_output_file (str/NoneType): path to instrument units output_file
        time_format (str): string defining the time_format to format into
        object_name (str/NoneType): string to be used as file identifier
        lineage_root (str/NoneType): name of the root file associated to the one, from which metadata is queried
    """
    rmd_df, imu_df, smd_df, sn_df = extract_meta_data(
        input_file=input_file,
        time_format=time_format,
        object_name=object_name,
        lineage_root=lineage_root,
    )

    if run_output_file is not None:
        rmd_df.to_csv(run_output_file, index=False)

    if spec_output_file is not None:
        smd_df.to_csv(spec_output_file, index=False)

    if unit_output_file is not None:
        imu_df.to_csv(unit_output_file, index=False)

    if noise_output_file is not None:
        sn_df.to_csv(noise_output_file, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_file",
        dest="input_file",
        help="input file",
    )
    parser.add_argument(
        "-ro",
        "--run_output_file",
        type=lambda x: None if x == "None" else str(x),
        dest="run_output_file",
        help="run metadata output file",
    )
    parser.add_argument(
        "-so",
        "--spec_output_file",
        type=lambda x: None if x == "None" else str(x),
        dest="spec_output_file",
        help="spectrum metadata output file",
    )
    parser.add_argument(
        "-sno",
        "--spec_noise_output_file",
        type=lambda x: None if x == "None" else str(x),
        dest="spec_noise_output_file",
        help="spectrum noise output file",
    )
    parser.add_argument(
        "-iuo",
        "--instrument_unit_output_file",
        type=lambda x: None if x == "None" else str(x),
        dest="instrument_unit_output_file",
        help="instrument unit output file",
    )
    parser.add_argument(
        "-tf",
        "--time_format",
        default="%Y-%m-%d %H:%M:%S",
        dest="time_format",
        help="time format string",
    )
    parser.add_argument(
        "-on",
        "--object_name",
        dest="object_name",
        help="object_name string",
    )
    parser.add_argument(
        "-lr",
        "--lineage_root",
        dest="lineage_root",
        help="lineage_root string",
    )
    args = parser.parse_args()

    if args.object_name is None:
        args.object_name = Path(args.input_file).name
    if args.lineage_root is None:
        args.lineage_root = Path(args.input_file).resolve()
    main(
        input_file=args.input_file,
        run_output_file=args.run_output_file,
        spec_output_file=args.spec_output_file,
        noise_output_file=args.spec_noise_output_file,
        unit_output_file=args.instrument_unit_output_file,
        time_format=args.time_format,
        object_name=args.object_name,
        lineage_root=args.lineage_root,
    )
