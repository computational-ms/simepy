import argparse
from pathlib import Path
from spectrum_data.extract_meta_data import extract_meta_data


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
        dest="run_output_file",
        help="run metadata output file",
    )
    parser.add_argument(
        "-so",
        "--spec_output_file",
        dest="spec_output_file",
        help="spectrum metadata output file",
    )
    parser.add_argument(
        "-sno",
        "--spec_noise_output_file",
        dest="spec_noise_output_file",
        help="spectrum noise output file",
    )
    parser.add_argument(
        "-iuo",
        "--instrument_unit_output_file",
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
    extract_meta_data(
        input_file=args.input_file,
        run_output_file=args.run_output_file,
        spec_output_file=args.spec_output_file,
        noise_output_file=args.spec_noise_output_file,
        unit_output_file=args.instrument_unit_output_file,
        time_format=args.time_format,
        object_name=args.object_name,
        lineage_root=args.lineage_root,
    )
