import argparse
from simepy.extract_scans import extract_scan_data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_file",
        dest="input_file",
        help="mzml input file",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        dest="output_file",
        help="scans output file",
    )
    args = parser.parse_args()

    extract_scan_data(input_file=args.input_file, output_file=args.output_file)
