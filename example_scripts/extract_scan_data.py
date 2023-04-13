import argparse

from simepy.extract_scans import extract_scan_data


def main(input_file, output_file):
    """
    Extract scan info from mzml file using extract_scan_data().

    Args:
        input_file (PosixPath): path to input file
        output_file (PosixPath): path to scans output file
    """
    scans_out = extract_scan_data(input_file)
    scans_out.to_csv(output_file, index=False)


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

    main(input_file=args.input_file, output_file=args.output_file)
