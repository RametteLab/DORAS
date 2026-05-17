#!/usr/bin/env python
import argparse
import sys
import gzip
from datetime import datetime, timedelta
import logging
import os

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def parse_timestamp(header):
    """Extracts the timestamp from a FASTQ header."""
    for part in header.split(" "):
        if "start_time=" in part:
            return datetime.fromisoformat(
                part.replace("start_time=", "").replace("+01:00", "")
            )
    raise ValueError


def filter_reads_by_time(
    input_fastq, output_file, start_time, end_time, gzip_input=False
):
    """Filters reads based on the specified time range and writes them to an output file."""

    # Determine the earliest timestamp
    min_timestamp = datetime.fromtimestamp(0)
    max_timestamp = datetime.fromtimestamp(0)
    if gzip_input:
        # with gzip.GzipFile(mode='r', filename=input_fastq.name) as infile:
        #     while True:
        #         try:
        #             header = next(infile).strip()
        #         except StopIteration:
        #             break

        #         if not header:
        #             break

        #         sequence = next(infile).strip()
        #         plus_line = next(infile).strip()  # This line is typically just a '+'
        #         quality_scores = next(infile).strip()

        #         try:
        #             timestamp = parse_timestamp(header)
        #             if min_timestamp is None or timestamp < min_timestamp:
        #                 min_timestamp = timestamp
        #         except ValueError as e:
        #             print(e, file=sys.stderr)
        with gzip.open(input_fastq.name, "rt") as infile:
            for line in infile:
                # header = next(line).strip()
                # print(header)
                if line.startswith("@"):
                    header = line.strip()
                    try:
                        timestamp = parse_timestamp(header)
                        if (
                            min_timestamp == datetime.fromtimestamp(0)
                            or timestamp < min_timestamp
                        ):
                            min_timestamp = timestamp
                        if (
                            max_timestamp == datetime.fromtimestamp(0)
                            or timestamp > max_timestamp
                        ):
                            max_timestamp = timestamp

                    except ValueError as e:
                        logging.debug(e)
    else:
        with open(input_fastq.name) as infile:
            for line in infile:
                if line.startswith("@"):
                    header = line.strip()
                    try:
                        timestamp = parse_timestamp(header)
                        if (
                            min_timestamp == datetime.fromtimestamp(0)
                            or timestamp < min_timestamp
                        ):
                            min_timestamp = timestamp
                        if (
                            max_timestamp == datetime.fromtimestamp(0)
                            or timestamp > max_timestamp
                        ):
                            max_timestamp = timestamp
                    except ValueError as e:
                        logging.debug(e)

    # Calculate the start and end timestamps based on the provided time range
    if args.info:
        print(
            f"Min timestamp: {min_timestamp} Max timestamp: {max_timestamp} , duration is {max_timestamp - min_timestamp}"
        )
        sys.exit(0)
    start_timestamp = min_timestamp + timedelta(hours=start_time)
    end_timestamp = min_timestamp + timedelta(hours=end_time)
    # verify whether max_timestamp is below the end_timestamp if so warning or exit
    if max_timestamp < end_timestamp:
        logging.error(
            f"Endtime is not in range of the data, max is {max_timestamp} and endpoint is {end_timestamp}\n The experiment ran from {min_timestamp} to {max_timestamp}"
        )
        end_timestamp = min(end_timestamp, max_timestamp)
        logging.warning(
            f"Endtime is not in range of the data, max is {max_timestamp} and endpoint is {end_timestamp}\n The experiment ran from {min_timestamp} to {max_timestamp}"
        )
        # sys.exit(1)
    logging.info(f"Start Time chosen: {start_timestamp}")
    logging.info(f"End Time chosen: {end_timestamp}")

    # Filter reads within this range and write to output file
    if gzip_input:
        with gzip.open(input_fastq.name, "rt") as infile:
            with gzip.GzipFile(fileobj=output_file, mode="wb", filename="") as outfile:
                while True:
                    try:
                        header = next(infile).strip()
                    except StopIteration:
                        break

                    if not header or not header.startswith("@"):
                        break

                    sequence = next(infile).strip()
                    plus_line = next(
                        infile
                    ).strip()  # This line is typically just a '+'
                    quality_scores = next(infile).strip()

                    if not (
                        header and sequence and plus_line == "+" and quality_scores
                    ):
                        raise ValueError("Unexpected FASTQ format")

                    timestamp = parse_timestamp(header)

                    if start_timestamp <= timestamp <= end_timestamp:
                        outfile.write(
                            f"{header}\n{sequence}\n+\n{quality_scores}\n".encode(
                                "utf-8"
                            )
                        )
    else:
        with open(input_fastq.name) as infile:
            with gzip.GzipFile(fileobj=output_file, mode="wb", filename="") as outfile:
                while True:
                    header = next(infile).strip()
                    if not header or not header.startswith("@"):
                        break

                    sequence = next(infile).strip()
                    plus_line = next(
                        infile
                    ).strip()  # This line is typically just a '+'
                    quality_scores = next(infile).strip()

                    if not (
                        header and sequence and plus_line == "+" and quality_scores
                    ):
                        raise ValueError("Unexpected FASTQ format")

                    timestamp = parse_timestamp(header)

                    if start_timestamp <= timestamp <= end_timestamp:
                        outfile.write(
                            f"{header}\n{sequence}\n+\n{quality_scores}\n".encode(
                                "utf-8"
                            )
                        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter reads from a FASTQ file based on timestamps. Requires the start_time flag in the header (ONT Only)"
    )

    # Add arguments for input and output files
    parser.add_argument(
        "--input",
        "-i",
        type=argparse.FileType("r"),
        default=sys.stdin,
        help="Input FASTQ file (can be gzipped)",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default="-",
        help="Output filtered FASTQ file (default: stdout.gz if output is a file, otherwise stdout)",
    )

    # Time range arguments
    parser.add_argument(
        "start_time",
        type=float,
        nargs="?",
        default=0.0,
        help="Initial time point (hours), accepts floats. Default is 0.0.",
    )
    parser.add_argument(
        "end_time",
        type=float,
        help="Last time point (hours), accepts floats. Default will go from start to this time point.",
    )

    # Optional argument to specify if the input is a gzipped file
    parser.add_argument(
        "--gzip-input",
        action="store_true",
        default=False,
        help="Indicate that the input FASTQ file is gzipped",
    )
    parser.add_argument(
        "--info",
        action="store_true",
        help="Return the min and max time and duration only.",
    )

    args = parser.parse_args()

    # Determine output file type
    if args.output == "-":
        output_file = sys.stdout.buffer  # Use stdout in binary mode for gzip
    else:
        output_dir = os.path.dirname(args.output)
        if output_dir and not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir, exist_ok=True)
            except OSError as e:
                print(f"Error creating directory {output_dir}: {e}", file=sys.stderr)
                sys.exit(1)
        output_file = open(
            args.output + ".gz", "wb"
        )  # Open the file in binary write mode

    # Debugging: Print the input and output paths
    print(f"Input file: {args.input.name}", file=sys.stderr)
    print(
        f"Output file: {output_file.name if hasattr(output_file, 'name') else 'stdout'}",
        file=sys.stderr,
    )

    filter_reads_by_time(
        args.input,
        output_file,
        args.start_time,
        args.end_time,
        gzip_input=args.gzip_input,
    )

    # Close the output file if it's not stdout
    if output_file is not sys.stdout.buffer:
        output_file.close()
