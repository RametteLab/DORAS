#!/usr/bin/env python
import gzip
from datetime import datetime, timedelta
import logging
import asyncio
from pathlib import Path

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


class FastQTimeExtractor:
    """
    A class to handle FASTQ files with timestamps and extract reads based on time ranges.
    """

    def __init__(self, input_fastq_path, end_point=None):
        """
        Initialize the extractor with a FASTQ file path.
        """
        self.input_fastq_path = input_fastq_path
        self.min_timestamp = None
        self.max_timestamp = None
        self.end_point = end_point
        self.total_reads = 0
        self.gzip_input = input_fastq_path.endswith(".gz")

    def parse_timestamp(self, header):
        """Extracts the timestamp from a FASTQ header."""
        for part in header.split(" "):
            if "start_time=" in part:
                return datetime.fromisoformat(
                    part.replace("start_time=", "").replace("+01:00", "")
                )
        raise ValueError(f"Timestamp not found in header: {header}")

    def get_run_info(self):
        """
        Analyze the FASTQ file and extract run information.
        """
        self.min_timestamp = datetime.fromtimestamp(0)
        self.max_timestamp = datetime.fromtimestamp(0)
        self.total_reads = 0

        if self.gzip_input:
            with gzip.open(self.input_fastq_path, "rt") as infile:
                for line in infile:
                    if line.startswith("@"):
                        header = line.strip()
                        try:
                            timestamp = self.parse_timestamp(header)
                            if (
                                self.min_timestamp == datetime.fromtimestamp(0)
                                or timestamp < self.min_timestamp
                            ):
                                self.min_timestamp = timestamp
                            if (
                                self.max_timestamp == datetime.fromtimestamp(0)
                                or timestamp > self.max_timestamp
                            ):
                                self.max_timestamp = timestamp
                            self.total_reads += 1
                        except ValueError as e:
                            logging.debug(f"Error parsing header: {e}")
        else:
            with open(self.input_fastq_path) as infile:
                for line in infile:
                    if line.startswith("@"):
                        header = line.strip()
                        try:
                            timestamp = self.parse_timestamp(header)
                            if (
                                self.min_timestamp == datetime.fromtimestamp(0)
                                or timestamp < self.min_timestamp
                            ):
                                self.min_timestamp = timestamp
                            if (
                                self.max_timestamp == datetime.fromtimestamp(0)
                                or timestamp > self.max_timestamp
                            ):
                                self.max_timestamp = timestamp
                            self.total_reads += 1
                        except ValueError as e:
                            logging.debug(f"Error parsing header: {e}")

        duration = self.max_timestamp - self.min_timestamp
        logging.info(
            f"Run info - Min timestamp: {self.min_timestamp}, Max timestamp: {self.max_timestamp}, Total reads: {self.total_reads}, Duration: {duration}"
        )
        return {
            "min_timestamp": self.min_timestamp,
            "max_timestamp": self.max_timestamp,
            "total_reads": self.total_reads,
            "duration": duration,
        }

    def extract_by_time_range(self, start_hours=0.0, end_hours=None, output_file=None):
        """
        Extract reads from the specified time range.
        """
        if self.min_timestamp is None:
            self.get_run_info()

        start_timestamp = self.min_timestamp + timedelta(hours=start_hours)
        end_timestamp = (
            self.min_timestamp + timedelta(hours=end_hours)
            if end_hours is not None
            else self.max_timestamp
        )

        if (
            end_hours is not None
            and end_hours
            > (self.max_timestamp - self.min_timestamp).total_seconds() / 3600
        ):
            logging.warning(
                f"End hours {end_hours} exceeds experiment duration. Using max timestamp instead."
            )
            end_timestamp = self.max_timestamp

        extracted_reads = []

        if self.gzip_input:
            with gzip.open(self.input_fastq_path, "rt") as infile:
                for line in infile:
                    if line.startswith("@"):
                        header = line.strip()
                        try:
                            timestamp = self.parse_timestamp(header)
                            if start_timestamp <= timestamp <= end_timestamp:
                                sequence = next(infile).strip()
                                plus_line = next(infile).strip()
                                quality_scores = next(infile).strip()

                                if not (
                                    header
                                    and sequence
                                    and plus_line == "+"
                                    and quality_scores
                                ):
                                    raise ValueError("Unexpected FASTQ format")

                                extracted_reads.append(
                                    (header, sequence, plus_line, quality_scores)
                                )
                        except ValueError as e:
                            logging.debug(f"Error parsing header: {e}")
        else:
            with open(self.input_fastq_path) as infile:
                for line in infile:
                    if line.startswith("@"):
                        header = line.strip()
                        try:
                            timestamp = self.parse_timestamp(header)
                            if start_timestamp <= timestamp <= end_timestamp:
                                sequence = next(infile).strip()
                                plus_line = next(infile).strip()
                                quality_scores = next(infile).strip()

                                if not (
                                    header
                                    and sequence
                                    and plus_line == "+"
                                    and quality_scores
                                ):
                                    raise ValueError("Unexpected FASTQ format")

                                extracted_reads.append(
                                    (header, sequence, plus_line, quality_scores)
                                )
                        except ValueError as e:
                            logging.debug(f"Error parsing header: {e}")

        if output_file:
            with gzip.open(output_file, "wb") as outfile:
                for read in extracted_reads:
                    header, sequence, plus_line, quality_scores = read
                    outfile.write(
                        f"{header}\n{sequence}\n{plus_line}\n{quality_scores}\n".encode(
                            "utf-8"
                        )
                    )

        return extracted_reads

    async def create_mock_files(
        self,
        output_dir,
        interval_hours=1.0,
        start_time_point=0,
        max_records_per_file=0,
        stop_event=None,
    ):
        """
        Create mock FASTQ files with time-based windows from the original experiment.
        """
        if self.min_timestamp is None:
            self.get_run_info()

        if not Path(output_dir).exists():
            Path(output_dir).mkdir(parents=True, exist_ok=True)
        if start_time_point > 0:
            logging.info(
                f"Adjusted min timestamp to start time point: {self.min_timestamp + timedelta(hours=start_time_point)}"
            )
        # Read and parse all records with timestamps
        timestamped_records = []
        if self.gzip_input:
            with gzip.open(self.input_fastq_path, "rt") as infile:
                while True:
                    header_line = infile.readline()
                    if not header_line:
                        break
                    if header_line.startswith("@"):
                        # print(header_line)
                        header = header_line.strip()
                        sequence = infile.readline()
                        sequence = sequence.strip()
                        plus_line = infile.readline()
                        plus_line = plus_line.strip()
                        quality_scores = infile.readline()
                        quality_scores = quality_scores.strip()

                        if not (
                            header and sequence and plus_line == "+" and quality_scores
                        ):
                            raise ValueError("Unexpected FASTQ format")

                        try:
                            timestamp = self.parse_timestamp(header)
                            timestamped_records.append(
                                (
                                    timestamp,
                                    (
                                        header.strip(),
                                        sequence,
                                        plus_line,
                                        quality_scores,
                                    ),
                                )
                            )
                        except ValueError as e:
                            logging.debug(f"Error parsing header: {e}")
        else:
            with open(self.input_fastq_path) as infile:
                while True:
                    header_line = infile.readline()
                    if not header_line:
                        break
                    if header_line.startswith("@"):
                        header = header_line.strip()
                        sequence = infile.readline().strip()
                        plus_line = infile.readline().strip()
                        quality_scores = infile.readline().strip()

                        if not (
                            header
                            and sequence
                            and plus_line == "+"
                            and len(sequence) == len(quality_scores)
                        ):
                            raise ValueError("Unexpected FASTQ format")

                        try:
                            timestamp = self.parse_timestamp(header)
                            timestamped_records.append(
                                (
                                    timestamp,
                                    (header, sequence, plus_line, quality_scores),
                                )
                            )
                        except ValueError as e:
                            logging.debug(f"Error parsing header: {e}")

        total_duration = (
            self.max_timestamp - self.min_timestamp
        ).total_seconds() / 3600
        files_created = 0
        max_intervals = max(int(total_duration / interval_hours) + 1, 1)
        # monitor how many files are created and how many records are in each file, and log this information
        logging.info(
            f"Total duration: {total_duration} hours, Max intervals: {max_intervals}, Total records: {len(timestamped_records)}"
        )
        # Create time-based mock files
        # interval_hours = min(interval_hours, total_duration / max_intervals)

        # interval_hours = min(interval_hours, total_duration / max_intervals)
        current_time = self.min_timestamp + timedelta(hours=start_time_point)
        stop_timestamp = self.min_timestamp + timedelta(hours=self.end_point)
        logging.info(
            f"Using interval hours: {interval_hours}, Stop timestamp: {stop_timestamp}"
        )
        while current_time < stop_timestamp:
            if stop_event and stop_event.is_set():
                # raise TerminateTaskGroup()
                raise asyncio.CancelledError("Mock file creation stopped by event")
                # break

            end_time = current_time + timedelta(hours=interval_hours)

            # Filter records for this time window
            window_records = [
                record
                for timestamp, record in timestamped_records
                if current_time <= timestamp < end_time
            ]

            # Limit number of records if needed
            if len(window_records) > max_records_per_file and max_records_per_file > 0:
                window_records = window_records[:max_records_per_file]

            if not window_records:
                logging.info(
                    f"No records found for time window {current_time} to {end_time}"
                )
                break

            # Create filename with time information
            time_str = current_time.strftime("%Y%m%d_%H%M%S")
            file_path = (
                Path(output_dir)
                / f"timepoint_{files_created}_{time_str}_{interval_hours}h.fastq.gz"
            )

            logging.info(
                f"Creating mock file for time window {current_time} to {end_time}..."
            )
            logging.info(f"  Records: {len(window_records)}")

            ## aiogzip does not support writing in text mode, so we need to encode the strings to bytes before writing
            # async with AsyncGzipTextFile(file_path, "wb",newline="") as f:
            #     for record in window_records:
            #         header, sequence, plus_line, quality_scores = record
            #         if not (header and sequence and plus_line == '+' and quality_scores):
            #             raise ValueError("Unexpected FASTQ format")
            #         await f.write(f"{header}\n{sequence}\n{plus_line}\n{quality_scores}\n".encode('utf-8'))
            # Write records to file

            with gzip.open(file_path, "wb") as f:
                for record in window_records:
                    header, sequence, plus_line, quality_scores = record
                    if not (
                        header and sequence and plus_line == "+" and quality_scores
                    ):
                        raise ValueError("Unexpected FASTQ format")

                    f.write(
                        f"{header}\n{sequence}\n{plus_line}\n{quality_scores}\n".encode(
                            "utf-8"
                        )
                    )

            logging.info(f"Created mock file: {file_path}")
            files_created += 1
            # Increment the current time by the interval for the next window

            # Move to next time window
            current_time = end_time

            # Small delay between file creation
            await asyncio.sleep(60)

        # Create additional random files if needed
        # while files_created < max_files:
        #     if stop_event and stop_event.is_set():
        #         break

        #     # Create random time window within the experiment duration
        #     start_fraction = random.uniform(0.0, 0.8)
        #     end_fraction = random.uniform(start_fraction + 0.1, 1.0)

        #     start_hours = start_fraction * total_duration
        #     end_hours = end_fraction * total_duration

        #     start_timestamp = self.min_timestamp + timedelta(hours=start_hours)
        #     end_timestamp = self.min_timestamp + timedelta(hours=end_hours)

        #     # Filter records for this random time window
        #     window_records = [
        #         record for timestamp, record in timestamped_records
        #         if start_timestamp <= timestamp <= end_timestamp
        #     ]

        #     # Limit number of records if needed
        #     if len(window_records) > max_records_per_file:
        #         window_records = window_records[:max_records_per_file]

        #     if not window_records:
        #         continue

        #     # Create filename with time information
        #     time_str = start_timestamp.strftime('%Y%m%d_%H%M%S')
        #     file_path = Path(output_dir) / f"random_{files_created}_{time_str}.fastq.gz"

        #     logging.info(f"Creating random mock file for time window {start_timestamp} to {end_timestamp}...")
        #     logging.info(f"  Records: {len(window_records)}")

        #     # Write records to file
        #     with gzip.open(file_path, "wb") as f:
        #         for record in window_records:
        #             header, sequence, plus_line, quality_scores = record
        #             if not (header and sequence and plus_line == '+' and quality_scores):
        #                 raise ValueError("Unexpected FASTQ format")

        #             f.write(f"{header}\n{sequence}\n+\n{quality_scores}\n".encode('utf-8'))

        #     logging.info(f"Created mock file: {file_path}")
        #     files_created += 1

        #     # Small delay between file creation
        #     await asyncio.sleep(0.1)

        # Verify all created files
        for file_path in Path(output_dir).iterdir():
            if file_path.suffix == ".fastq.gz":
                with gzip.open(file_path, "rb") as f:
                    content = f.read()
                    logging.info(f"File {file_path} size: {len(content)} bytes")
                    if len(content) == 0:
                        logging.error(f"File {file_path} is empty!")

                    # Verify FASTQ format
                    try:
                        lines = content.decode("utf-8").strip().split("\n")
                        if len(lines) % 4 != 0:
                            logging.error(f"File {file_path} has improper FASTQ format")
                        else:
                            for i in range(0, len(lines), 4):
                                header, seq, plus, qual = lines[i : i + 4]
                                if not (
                                    header.startswith("@")
                                    and plus == "+"
                                    and len(seq) == len(qual)
                                ):
                                    logging.error(
                                        f"File {file_path} has improper FASTQ format at record {i // 4}"
                                    )
                    except Exception as e:
                        logging.error(f"Error verifying file {file_path}: {e}")

        return files_created

    def extract_by_time_window(self, start_time, end_time, output_file=None):
        """
        Extract reads between specific datetime objects.
        """
        if self.min_timestamp is None:
            self.get_run_info()

        if start_time > end_time:
            raise ValueError("Start time must be before end time")

        extracted_reads = []

        if self.gzip_input:
            with gzip.open(self.input_fastq_path, "rt") as infile:
                for line in infile:
                    if line.startswith("@"):
                        header = line.strip()
                        try:
                            timestamp = self.parse_timestamp(header)
                            if start_time <= timestamp <= end_time:
                                sequence = next(infile).strip()
                                plus_line = next(infile).strip()
                                quality_scores = next(infile).strip()

                                if not (
                                    header
                                    and sequence
                                    and plus_line == "+"
                                    and quality_scores
                                ):
                                    raise ValueError("Unexpected FASTQ format")

                                extracted_reads.append(
                                    (header, sequence, plus_line, quality_scores)
                                )
                        except ValueError as e:
                            logging.debug(f"Error parsing header: {e}")
        else:
            with open(self.input_fastq_path) as infile:
                for line in infile:
                    if line.startswith("@"):
                        header = line.strip()
                        try:
                            timestamp = self.parse_timestamp(header)
                            if start_time <= timestamp <= end_time:
                                sequence = next(infile).strip()
                                plus_line = next(infile).strip()
                                quality_scores = next(infile).strip()

                                if not (
                                    header
                                    and sequence
                                    and plus_line == "+"
                                    and quality_scores
                                ):
                                    raise ValueError("Unexpected FASTQ format")

                                extracted_reads.append(
                                    (header, sequence, plus_line, quality_scores)
                                )
                        except ValueError as e:
                            logging.debug(f"Error parsing header: {e}")

        if output_file:
            with gzip.open(output_file, "wb") as outfile:
                for read in extracted_reads:
                    header, sequence, plus_line, quality_scores = read
                    outfile.write(
                        f"{header}\n{sequence}\n+\n{quality_scores}\n".encode("utf-8")
                    )

        return extracted_reads

    def get_time_windows(self, window_size_hours=1.0):
        """
        Get time windows for the experiment.
        """
        if self.min_timestamp is None:
            self.get_run_info()

        (self.max_timestamp - self.min_timestamp).total_seconds() / 3600
        windows = []

        current_time = self.min_timestamp
        window_num = 0

        while current_time < self.max_timestamp:
            end_time = current_time + timedelta(hours=window_size_hours)
            if end_time > self.max_timestamp:
                end_time = self.max_timestamp

            windows.append(
                {
                    "window_num": window_num,
                    "start_time": current_time,
                    "end_time": end_time,
                    "start_hours": (current_time - self.min_timestamp).total_seconds()
                    / 3600,
                    "end_hours": (end_time - self.min_timestamp).total_seconds() / 3600,
                }
            )

            current_time = end_time
            window_num += 1

        return windows

    def extract_time_window_to_file(
        self, window_num, window_size_hours=1.0, output_file=None
    ):
        """
        Extract a specific time window to a file.
        """
        windows = self.get_time_windows(window_size_hours)
        if window_num >= len(windows):
            raise ValueError(f"Window number {window_num} exceeds available windows")

        window = windows[window_num]
        return self.extract_by_time_range(
            start_hours=window["start_hours"],
            end_hours=window["end_hours"],
            output_file=output_file,
        )


# Command line interface
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="FASTQ Time Extractor")
    parser.add_argument("input_file", help="Input FASTQ file")
    parser.add_argument("--start", type=float, help="Start time in hours")
    parser.add_argument("--end", type=float, help="End time in hours")
    parser.add_argument("--output", help="Output file")
    parser.add_argument(
        "--create-mock", help="Create mock files in specified directory"
    )
    parser.add_argument(
        "--interval", type=float, default=1.0, help="Interval in hours for mock files"
    )

    args = parser.parse_args()

    extractor = FastQTimeExtractor(
        args.input_file, end_point=args.end if args.end is not None else None
    )

    if args.create_mock:
        # Run the async method in a new event loop
        async def run_create_mock():
            await extractor.create_mock_files(
                args.create_mock, interval_hours=args.interval
            )

        asyncio.run(run_create_mock())
    elif args.start is not None:
        extractor.extract_by_time_range(
            start_hours=args.start, end_hours=args.end, output_file=args.output
        )
    else:
        print("Run with --help for usage information")
