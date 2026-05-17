import gzip
import asyncio
import subprocess
import logging
from pathlib import Path
from datetime import datetime, timedelta
from typing import Any

import numpy as np
import random


def verify_with_blast(input_file, outfile, subject_db, format=6):
    """Run BLAST to verify the sequence against a reference database."""
    cmd = [
        "blastn",
        "-query",
        str(input_file),
        "-subject",
        subject_db,
        "-outfmt",
        str(format),
    ]

    result_path = Path(outfile)

    try:
        # Use asyncio subprocess for async execution
        with open(result_path, "w") as out_file:
            subprocess.run(
                cmd,
                stdout=out_file,
                check=True,
                stderr=subprocess.PIPE,
            )
        # with open(result_path, "w") as out_file:
        #     process = await asyncio.create_subprocess_exec(
        #         *cmd, stdout=out_file, stderr=asyncio.subprocess.PIPE
        #     )
        #     await process.communicate()

        logging.info(f"BLAST completed successfully. Results saved to {result_path}")

    except Exception as e:
        logging.error(f"Error running BLAST: {e}")


def compare_to_real_ref(real_ref, doras_ref):
    real_ref = Path(real_ref)
    doras_ref = Path(doras_ref)

    # Ensure the parent directory exists for the BAM file
    output_dir = doras_ref.parent / "mapping_verif_2"
    if not output_dir.exists():
        logging.info(f"Creating output directory: {output_dir}")
        output_dir.mkdir(parents=True, exist_ok=True)

    # Extract the base name of the fastq file without extension
    real_ref_base_name = real_ref.stem

    doras_ref_base_name = doras_ref.stem
    # Step 1: Run minimap2 and samtools view, sort
    # bam_file_path = output_dir / f"{fastq_base_name}.{fasta_base_name}.bam"
    try:
        logging.debug(
            f"Running minimap2 with fasta: {doras_ref} and the original reference: {real_ref}"
        )
        result = subprocess.run(
            ["minimap2", "-ax", "map-ont", str(doras_ref), str(real_ref)],
            check=True,
            stdout=subprocess.PIPE,
        )
    except subprocess.SubprocessError as r:
        logging.error(f"Failed to map the DORAS ref to the original reference {r}")

    result = subprocess.run(
        args=["samtools", "view", "-b", "-"],
        stdout=subprocess.PIPE,
        input=result.stdout,
    )
    result = subprocess.run(
        args=["samtools", "sort", "-"], input=result.stdout, stdout=subprocess.PIPE
    )

    coverage_real_ref = (
        output_dir / f"{real_ref_base_name}.{doras_ref_base_name}.quality.coverage"
    )
    try:
        logging.info(f"Saving the results of the coverage to {str(coverage_real_ref)}")
        subprocess.run(
            args=["samtools", "coverage", "-", "-o", coverage_real_ref],
            input=result.stdout,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        return coverage_real_ref
    except subprocess.CalledProcessError as r:
        logging.error(f"Error when writing the coverage results to disk {r}")


def run_mapping_pipeline(fasta_path, fastq_path, from_sr=True):

    # Convert input paths to Path objects
    fasta_path = Path(fasta_path)
    fastq_path = Path(fastq_path)

    # Ensure the parent directory exists for the BAM file
    output_dir = fasta_path.parent / "evaluation_recruitment_reads"
    if not output_dir.exists():
        logging.info(f"Creating output directory: {output_dir}")
        output_dir.mkdir(parents=True, exist_ok=True)

    # Extract the base name of the fastq file without extension
    fastq_base_name = fastq_path.stem

    fasta_base_name = fasta_path.stem
    if from_sr:  ## if you have reads that were already sequenced using as and you want to compare your sequence to the one that you created (e.g. during DORAS)
        fastq_output_path = output_dir / f"{fastq_base_name}.first500.fastq"
        try:
            with open(fastq_output_path, "w") as fastq_file:
                logging.debug("Extracting the first 500 reads")
                subprocess.run(
                    ["seqkit", "subseq", "-r", "1:500", str(fastq_path)],
                    stdout=fastq_file,
                    check=True,
                    stderr=subprocess.PIPE,
                )
        except subprocess.CalledProcessError as e:
            logging.error(f"Command failed: {' '.join(e.cmd)}")
            logging.error(f"Error output: {e.stderr.decode('utf-8')}")
            return None
    else:
        # Step 1: Run minimap2 and samtools view, sort
        bam_file_path = output_dir / f"{fastq_base_name}.{fasta_base_name}.bam"
        try:
            logging.debug(
                f"Running minimap2 with fasta: {fasta_path} and fastq: {fastq_path}"
            )
            result = subprocess.run(
                ["minimap2", "-ax", "map-ont", str(fasta_path), str(fastq_path)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )

            logging.debug("Converting to BAM and sorting")
            with open(bam_file_path, "wb") as bam_file:
                subprocess.run(
                    ["samtools", "view", "-b", "-F", "4"],
                    input=result.stdout,
                    stdout=bam_file,
                    check=True,
                    stderr=subprocess.PIPE,
                )

            sorted_bam_file_path = (
                output_dir / f"{fastq_base_name}.{fasta_base_name}.sorted.bam"
            )
            subprocess.run(
                ["samtools", "sort", "-o", str(sorted_bam_file_path), bam_file_path],
                check=True,
                stderr=subprocess.PIPE,
            )

            # Index the sorted BAM file
            subprocess.run(
                ["samtools", "index", str(sorted_bam_file_path)],
                check=True,
                stderr=subprocess.PIPE,
            )

            # Output the coverage
            output_coverage = (
                output_dir / f"{fastq_base_name}.{fasta_base_name}.coverage"
            )
            with open(output_coverage, "w") as cov:
                subprocess.run(
                    ["samtools", "coverage", str(sorted_bam_file_path)],
                    check=True,
                    stdout=cov,
                )
        except subprocess.CalledProcessError as e:
            logging.error(f"Command failed: {' '.join(e.cmd)}")
            logging.error(f"Error output: {e.stderr.decode('utf-8')}")
            return None

        # Step 2: Convert BAM to FASTQ and extract the first 400 reads
        fastq_output_path = output_dir / f"{fastq_base_name}.first500.fastq"
        try:
            logging.debug(f"Converting BAM to FASTQ: {sorted_bam_file_path}")
            result = subprocess.run(
                ["samtools", "fastq", str(sorted_bam_file_path)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )

            with open(fastq_output_path, "w") as fastq_file:
                logging.debug("Extracting the first 500 reads")
                subprocess.run(
                    ["seqkit", "subseq", "-r", "1:500"],
                    input=result.stdout,
                    stdout=fastq_file,
                    check=True,
                    stderr=subprocess.PIPE,
                )
        except subprocess.CalledProcessError as e:
            logging.error(f"Command failed: {' '.join(e.cmd)}")
            logging.error(f"Error output: {e.stderr.decode('utf-8')}")
            return None

        # Step 3: Run minimap2 and samtools view, sort, index
        genes_bam_file_path = (
            output_dir / f"{fastq_base_name}.mlst_genes_alone.sorted.bam"
        )
    try:
        logging.debug(f"Running minimap2 on FASTQ extracted reads: {fastq_output_path}")
        result = subprocess.run(
            [
                "minimap2",
                "-ax",
                "map-ont",
                str(fasta_path),
                str(fastq_output_path),
            ],  # map-ont is the readfish standard and sr is the preset used by nanopore
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        logging.debug("Converting to BAM and sorting")
        result = subprocess.run(
            [
                "samtools",
                "view",
                "-b",
            ],
            input=result.stdout,
            stdout=subprocess.PIPE,
            check=True,
            # stderr=subprocess.PIPE
        )

        # Sort the BAM file
        sorted_genes_bam_file_path = (
            output_dir / f"{fastq_base_name}.{fasta_base_name}.sorted.bam"
        )
        subprocess.run(
            ["samtools", "sort", "-o", str(sorted_genes_bam_file_path)],
            check=True,
            input=result.stdout,
            stderr=subprocess.PIPE,
        )

        # Index the sorted BAM file
        subprocess.run(
            ["samtools", "index", str(sorted_genes_bam_file_path)],
            check=True,
            stderr=subprocess.PIPE,
        )
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {' '.join(e.cmd)}")
        logging.error(f"Error output: {e.stderr.decode('utf-8')}")
        return None

    # Get flagstats
    command_5 = ["samtools", "flagstats", str(sorted_genes_bam_file_path)]
    try:
        logging.debug("Running samtools flagstats to get mapping statistics")
        result = subprocess.run(
            command_5,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        logging.debug(msg=f"Flagstats : {result}")
        total_reads_line = result.stdout.splitlines()[0]
        total_reads = int(total_reads_line.split()[0])
        mapped_reads_line = result.stdout.splitlines()[6]
        mapped_reads = int(mapped_reads_line.split()[0])
        prim_mapped_reads_line = result.stdout.splitlines()[7]
        prim_mapped_reads = int(prim_mapped_reads_line.split()[0])
        if total_reads > 0:
            percentage_primary_mapped = (prim_mapped_reads / total_reads) * 100
            percentage_mapped = (mapped_reads / total_reads) * 100
            logging.info(f"mapped_reads: {mapped_reads}, total_reads: {total_reads}")
        else:
            logging.warning(
                f"Total mapped reads is zero, cannot calculate percentage {fasta_base_name}."
            )
            percentage_primary_mapped = 0.0
            percentage_mapped = 0.0

        logging.info(
            f"percentage mapped {percentage_mapped:.2f}, total_reads: {total_reads}"
        )
        command_6 = ["samtools", "coverage", str(sorted_genes_bam_file_path)]
        try:
            logging.info("Running samtools coverage to get coverage statistics")
            result_cov = subprocess.run(
                command_6,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            logging.info(msg=f"Coverage : {result_cov.stdout}")
        except subprocess.CalledProcessError as e:
            logging.error(f"Command failed: {' '.join(e.cmd)}")
            logging.error(f"Error output: {e.stderr.decode('utf-8')}")
            # return None

        return f"{percentage_mapped:.2f}\t{percentage_primary_mapped:.2f}\t{total_reads}\t{fasta_base_name}"
        # for line in result.stdout.splitlines():

        #     if "primary" in line:
        #         parts = line.split()
        #         mapped_reads = int(parts[0])
        #         total_reads = int(parts[2])
        #         logging.debug(f"mapped_reads: {mapped_reads}, total_reads: {total_reads}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {' '.join(e.cmd)}")
        logging.error(f"Error output: {e.stderr.decode('utf-8')}")
        return None


## functions to be used in case


def run_minimap2_map_ont(reference_fasta, fastq_file, output_bam):
    cmd = ["minimap2", "-x", "map-ont", reference_fasta, fastq_file, "-o", output_bam]
    subprocess.run(cmd, check=True)


def sort_samtools(input_bam, sorted_bam):
    cmd = ["samtools", "sort", input_bam, "-o", sorted_bam]
    subprocess.run(cmd, check=True)


def evaluate_coverage(sorted_bam, output_prefix):
    depth_file = f"{output_prefix}_coverage.txt"
    cmd = ["samtools", "depth", sorted_bam, "-o", depth_file]
    subprocess.run(cmd, check=True)
    return depth_file


# Example usage
def polish_with_medaka(input_file, reference_fasta, output_dir, bacteria):
    """Run medaka_consensus with proper error handling.

    Args:
        input_file: Path to input reads file
        reference_fasta: Path to reference FASTA file
        output_dir: Output directory path
        bacteria: Boolean flag for bacterial analysis mode

    Raises:
        subprocess.CalledProcessError: If medaka command fails
        ValueError: For invalid inputs or incompatible data
    """

    # Input validation before running the command
    if not Path(input_file).exists():
        raise FileNotFoundError(f"Input file {input_file} does not exist")

    if not Path(reference_fasta).exists():
        raise FileNotFoundError(f"Reference file {reference_fasta} does not exist")

    # Define the medaka_consensus command
    cmd_base = [
        "medaka_consensus",
        "-x",
        "-b",
        "50",
        "-f",
        "-i",
        str(input_file),
        "-d",
        str(reference_fasta),
        "-o",
        str(output_dir),
    ]

    # Add bacteria flag if specified
    cmd = cmd_base + ["--bacteria"] if bacteria else cmd_base

    try:
        logging.debug(f"Running medaka command: {' '.join(cmd)}")

        # Use check=True to raise CalledProcessError on failure
        subprocess.run(cmd, check=True, stderr=subprocess.DEVNULL)

        logging.info(
            f"Medaka polishing completed successfully. Output saved to {output_dir}"
        )

    except subprocess.CalledProcessError:
        error_msg = "ERROR: --bacteria was specified but input data is not compatible. If you wish to proceed anyway, please provide a bacterial model explicitly using -m."

        logging.error(error_msg)

        # Raise the specific exception with your exact message
        raise ValueError(error_msg) from None
    except FileNotFoundError as e:
        # This handles case where medaka executable isn't found
        error_msg = f"Medaka executable not found: {e}"
        logging.error(error_msg)
        raise RuntimeError(error_msg) from e

    except Exception as e:
        # Catch any other unexpected errors
        error_msg = f"Unexpected error running Medaka: {e}"
        logging.error(error_msg)
        raise RuntimeError(error_msg) from e


async def create_mock_files_by_timepoints(
    stop_event,
    existing_fastq_path,
    testdir: Path,
    interval_hours=1.0,
    max_files=24,
    max_records_per_file=1000,
    max_files_to_create=30,
):
    """
    Create mock FASTQ files organized by time points (e.g., hourly intervals).
    Ensures output is compatible with FastQTimeExtractor.

    Args:
        stop_event: Event to signal when to stop creating files
        existing_fastq_path: Path to existing FASTQ file to use as source
        testdir: Directory to create mock files in
        interval_hours: Time interval in hours for each file
        max_files: Maximum number of time-based files to create
        max_records_per_file: Maximum records per output file
        max_files_to_create: Maximum total files to create
    """
    # Read and parse the source FASTQ file
    with gzip.open(str(existing_fastq_path), "rb") as f:
        fastq_file = f.read().decode("utf-8")

    # Parse FASTQ records and extract timestamps
    splitted = fastq_file.split("\n")
    records = []
    for i in range(0, len(splitted), 4):
        record = splitted[i : i + 4]
        if len(record) == 4:
            records.append(record)

    if not records:
        logging.error("No valid FastQ records found in the input file.")
        return

    # Extract timestamps from records (same method as FastQTimeExtractor)
    timestamped_records = []
    for record in records:
        header = record[0]
        try:
            # Parse timestamp from header (ONT format) - same as FastQTimeExtractor
            timestamp = None
            for part in header.split(" "):
                if "start_time=" in part:
                    timestamp = datetime.fromisoformat(
                        part.replace("start_time=", "").replace("+01:00", "")
                    )
                    break

            if timestamp:
                timestamped_records.append((timestamp, record))
        except ValueError as e:
            logging.debug(f"Error parsing timestamp: {e}")

    if not timestamped_records:
        logging.error("No valid timestamps found in the input file.")
        return

    # Sort records by timestamp
    timestamped_records.sort(key=lambda x: x[0])

    # Determine time range
    min_time = timestamped_records[0][0]
    max_time = timestamped_records[-1][0]
    total_duration = (max_time - min_time).total_seconds() / 3600  # in hours

    logging.info(
        f"Source file time range: {min_time} to {max_time} ({total_duration:.2f} hours)"
    )

    # Create test directory if it doesn't exist
    test_folder = Path(testdir)
    if not test_folder.exists():
        test_folder.mkdir(parents=True)

    total_files_created = 0
    current_time = min_time

    # Create files for each time interval
    while total_files_created < max_files and total_files_created < max_files_to_create:
        # Calculate time range for this file
        end_time = current_time + timedelta(hours=interval_hours)

        # Filter records for this time window
        window_records = [
            record
            for timestamp, record in timestamped_records
            if current_time <= timestamp < end_time
        ]

        # Limit number of records if needed
        if len(window_records) > max_records_per_file:
            window_records = window_records[:max_records_per_file]

        if not window_records:
            logging.info(
                f"No records found for time window {current_time} to {end_time}"
            )
            break

        # Create filename with time information
        time_str = current_time.strftime("%Y%m%d_%H%M%S")
        file_path = (
            test_folder
            / f"timepoint_{total_files_created}_{time_str}_{interval_hours}h.fastq.gz"
        )

        logging.info(
            f"Creating mock file for time window {current_time} to {end_time}..."
        )
        logging.info(f"  Records: {len(window_records)}")

        # Write records to file - ensure proper FASTQ format
        with gzip.open(file_path, "wb") as f:
            for record in window_records:
                # Ensure each record has proper FASTQ format
                header, sequence, plus_line, quality_scores = record
                if not (header and sequence and plus_line == "+" and quality_scores):
                    raise ValueError("Unexpected FASTQ format")

                # Write in proper FASTQ format (same as extractor)
                f.write(f"{header}\n{sequence}\n+\n{quality_scores}\n".encode("utf-8"))

        logging.info(f"Created mock file: {file_path}")
        total_files_created += 1

        # Move to next time window
        current_time = end_time

        # Check if we should stop
        if stop_event.is_set():
            break

        await asyncio.sleep(0.1)  # Small delay between file creation

    # Create additional random files if needed
    while total_files_created < max_files_to_create:
        if stop_event.is_set():
            break

        # Create random time window within the experiment duration
        start_fraction = random.uniform(0.0, 0.8)  # Start somewhere in first 80%
        end_fraction = random.uniform(start_fraction + 0.1, 1.0)  # End after start

        start_hours = start_fraction * total_duration
        end_hours = end_fraction * total_duration

        start_timestamp = min_time + timedelta(hours=start_hours)
        end_timestamp = min_time + timedelta(hours=end_hours)

        # Filter records for this random time window
        window_records = [
            record
            for timestamp, record in timestamped_records
            if start_timestamp <= timestamp <= end_timestamp
        ]

        # Limit number of records if needed
        if len(window_records) > max_records_per_file:
            window_records = window_records[:max_records_per_file]

        if not window_records:
            continue

        # Create filename with time information
        time_str = start_timestamp.strftime("%Y%m%d_%H%M%S")
        file_path = test_folder / f"random_{total_files_created}_{time_str}.fastq.gz"

        logging.info(
            f"Creating random mock file for time window {start_timestamp} to {end_timestamp}..."
        )
        logging.info(f"  Records: {len(window_records)}")

        # Write records to file
        with gzip.open(file_path, "wb") as f:
            for record in window_records:
                header, sequence, plus_line, quality_scores = record
                if not (header and sequence and plus_line == "+" and quality_scores):
                    raise ValueError("Unexpected FASTQ format")

                f.write(f"{header}\n{sequence}\n+\n{quality_scores}\n".encode("utf-8"))

        logging.info(f"Created mock file: {file_path}")
        total_files_created += 1

        # Small delay between file creation
        await asyncio.sleep(0.1)

        if stop_event.is_set():
            break

    # Verify all created files are not empty and have proper format
    for file_path in test_folder.iterdir():
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


# Mock creation function


async def create_mock_files(
    stop_event,
    existing_fastq_path,
    testdir: Path,
    interval=2,
    num_files=10,
    max_records_per_file=1000,
    max_files_to_create=30,
):
    with gzip.open(str(existing_fastq_path), "rb") as f:
        fastq_file = f.read().decode("utf-8")

    splitted = fastq_file.split("\n")
    records = []
    for i in range(0, len(splitted), 4):
        record = splitted[i : i + 4]
        if len(record) == 4:
            records.append(record)

    if not records:
        logging.error("No valid FastQ records found in the input file.")
        return

    test_folder = Path(testdir)
    if not test_folder.exists():
        test_folder.mkdir(parents=True)
    total_files = 0
    while True:
        if total_files >= max_files_to_create:
            break
        for i in range(num_files):
            total_files += 1
            random.shuffle(records)
            num_records_to_include = min(max_records_per_file, len(records))
            selected_records = [records[j] for j in range(num_records_to_include)]

            file_path = (
                test_folder
                / f"test_{i}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.fastq.gz"
            )

            logging.info(f"Creating mock file: {file_path}...")
            with gzip.open(file_path, "wb") as f:
                for record in selected_records:
                    fastq_content = "\n".join(record) + "\n"
                    f.write(fastq_content.encode("utf-8"))

            logging.info(f"Created mock file: {file_path}")

        await asyncio.sleep(interval)

        if stop_event.is_set():
            await stop_event.wait()
        else:
            stop_event.clear

    # Check if the files are not empty
    for file_path in test_folder.iterdir():
        if file_path.suffix == ".fastq.gz":
            with gzip.open(file_path, "rb") as f:
                content = f.read()
                logging.info(f"File {file_path} size: {len(content)} bytes")
                if len(content) == 0:
                    logging.error(f"File {file_path} is empty!")


def rev_com(seq: str):
    seq_inv = seq[::-1]
    bases_dict: dict[str, str] = {
        "A": "T",
        "C": "G",
        "T": "A",
        "G": "C",
        "N": "N",
        "a": "t",
        "c": "g",
        "t": "a",
        "g": "c",
        "n": "n",
    }
    return "".join([bases_dict.get(nuc, nuc) for nuc in list(seq_inv)])


def get_depth(bam_file: str, reference_fasta: str) -> dict:
    """
    Calculate the depth of coverage for each position in the reference sequences.
    """
    import pysam
    from .io import FastaHolder

    samfile = pysam.AlignmentFile(bam_file, "rb")
    fasta_dict = FastaHolder(reference_fasta).get_dict_fasta()

    depth_dict = {}
    for record_id, sequence in fasta_dict.items():
        seq_length = len(sequence)
        depth_array = np.zeros(seq_length)
        for pileupcolumn in samfile.pileup(record_id, 0, seq_length):
            depth_array[pileupcolumn.pos] = pileupcolumn.nsegments
        depth_dict[record_id] = depth_array

    samfile.close()
    return depth_dict


def verify_depth(depth_dict: dict, threshold: int) -> tuple:
    """
    Verify if the depth of coverage is above a certain threshold and return the range.
    """
    # Assuming single record for simplicity as per existing test structure
    for record_id, depth_array in depth_dict.items():
        above_threshold = np.where(depth_array >= threshold)[0]
        if len(above_threshold) == 0:
            return (0, 0)
        return (int(above_threshold[0]), int(above_threshold[-1]) + 1)
    return (0, 0)


# PAF and alignment utility functions from pseudoasm
def clean_paf(paf, id_to_int=True):
    """
    This function will change the col names and mutate
    the read ids from strings to integers it takes a df pandas ans returns a df pandas
    #TODO clean this functions out, rething the way the strand are passed to the down functions
    """
    paf_subset = paf.iloc[:, 0:13]

    paf_subset.columns = [
        "query_id",
        "query_len",
        "query_start",
        "query_end",
        "strand",
        "target_id",
        "target_len",
        "target_start",
        "target_end",
        "aln_len",
        "match",
        "mapping_quality",
        "cigar",
    ]
    ids = paf_subset["query_id"].to_numpy()
    # return the strand of each query id

    strand = paf_subset["strand"].to_numpy()
    ids_strand = {id: strand[i] for i, id in enumerate(ids)}
    ids = np.append(ids, paf_subset["target_id"].to_numpy())
    # logging.debug(f"query id and target ids for debug{ids}")
    ids = np.unique(ids)
    ## dict to get the int corresponding to id
    id_to_int = {id: int for int, id in enumerate(ids)}
    ## adding the int_to_id for future purposes
    paf_subset["strand"] = paf_subset["strand"].str.strip()
    if id_to_int:
        int_to_id = {int: str.strip(id) for int, id in enumerate(ids)}
        ## replaced the ids by an int to improbve performance
        paf_subset["target_id"] = paf_subset["target_id"].map(id_to_int)
        paf_subset["query_id"] = paf_subset["query_id"].map(id_to_int)
        # remove white space in strand col
        return paf_subset, int_to_id, ids_strand
    else:
        return paf_subset


def index_rev_reads(paf, mode="query") -> Any:
    """
    Reverse reads have different index
    this functions aims to change the read as if it were positive strand to make it
    easier to assemble the final sequence
    modified to enable a mode that choses which of the target of the query to modify default query (ASMLST)
    """

    strand = paf["strand"].to_numpy()
    start_seq_tar = paf[f"{mode}_start"].to_numpy()
    end_seq_tar = paf[f"{mode}_end"].to_numpy()
    # end_seq_q = paf["query_end"].to_numpy()
    len_seq_tar = paf[f"{mode}_len"].to_numpy()
    strand_mask_neg = strand == "-"
    # print(strand)
    strand_mask_pos = strand == "+"
    start_correct = len_seq_tar - end_seq_tar
    end_correct = len_seq_tar - start_seq_tar
    # print(start_correct)
    # print(strand_mask_neg)
    paf[f"{mode}_start"] = (
        start_seq_tar * strand_mask_pos + start_correct * strand_mask_neg
    )  ## we take the original start for the positive strand and the modified ones for the rev strand
    paf[f"{mode}_end"] = end_seq_tar * strand_mask_pos + end_correct * strand_mask_neg
    return paf
