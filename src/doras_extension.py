import asyncio
import aiofiles
import threading
import pysam

import logging
import os
import time
from pathlib import Path
from typing import (
    Any,
    Dict,
    List,
    Tuple,
    Union,
)

import numpy as np
import pandas as pd
import mappy as mp

from .config import Phase
from .config import DorasParams

from src.bigsdb_tools import (
    query_st_database,
    read_encode_fasta_file,
)
import src.mapping as mpaf
from src.io import FastaHolder, FastqHolder, save_list_nodes_to_file

from src.utils import polish_with_medaka
from src.utils import rev_com, clean_paf, index_rev_reads
from src.simulation_model import find_optimal_quantile

# Setup logging configuration
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def compute_max_extensions(
    paf,
    min_ext=0,
    max_ext=3000,
    min_aln=80,
    overhang_ratio=0.8,
    mapq_filter=0,
    up_trim=0,
    down_trim=0,
):
    """
    Compute the max extension for the upstream and downstream flanking regions of a PAF alignment.

    Parameters:
    paf (pd.DataFrame): The PAF alignment DataFrame after cleaning.

    Returns:
    dict: A dictionary where keys are target IDs, and values are dictionaries containing the maximum
          extensions in the upstream ('up') and downstream ('down') directions. Each sub-dictionary
          contains the query ID and corresponding extension coordinates.

    Notes:
    - The function removes reads that are entirely within the target region.
    - It calculates the upstream (`up_ext`) and downstream (`down_ext`) extensions for each alignment.
    - For each target, it identifies the alignments with the maximum upstream and downstream extensions.
    - The output includes query IDs and relevant extension coordinates.

    Example:
    {
        "target1": {
            "up": ["query1", 0, 3],
            "down": ["query2", 5, 8]
        },
        "target2": {
            "up": ["query3", 2, 4],
            "down": ["query4", 7, 10]
        }
    }
    """

    paf["down_ext"] = paf["query_len"] - (
        paf["query_end"] + (paf["target_len"] - paf["target_end"])
    )  # Updated last version was wrong

    # Calculate upstream extension
    paf["up_ext"] = paf["query_start"] - paf["target_start"]
    # Ratio of the down_ext and up_ext vs the aln length inspired by the algorithm used in miniasm (Li. et al 2016)
    paf["ratio_aln_down_ext"] = paf["down_ext"] / paf["aln_len"]
    paf["ratio_aln_up_ext"] = paf["up_ext"] / paf["aln_len"]

    paf["down_ext_trim"] = (
        down_trim
        # setting the trimming using the parameters of the functions in case you want to do basic trimming you can
    )
    paf["up_ext_trim"] = up_trim

    paf.loc[paf["down_ext"] <= min_ext, "down_ext"] = 0

    paf.loc[paf["up_ext"] <= min_ext, "up_ext"] = 0

    paf = paf.dropna(subset=["up_ext", "down_ext"])

    # filter out all low quality mapping
    if mapq_filter > 0:
        paf = paf.loc[paf["mapping_quality"] >= mapq_filter]

    # paf = paf.loc[(paf["down_ext"] >= min_ext) | (paf["up_ext"] >= min_ext) ]
    # TODO separate these two parts into two functions for clarity
    coord_targets = {}
    for target in paf["target_id"].unique():
        # Find the index of the maximum downstream extension for the current target
        paf_target = paf[paf["target_id"] == target]

        paf_target = paf_target.reset_index(drop=True)

        try:
            # branch aln_len here we filter out aln that contain a long overhang'
            # an overhang consists of the unaligned region situated at the both ends of an alignment
            # a too long overhang (when compared to the aln length) suggests a bad alignment and thus the read
            # is not kept

            # paf_target["target_len"]*min_aln_ext_ratio
            paf_target["overhang"] = np.minimum(
                paf_target["target_start"], paf_target["query_start"]
            ) + np.minimum(
                paf_target["query_len"] - paf_target["query_end"],
                paf_target["target_len"] - paf_target["target_end"],
            )
            threshold = np.minimum(
                paf_target["aln_len"] * overhang_ratio,
                np.minimum(paf_target["target_len"] * overhang_ratio, 1000),
            )
            paf_target = paf_target.loc[paf_target["overhang"] <= threshold]
        except Exception as r:
            logging.error(f"Problem when filtering the overhang length{r}")

        paf_target = paf_target.reset_index(drop=True)
        down_index_max = int(
            np.argmax(paf_target["down_ext"])
        )  # TODO Urgent, change that part to accept to take the second longest extension if it matches the criteria, Why? because now if the longest is not acceptable
        # Find the index of the maximum upstream extension for the current target !warning the data is filterd by target so it is not the same looking at the entire paf  (Debug)
        up_index_max = int(np.argmax(paf_target["up_ext"]))

        # Get the query ID, start position (0), and upstream next coordinate for the max upstream extension
        # up_ext_coord = paf.loc[up_index_max, ["query_id", 0, "up_ext"]].tolist()
        # verify the strand if + do that:

        up_ext_coord = (
            paf_target.loc[up_index_max, "query_id"],
            int(paf_target.loc[up_index_max, "up_ext_trim"]),
            int(
                paf_target.loc[up_index_max, "up_ext"]
                # + paf_target.loc[up_index_max, "up_ext_trim"] #Removed because no needd
            ),
        )
        down_ext_coord = (
            paf_target.loc[down_index_max, "query_id"],
            int(
                paf_target.loc[down_index_max, "query_len"]
                - paf_target.loc[down_index_max, "down_ext"]
            ),
            int(
                paf_target.loc[down_index_max, "query_len"]
                - paf_target.loc[down_index_max, "down_ext_trim"]
            ),
        )

        logging.debug(
            f"The coordinates for target {target} has the largest upstream ext is {up_ext_coord} and the downstream {down_ext_coord}"
        )
        # Get the query ID, downstream next coordinate, and query length for the max downstream extension
        # down_ext_coord = paf.iloc[down_index_max, ["query_id", "down_ext", "query_len"]].tolist()
        # the extension lengths can be calculated : up_ext = dict_ext_max["up"][2], downext_length = dict_ext_max["down"][2]- dict_ext_max["down"][1]
        coords = {"up": up_ext_coord, "down": down_ext_coord}

        coord_targets[target] = coords

    return coord_targets


def get_extensions_by_mapping(
    query_base_ref, extended_ref
):  # -> dict[Any, dict[str, type[tuple]]]:
    "The extension sizes will be inferred by mapping"
    "the base references onto the extended reference"
    "the paf will be used where the query id will contain the base ref with the loci"
    ":params: query_base_ref the sequences of the initial reference that were extended"
    ":params: the reference being extended that we need to evaluate the position of the genes from"
    ":return: dict with bed-like coord of the extensions to the seq"
    ""
    aligner = mp.Aligner(extended_ref, preset="map-ont")  # load reference

    if not aligner:
        raise Exception("ERROR: failed to load/build index")

    extensions_positions = {
        name: {"up": None, "down": None} for name, _, _ in mp.fastx_read(query_base_ref)
    }
    logging.debug(
        f"Verification of extensions by  Performing alignment {query_base_ref} to {extended_ref}"
    )
    for name, seq, qual in mp.fastx_read(query_base_ref):
        # print(seq)
        prev_max = 0
        for hit in aligner.map(seq):  # traverse alignments
            if (
                hit.mlen > prev_max
            ):  # verify that we do not add small random alignment but in a real case the id should match as well
                up_ext = (
                    0,
                    hit.r_st,
                    hit.r_st,
                )  # tupe with the coord of the up ext and the length thereof
                down_ext = (
                    hit.r_en,
                    hit.ctg_len,
                    hit.ctg_len - hit.r_en,
                )  # tuple with the coord of the down ext and the length thereof
                prev_max = hit.mlen
                ## Add the current extension to the dict
                extensions_positions.get(name)["up"] = up_ext
                extensions_positions.get(name)["down"] = down_ext
                # extensions_infos[name] = up_down
    return extensions_positions


def trim_sequences(
    extension_positions: Any,
    fasta_dict: Any,
    q50_extension: int,
    target_size: int,
    barcode: str,
):
    """
    The final sequences that reached the ideal size are trimmed on each end to
    match the self.target size in case they overextended (common case). The calculation is
    based on the position of the mlst genes in the sequence
    gene    start   end
    adk_up     0    43000
    adk_down    45000 70000
    target_size = 15000 (upstream and downstream)
    return adk[43000-target_size:45000+target_size]

    """
    trimmed_fastas = {}
    for header in fasta_dict:
        up_pos, down_pos = (
            extension_positions.get(header)["up"][1],
            extension_positions.get(header)["down"][0],
        )
        seq = fasta_dict.get(header)
        up_target_size = max(up_pos - target_size, 0)
        up_q50_size = max(up_pos - q50_extension, 0)
        up_cut_pos = min(up_q50_size, up_target_size)
        down_target_size = min(down_pos + target_size, len(seq))
        down_q50_size = min(down_pos + q50_extension, len(seq))
        down_cut_pos = max(down_q50_size, down_target_size)
        trimmed_fastas[f"{header}"] = seq[
            up_cut_pos:down_cut_pos
        ]  # TODO Modiifcation to prevent taking a seq with minus coordinates and allows to force the trimming to finalize

        # trimmed_fastas[f"{header}"] = fasta_dict.get(header)[
        #     max(up_pos - target_size,0) : min(down_pos + target_size,len(seq)) #TODO Modiifcation to prevent taking a seq with minus coordinates and allows to force the trimming to finalize
        # ]
    return trimmed_fastas


async def process_paf_and_compute_extensions(
    paf: pd.DataFrame,
    mapq_filter=10,
    up_trim=0,
    down_trim=0,
    min_aln_ext_ratio=0.8,
) -> Dict[str, List[Tuple[int, int]]]:
    """
       Wrapper function, process the PAF DataFrame to clean alignments and compute maximum extensions.

    , min_ext=0, max_ext = 3000,min_aln=80, min_aln_ext_ratio=1.2, mapq_filter = 0, up_trim = 0, down_trim = 0):
        :param paf: DataFrame containing alignment information.
        :param int_to_id: Dictionary mapping integer IDs to sequence identifiers.
        :param ref_file: Path to the reference file.
        :return: Dictionary of sequences with their max extension coordinates.
    """
    try:
        # Clean PAF
        paf_cleaned, int_to_id, ids_strand = clean_paf(paf)
        logging.debug("Cleaned PAF loaded successfully")

        # Index reverse reads
        paf_indexed = index_rev_reads(paf_cleaned)
        # logging.debug(f"paf after index_rev_reads:\n{paf_indexed}")
        # Compute maximum extensions
        dict_max_ext = compute_max_extensions(
            paf_indexed,
            mapq_filter=mapq_filter,
            overhang_ratio=min_aln_ext_ratio,
            up_trim=up_trim,
            down_trim=down_trim,
        )
        if not dict_max_ext:
            logging.warning(f"No Max extensions found {dict_max_ext}")
        else:
            logging.debug(f"Number of extensions found : {len(dict_max_ext)}")
        return dict_max_ext, int_to_id, ids_strand

    except Exception as e:
        logging.error(f"Processing PAF failed: {e}")
        return None


# function to
def build_ref_from_ext_coord(
    coords: dict, fasta_dict: dict, int_to_id: dict, strand_ids
):
    """
    The function will take the coordinate (output of compute_max_extensions) and construct a new reference sequence
    sequence from the query sequences and target sequences.
    Negative strand reads should be reverse complemented !!!
    e.g up_ext+target+down_ext
    Input ex (one isolate i.e barcode):
    {
        "target1": {
            "up": ["query1", 0, 3],
            "down": ["query2", 5, 8]
        },
        "target2": {
            "up": ["query3", 2, 4],
            "down": ["query4", 7, 10]
        }
    }
    Output ex:
    {
        "target1": ["ATCT"],
        "target2": ["ATATCT"]
    }

    """
    new_ref_seq = (
        [""] * 3
    )  # upstream extension at pos 0 , current ref at 1, and downstream ext at postiion 2
    target_ref_seqs = {}
    # TODO  should I transform all the int to id instead of doing a transformation each time
    for coord in coords:  # TODO why not use .items()
        target_coord = coords[coord]
        q_id, up_ext_coord1, up_ext_coord2 = target_coord["up"]
        up_seq, down_seq = "", ""
        # add up extension
        id_up_seq = int_to_id[q_id]
        if up_ext_coord1 >= 0 and up_ext_coord2 > 0:
            up_seq = fasta_dict[id_up_seq]
            if strand_ids[id_up_seq] == "-":
                up_seq = rev_com(up_seq)
            new_ref_seq[0] = up_seq[up_ext_coord1:up_ext_coord2]

        # add the reference sequence either the base sequence or the previously (growing ) pseudoreference
        else:  # setting to 0 if the values are eventually negative for the scoring
            up_ext_coord1, up_ext_coord2 = 0, 0
        new_ref_seq[1] = fasta_dict[int_to_id[coord]]

        q_id, down_ext_coord1, down_ext_coord2 = target_coord["down"]

        # add down extension if it exists
        id_down_seq = int_to_id[q_id]
        if (
            down_ext_coord1 > 0
            and down_ext_coord2 > 0
            and down_ext_coord1 != down_ext_coord2
        ):
            down_seq = fasta_dict[id_down_seq]
            if strand_ids[id_down_seq] == "-":
                down_seq = rev_com(down_seq)
            new_ref_seq[2] = down_seq[down_ext_coord1:down_ext_coord2]
        else:  # set to 0 to prevent to decrease the score
            down_ext_coord1, down_ext_coord2 = 0, 0
        if new_ref_seq[0] != "" or new_ref_seq[2] != "":
            target_ref_seqs[int_to_id[coord]] = "".join(new_ref_seq)
            logging.debug(
                f"Read ids used to extend:\n upstream: {id_up_seq} at pos {target_coord['up']} \n downstream {id_down_seq} at pos {target_coord['down']}"
            )
        new_ref_seq = [""] * 3

    return target_ref_seqs


def save_dict_to_fasta(data_dict, output_file):
    """
    Write a dictionary of sequences to a FASTA file.

    :param data_dict: Dictionary where keys are sequence identifiers and values are sequences.
    :param output_file: Path to the output FASTA file.
    """
    with open(output_file, "w") as fasta_handle:
        for identifier, sequence in data_dict.items():
            logging.debug(msg=f"the ref {identifier} is of length {len(sequence)}")
            fasta_handle.write(f">{identifier}\n")
            fasta_handle.write(f"{sequence}\n")


# Example usage within your class or function
def construct_reference_from_extensions(
    dict_max_ext: dict,
    ref_dict: dict,
    int_to_id: dict,
    strand_ids: dict,
) -> bool:
    """
    Construct and build the new reference sequences based on max extensions.

    :param dict_max_ext: Dictionary of maximum extension coordinates.
    :param ref_dict: Dictionary containing reference sequences.
    :param int_to_id: Dictionary mapping integer IDs to sequence identifiers.
    :param strand_ids: Dictionary containing strand information.
    :param output_path: Path where the new FASTA file will be saved.
    :return: True if successful, False otherwise.
    """
    try:
        # logging.info(f"paf after clean_paf:\n{paf.head()}")
        # Construct new reference sequences
        target_ref_seqs = build_ref_from_ext_coord(
            dict_max_ext, ref_dict, int_to_id=int_to_id, strand_ids=strand_ids
        )

        return target_ref_seqs

    except Exception as e:
        logging.error(f"Building reference failed: {e}")
        return False


## Main class to process files for a given barcode
class ExtensionQuery:
    def __init__(
        self,
        barcode,
        params: DorasParams,
        folder_path: Path,  ## where the output of DORAS will be stored
        raw_fastq_path: Path,  ## fastq_pass
        mlst_ref_genes_path: Path,
        target_size: int,
        parent_logger: logging.Logger,
        pause_event=asyncio.Event(),
        mapq=60,
        starting_phase=Phase.EXTENSION,
        st=None,
        ratio_map_to_ext=0.8,
        overwrite=False,
        clean_up=False,
        quantile_length=0.90,
        min_depth_consensus=20,
        polishing_interval=60,  # Time in seconds
    ):
        if st is None:
            st = {"ST": "Unknown", "final_timepoint": None}
        self.extension_length_set = False
        self.barcode = barcode
        self.phase = starting_phase
        self.completed = False  # The isolate ST is solved
        self.logging = parent_logger.getChild(f"processor_{self.barcode}")
        self.fastqs_path = Path(raw_fastq_path) / barcode
        self.output_dir: Path = (
            folder_path / barcode
        )  # the files will be found the folder containing the barcode
        self.postprocessing_dir = (
            self.output_dir / f"{barcode}_extended_ref_postprocessing"
        )
        self.concatenated_file = (
            self.output_dir / f"{self.barcode}_concatenated.fastq.gz"
        )  # file containing all the
        self.base_ref = FastaHolder(
            mlst_ref_genes_path
        )  # TODO this could take a different ref for each barcode if each barcode is a different species
        self.query_ref = None  # The reference used to query the st database'
        self.number_of_references = sum(1 for _ in self.base_ref.get_dict_fasta())
        self.current_extended_sizes = {
            header: {"up": 0, "down": 0}
            for header in self.base_ref.get_dict_fasta()  # Used to measure the extension of the two flanking regions
        }
        self.extension_positions = None  # Dict to store the coordinates of the extensions for each gene {adk:{"up":(start,end), "down":(start,end)}}
        self.polishing_interval = polishing_interval
        self.previously_processed_files = (
            0  # To keep track of the number of files processed during the last polish
        )
        self.last_medaka_polish_time = time.time()
        self.params = params
        # self.mapq = mapq/
        self.mapq_filter = mapq
        self.clean_up = clean_up
        self.ratio_map_to_ext = ratio_map_to_ext
        self.extended_ref = None  # Initial ref is the base ref of the mlst gene associated with the barcode
        self.overwrite_extended_ref = overwrite
        self.overwrite = overwrite
        self.quantile_length = quantile_length  # The quantile used to determine the longest extension to consider (up and down stream the ROI)
        # self.min_quantile = self.params.min_quantile # The quantile used to determine the longest extension to consider (up and down stream the ROI)
        self.extension_length = target_size  # Base extension long enough to prevent the process to stop, but should be updated during the run
        self.q50_extension_length = target_size
        self.st = st
        self.min_consensus_depth = (
            min_depth_consensus  # Minimum depth threshold for consensus verification
        )
        self.dict_max_ext = None  # Dict to store the coordinates of the max extension for each gene {adk:{"up":[query_id, up_ext_coord1, up_ext_coord2], "down":[query_id, down_ext_coord1, down_ext_coord2]}}
        self.store_extensions = set()
        self.target_size = {
            header: self.extension_length - length
            for header, length in self.base_ref.get_sequence_lengths().items()
        }
        self.roi_positions = {
            header: None for header in self.target_size
        }  # Dict to store the positions of the regions of interest in the extended ref {adk:(start,end)} used to obtain the depth of the region
        self.processed_files = set()  # TODO to remove
        self.pause_event = pause_event
        self.query_event = (
            asyncio.Event()
        )  # Event to signal when the query phase can start
        self.lock = threading.Lock()
        self.complete_headers: Dict[str, Union[bool, str]] = {
            header: (False, time.ctime()) for header in self.target_size
        }
        self.list_nodes = {
            header: [] for header in self.target_size
        }  # List of the coordinates of all reads used to create the new ref (Used to be able to track back to identify bad reads or chimeras)
        self.extension_complete = False

    def initialize_extended_ref(self):
        # If a final sequence is already present the process can be terminated, in case a certain gene cannot reach the goal
        # we can restart with a lower quantile to finalize the process
        # if Path(self.output_dir.res)
        try:
            Path(self.output_dir).mkdir(parents=True, exist_ok=True)
        except FileExistsError:
            self.logging.error("")
            raise asyncio.CancelledError()
        if (
            Path(self.output_dir / f"{self.barcode}_final_ref.fasta").exists()
            and not self.overwrite
        ):
            self.logging.warning(
                f"A final sequence was found, terminating the process of extension for {self.barcode}"
            )
            self.extension_complete = True
            return
        elif (
            not Path(self.output_dir / f"{self.barcode}_final_ref.fasta").exists()
            and self.phase == Phase.QUERY_ST
        ):
            self.logging.error("No final ref available, cancelling run")
            raise asyncio.CancelledError()

        if Path(self.concatenated_file).exists():
            if self.clean_up:
                self.logging.warning(
                    f"Clean up is set to :{self.clean_up}, Erasing previous concatenated file : {self.concatenated_file}"
                )
                self.concatenated_file.unlink()

        output_path = self.output_dir / f"{self.barcode}_extended_ref.fasta"
        self.logging.info(
            f"Initializing the extended reference, verifying if a previous sequence was started at {output_path} if overwrite is ({self.overwrite}) "
        )
        if output_path.exists() and not self.overwrite:
            self.extended_ref = FastaHolder(str(output_path))
            self.logging.info(
                f"Previous extended sequence found, continuing the work: {output_path}"
            )
        elif self.overwrite_extended_ref and output_path.exists():
            # Handle the case where the file does not exist, e.g., set to None or raise an exception
            self.logging.info(f"Erasing previous extended: {output_path}")
            output_path.unlink()
            self.extended_ref = (
                self.base_ref
            )  # or raise FileNotFoundError(f"The file {output_path} does not exist.")
        else:
            self.logging.info("Initializing the extended reference ")
            self.extended_ref = (
                self.base_ref
            )  # or raise FileNotFoundError(f"The file {output_path} does not exist.")

        if self.params.min_consensus_depth:
            self.min_consensus_depth = self.params.min_consensus_depth

    async def concatenate_fastq_files_async(self):
        try:
            concatenated_file = (
                self.output_dir / f"{self.barcode}_concatenated.fastq.gz"
            )
            if not self.concatenated_file:
                self.concatenated_file = concatenated_file
            fastq_files = [
                f
                for f in Path(self.fastqs_path).rglob("*.fastq.gz")
                if f != concatenated_file
            ]
            if not fastq_files:
                self.logging.critical(
                    f"No FASTQ files found in {self.output_dir}, .fastq.gz only will be detected"
                )
                return

            # Prevent to concatenate when files were already concatenated
            fastq_files_set = set(fastq_files)
            self.logging.debug(
                f"current files {fastq_files_set} and processed files: {self.processed_files}"
            )
            if fastq_files_set.issubset(self.processed_files):
                self.logging.warning("All files are concatenated")
                return

            # Determine file mode: 'wb' to overwrite, 'ab' to append
            file_mode = (
                "wb" if self.overwrite or not concatenated_file.exists() else "ab"
            )

            if self.overwrite and concatenated_file.exists():
                self.logging.info(f"Overwriting existing file: {concatenated_file}")
            async with aiofiles.open(concatenated_file, file_mode) as outfile:
                for fastq in fastq_files:
                    if str(fastq) not in self.processed_files:
                        async with aiofiles.open(fastq, "rb") as infile:
                            content = await infile.read()
                            await outfile.write(content)
                        self.processed_files.add(str(fastq))

                self.logging.info(
                    f"Concatenated {len(fastq_files)} FASTQ files into {concatenated_file}"
                )
        except asyncio.CancelledError:
            self.logging.warning("Concatenation cancelled.")
            raise
        except FileNotFoundError:
            self.logging.warning("No FASTQ files found for {self.barcode}")
        except Exception as e:
            self.logging.error(f"Error concatenating FASTQ files: {e}")
            raise

    async def extension_process(self, file_path: str):
        """
        asynchronous function to control the extension of the reference from the base reference
        """
        ##########################################################################
        ##########################################################################
        # Define the extension length if not defined yet #TODO create a method for that
        try:
            if self.extension_complete:
                self.logging.info(
                    f"Processing for {self.barcode} is complete. Skipping."
                )
                return

            # Step 1: Alignment
            try:
                paf = await mpaf.perform_alignment_async(
                    self.concatenated_file,
                    self.extended_ref.get_filepath(),  # TODO Modified to take only the concatenated file
                )
                if paf is None or paf.empty:
                    self.logging.error(
                        f"No alignments found for {self.concatenated_file} against {self.extended_ref.get_filepath()}"
                    )
                    return False

            except Exception as r:
                self.logging.error(
                    f"Alignment failed for {self.concatenated_file}: {r}"
                )
                return False
            if not self.extension_length_set:
                try:
                    await self.infer_optimal_ref_length(paf)
                    self.extension_length_set = True
                except Exception as r:
                    self.logging.error(f"Failed to initialize extension length: {r}")
            ##########################################################################
            ##########################################################################
            # Extension process itself

            ## Exclude the targets (asmlst genes) that are already completed (fully extended)
            complete_headers = [
                header
                for header, value in self.complete_headers.items()
                if value[0] == True
            ]
            paf = paf[~paf["target_name"].isin(complete_headers)]
            # Step 3: Process PAF and Compute Extensions

            (
                dict_max_ext,
                int_to_id,
                ids_strand,
            ) = await process_paf_and_compute_extensions(
                paf,
                mapq_filter=self.mapq_filter,
                up_trim=500,
                down_trim=500,
                min_aln_ext_ratio=self.ratio_map_to_ext,  # TODO add these parameters to a higher level High Priority
            )
            ### Put the dict_max_extension in a variable of the class
            self.dict_max_ext = dict_max_ext
            if self.dict_max_ext is None:
                self.logging.error("Failed to compute max extensions")
                return False
            ## fastq and ref fasta into one single dict
            fastq_dict = self.extended_ref.get_dict_fasta()

            fastq_dict.update(FastqHolder(filepath=file_path).get_dict_fasta())

            # Step 3: Build Reference

            result = construct_reference_from_extensions(
                self.dict_max_ext,
                ref_dict=fastq_dict,
                int_to_id=int_to_id,
                strand_ids=ids_strand,
            )
            if not result:
                return False
            ##
            ## Add the current extension to the dictionary if (usually the case) not all are present in the result (the new extended refs)
            number_refs_exten = sum(1 for _ in result)
            if number_refs_exten < self.number_of_references:
                base_ref_dict = self.extended_ref.get_dict_fasta()
                for key, value in base_ref_dict.items():
                    if key not in result:
                        result[key] = value  # TODO change to default dict

            # Write the new ref to disk
            output_path = str(self.output_dir / f"{self.barcode}_extended_ref.fasta")
            save_dict_to_fasta(result, output_path)
            self.extended_ref = FastaHolder(filepath=output_path)
            # Write a method with this part

            ## Update the dict_max_ext to have the header name instead of the int id for easier manipulation in the next steps
            self.dict_max_ext = {
                int_to_id[header]: values
                for header, values in self.dict_max_ext.items()
            }
            # dict_extensions = get_extensions_by_mapping(query_base_ref=self.base_ref.get_filepath(),extended_ref=self.extended_ref.get_filepath())
            ## Tracking the sizes of the extensions #TODO Change this part to integrate mapping of the original base_ref
            for header in self.current_extended_sizes:
                # cur = self.current_extended_sizes[header]
                if header in dict_max_ext:
                    new = dict_max_ext.get(header)
                    # track the nodes used and their coordinates to verify at the end
                    # TODO Add a storage of the reads used to extend to eliminate the ones already used High priority
                    self.list_nodes[header].append(
                        {
                            direction: (int_to_id[coord[0]],) + coord[1:]
                            for direction, coord in new.items()
                        }
                    )
            all_done = True
            #################################################################################
            #################################################################################
            # Verify the extension status
            self.extension_positions = get_extensions_by_mapping(
                self.base_ref.get_filepath(), self.extended_ref.get_filepath()
            )
            # for head, ext in self.current_extended_sizes.items():
            for header, ext in self.extension_positions.items():
                cur_up_ext_length = ext["up"][2]
                self.roi_positions[header] = (ext["up"][1], ext["down"][0])
                cur_down_ext_length = ext["down"][2]
                self.current_extended_sizes[header]["up"] = cur_up_ext_length
                self.current_extended_sizes[header]["down"] = cur_down_ext_length

                if (
                    cur_up_ext_length < self.extension_length
                    or cur_down_ext_length < self.extension_length
                ):
                    self.logging.info(
                        f"{self.barcode}: The upstream region {header} has reached {round(cur_up_ext_length / self.extension_length, 1) * 100}% and downstream {round(cur_down_ext_length / self.extension_length, 1) * 100}% of the target length ({self.extension_length})"
                    )
                    all_done = False

                    if self.params.force_finalize:
                        self.logging
                        all_done = True

                else:
                    self.complete_headers[header] = (
                        True,
                        time.ctime(),
                    )  # Keep track of the complete headers
                    logging.info(
                        f"Sequence {header} for barcode {self.barcode} has reached the required extension of {self.target_size.get(header)}"
                    )
            # Update base_ref to the new extended reference after the first sequence was successfully added
            self.extended_ref = FastaHolder(
                self.output_dir / f"{self.barcode}_extended_ref.fasta"
            )
            save_list_nodes_to_file(
                self.list_nodes, self.output_dir / f"{self.barcode}_nodes.json"
            )
            # self.processed_files.add(file_path)
            self.logging.info(
                f"Processed {file_path} for barcode {self.barcode}: Successfully updated extended reference"
            )

            ##########################################################################################
            ##########################################################################################

            # If all sequences have been sufficiently extended, mark as self.extension_complete = True
            if all_done:
                try:
                    await self.finalize_extension_process()
                except asyncio.CancelledError as r:
                    self.logging.error(
                        f"Finalization of the extension process failed {r}"
                    )

        except asyncio.CancelledError:
            self.logging.error(
                msg=f"Processing for barcode {self.barcode} was canceled."
            )

        except Exception as e:
            self.logging.error(msg=f"Error processing file: {e}")

    async def infer_optimal_ref_length(self, paf):
        """
        method to infer the optimal extension length based on the median of the lengths of the longest extensions

        """
        if self.extension_length == 0:
            self.logging.info(
                f"{self.barcode}: Inferring the optimal size from the reads"
            )
            if paf["query_length"].shape[0] < 500:
                self.logging.warning("The first batch of read was below 500")
            full_length_reads = paf["query_length"]
            # remove all reads ablow 50kb
            reads_sample = full_length_reads[full_length_reads <= 50000]
            try:
                _r = find_optimal_quantile(
                    reads_sample.to_numpy(), 400, genome_size=self.params.genome_size
                )
            except Exception as e:
                self.logging.error(f"Failed to infer optimal extension length: {e}")
                raise asyncio.CancelledError()

            self.extension_length = _r["ext_length"]
            # self.q50_extension_length = _r["q50"]#TODO Added
            self.quantile_length = _r["q_max"]

            self.logging.info(f"New extension length is {self.extension_length}")
            self.logging.info(f"Quantile chosen was {self.quantile_length}")
            self.target_size = {
                header: self.extension_length + length + self.extension_length
                for header, length in self.base_ref.get_sequence_lengths().items()
            }
        else:
            self.logging.info(
                f"{self.barcode}:Size of the extension chosen by user ({self.extension_length})"
            )
            self.target_size = {
                header: self.extension_length + length + self.extension_length
                for header, length in self.base_ref.get_sequence_lengths().items()
            }
        try:
            with (
                open(
                    f"{self.output_dir}/{self.barcode}_target_extension_size.bed",
                    "w",
                ) as bedfile
            ):  # Write the extensions of each target to a bed file heade, q99 length, final size of the
                bedfile.write(
                    f"header_barcode\t{self.quantile_length}_extension\tq50_extension\ttotal_size_after_extension\n"
                )

            with (
                open(
                    f"{self.output_dir}/{self.barcode}_target_extension_size.bed",
                    "a",
                ) as bedfile
            ):  # Write the extensions of each target to a bed file heade, q99 length, final size of the
                for header, size in self.target_size.items():
                    self.logging.info(
                        f"Target {header} for barcode {self.barcode} has an extension of {size} based on the quantile {self.quantile_length} of the read length distribution"
                    )
                    bedfile.write(
                        f"{header}_{self.barcode}\t{str(self.extension_length)}\t{self.q50_extension_length}\t{str(size)}\n"
                    )
        except Exception as r:
            self.logging.error(f"error write bedfile {r}")

    async def finalize_extension_process(self):
        """
        Function to finalize the extension process by trimming the extended reference and polishing it before moving to the query phase
        """
        self.extension_complete = True
        # Once all sequences have been extended they are trimmed to the right size upstream and downstream
        output_path = str(self.output_dir / f"{self.barcode}_final_ref.fasta")
        self.logging.info(
            f"Finalizing process, trimming the extended reference of {self.barcode} and saving it to {output_path}"
        )

        try:
            trimmed_final_extended_reference = trim_sequences(
                extension_positions=self.extension_positions,
                fasta_dict=self.extended_ref.get_dict_fasta(),
                target_size=self.extension_length,
                q50_extension=self.q50_extension_length,
                barcode=self.barcode,
            )
            save_dict_to_fasta(
                data_dict=trimmed_final_extended_reference,
                output_file=output_path,
            )

        except Exception as r:
            logging.error(f"Trimming for {self.barcode} failed {r}")
        save_list_nodes_to_file(
            self.list_nodes, self.output_dir / f"{self.barcode}_nodes.json"
        )
        # Polishing of the new ref before moving to the QUERY_PHASE
        try:
            polish_with_medaka(
                input_file=self.output_dir / f"{self.barcode}_concatenated.fastq.gz",
                reference_fasta=str(
                    self.output_dir / f"{self.barcode}_final_ref.fasta"
                ),
                output_dir=self.postprocessing_dir,
                bacteria=True,
            )
        except ValueError as b:
            self.logging.error(b)
            self.logging.warning(
                "Medaka could not be polished with `--bacteria`flag falling back to standard mode."
            )

            polish_with_medaka(
                input_file=self.output_dir / f"{self.barcode}_concatenated.fastq.gz",
                reference_fasta=str(
                    self.output_dir / f"{self.barcode}_final_ref.fasta"
                ),
                output_dir=self.postprocessing_dir,
                bacteria=False,
            )
        except RuntimeError as r:
            logging.error(
                f"Polishing of the final extended seq {self.barcode} failed {r}"
            )
            raise r
        consensus_seq = self.postprocessing_dir / "consensus.fasta"
        # Copy content from consensus.fasta to final_ref.fasta #TODO use shutil
        with open(consensus_seq, "r") as src:
            with open(output_path, "w") as dest:
                dest.write(src.read())
        # return asyncio.exceptions.CancelledError
        self.logging.info(
            f"\nAll sequences for {self.barcode} have reached target. Processing complete.\n"
        )

        self.logging.warning(
            f"Extension for {self.barcode} is finished, user will switch to phase {Phase.QUERY_ST} as soon as all other samples are finished"
        )
        await asyncio.sleep(5)

    ## Function to handle the switch to the QUERY_PHASE that will take place in a new folder that will
    async def main(self):
        """
        Main asynchronous function to control the workflow.
        """
        try:
            self.completed = False
            # self.initialize_extended_ref()
            while not self.completed:
                if self.phase == Phase.EXTENSION and not self.extension_complete:
                    try:
                        await self.concatenate_fastq_files_async()

                    except asyncio.CancelledError:
                        self.logging.error(
                            msg=f"Concatenation process for barcode {self.barcode} was canceled."
                        )
                        raise
                    except Exception as e:
                        self.logging.error(f"Error Concatenating Files together {e}")

                    self.logging.debug(f"Currently in phase {self.phase}")
                    self.logging.debug(f"Performing extension of {self.barcode}")
                    try:
                        await self.extension_process(
                            self.output_dir / f"{self.barcode}_concatenated.fastq.gz"
                        )
                    except asyncio.CancelledError:
                        self.logging.error(
                            msg=f"Extension process for barcode {self.barcode} was canceled."
                        )
                        raise
                    except Exception as r:
                        self.logging.error(
                            f"Problem during the extension process of {self.barcode} : {r}"
                        )

                elif self.phase == Phase.QUERY_ST and not self.completed:
                    self.logging.info(f"Starting phase {self.phase}")

                    try:
                        if Path(self.concatenated_file).exists() and self.clean_up:
                            self.logging.warning(
                                f"Clean up is set to :{self.clean_up}, Erasing previous concatenated file : {self.concatenated_file}"
                            )
                            os.remove(self.concatenated_file)
                            self.clean_up = False  # To prevent to clean up again in the next loop and erase new files

                        await asyncio.sleep(5)

                        await self.concatenate_fastq_files_async()

                        # Trimming the sequences mainly
                        self.prepare_for_query_phase()

                        # Main function of the query process the trimmed final ref will be polished and queried to the database
                        await self.polish_query_ref()
                        # await self.schedule_polishing()
                    except asyncio.CancelledError:
                        self.logging.error(
                            msg=f"Query phase for barcode {self.barcode} was canceled."
                        )
                        raise
                    except Exception as r:
                        self.logging.error(
                            f"Problem during the query process of {self.barcode} : {r}"
                        )
                else:
                    ## In case you run only the extension phase it breaks if the extension is finished but not the whole run otherwise it loops forever #TODO find better fix for that
                    break

        except asyncio.CancelledError as r:
            logging.info(
                msg=f"Processing for barcode {self.barcode} has been canceled \n Due to {r}."
            )

        except Exception as r:
            logging.error(msg=f"Error processing files: {r}")

    def prepare_for_query_phase(self):
        """
        Wrapper that prepares the second phase of the process,
        the final ref will be trimmed to the right size before querying the BIGSdb database
        and more
        :param self
        """
        self.get_ext_and_trim_ref(trim_length=1000)
        polish_with_medaka(
            input_file=self.output_dir / f"{self.barcode}_concatenated.fastq.gz",
            reference_fasta=str(self.query_ref),
            output_dir=self.postprocessing_dir,
            bacteria=False,
        )

        self.logging.info(
            f"Finalizing process, trimming the extended reference of {self.barcode} and saving it to {str(self.query_ref)}"
        )

    def get_ext_and_trim_ref(self, trim_length: int):
        """
        Method to obtain the trimmed final reference after extension to query to the
        BIGSdb database, longer sequences are not accepted, therefore trimming is required
        Wrapper method to reduce the size of the final reference with respect to the ROI (default self.extension_size)

        :param trim_length: Length to trim the final reference with respect to the ROI
        :type trim_length: int
        """
        self.query_ref = self.output_dir / f"{self.barcode}_final_ref_trimmed.fasta"
        path_final_fasta = str(self.output_dir / f"{self.barcode}_final_ref.fasta")
        final_ref = FastaHolder(path_final_fasta)
        new_extensions = get_extensions_by_mapping(
            self.base_ref.get_filepath(), path_final_fasta
        )

        if trim_length == 0:
            trim_length = self.extension_length
            self.logging.debug(f"Trim length set to extension length {trim_length}")

        trimmed_final_extended_reference = trim_sequences(
            extension_positions=new_extensions,
            fasta_dict=final_ref.get_dict_fasta(),
            target_size=trim_length,
            q50_extension=self.q50_extension_length,
            barcode=self.barcode,
        )
        save_dict_to_fasta(
            data_dict=trimmed_final_extended_reference,
            output_file=self.query_ref,
        )

    async def schedule_polishing(self):

        # previously_processed_files = len(self.processed_files)
        # Get current processed files count
        current_processed_files = len(self.processed_files)

        # Check if new files have been processed since last polish
        new_files_processed = current_processed_files > self.previously_processed_files

        # current_time = time.time()
        # time_since_last_polish = current_time - self.last_medaka_polish_time
        # time_until_next_polish = max(
        #     0, self.polishing_interval - time_since_last_polish
        # )

        # self.logging.info(
        #     f"Time since last polish: {int(time_since_last_polish)}s, "
        #     f"Next polishing in {int(time_until_next_polish)} seconds"
        # )
        if new_files_processed:
            self.previously_processed_files = current_processed_files
            # if time_since_last_polish >= self.polishing_interval:
            try:
                polish_with_medaka(
                    input_file=self.output_dir
                    / f"{self.barcode}_concatenated.fastq.gz",
                    reference_fasta=str(self.query_ref),
                    output_dir=self.postprocessing_dir,
                    bacteria=True,
                )
                # self.last_medaka_polish_time = current_time
                self.logging.info("Medaka polishing completed successfully")
            except ValueError as b:
                self.logging.error(f"Value error during medaka polishing: {b}")
                self.logging.warning(
                    "Medaka could not be polished with `--bacteria`flag falling back to standard mode."
                )
                try:
                    polish_with_medaka(
                        input_file=self.output_dir
                        / f"{self.barcode}_concatenated.fastq.gz",
                        reference_fasta=str(self.query_ref),
                        output_dir=self.postprocessing_dir,
                        bacteria=False,
                    )
                    # self.last_medaka_polish_time = current_time
                    self.logging.info(
                        "Medaka polishing completed successfully in fallback mode"
                    )
                except Exception as e:
                    self.logging.error(
                        f"Error during medaka polishing in fallback mode: {e}"
                    )
                    # Wait before retrying
                    await asyncio.sleep(60)

            except Exception as e:
                self.logging.error(f"Error during medaka polishing: {e}")
                # Wait before retrying
                await asyncio.sleep(60)

        # Sleep until next check or next polishing time, whichever comes first
        # sleep_duration = (
        #     min(60, time_until_next_polish) if time_until_next_polish > 0 else 1
        # )
        await asyncio.sleep(5)

    async def run_polish_with_medaka(self):

        sequence_to_query = self.postprocessing_dir / "consensus.fasta"
        if not sequence_to_query.exists():
            sequence_to_query = self.query_ref
            self.logging.info(f"Polishing final ref {self.barcode}")
        self.logging.debug(
            f"Running medaka polishing for {self.barcode} with {sequence_to_query} "
        )
        try:
            polish_with_medaka(
                input_file=self.output_dir / f"{self.barcode}_concatenated.fastq.gz",
                reference_fasta=str(self.query_ref),
                output_dir=self.postprocessing_dir,
                bacteria=True,
            )
        except RuntimeError as e:
            # This handles other medaka-related errors (executable not found, etc.)
            print(f"Runtime error in Medaka: {e}")
            logging.error(f"Medaka runtime error: {e}")
            return False
        except Exception as e:
            logging.error(f"Error during medaka polishing: {e}")

    async def polish_query_ref(self):
        """
        The entire query phase with polishing logics takes place here

        """

        while not self.completed:
            try:
                await self.concatenate_fastq_files_async()
                await self.schedule_polishing()
                if not self.query_event.is_set():
                    if await self.verify_depth_coverage():
                        self.logging.info(
                            f"Depth coverage verification passed for {self.barcode}"
                        )
                        self.query_event.set()
                        await asyncio.sleep(5)
                else:
                    #  await self.query_db_print_result()
                    await asyncio.sleep(5)

            except asyncio.CancelledError:
                self.logging.error(
                    msg=f"Polishing and querying process for barcode {self.barcode} was canceled."
                )
                raise
            except Exception as r:
                self.logging.error(f"Problem with polishing and querying : {r}")
                return

    # Verify the depth of the consensus using pysam (bam file)
    async def verify_depth_coverage(self):
        try:
            bam_file = self.postprocessing_dir / "calls_to_draft.bam"
            if not bam_file.exists():
                self.logging.error(f"BAM file {bam_file} does not exist.")
                return False

            samfile = pysam.AlignmentFile(str(bam_file), "rb")

            consensus_file = self.postprocessing_dir / "consensus.fasta"
            if not consensus_file.exists():
                self.logging.error(f"Consensus file {consensus_file} does not exist.")
                return False
            consensus_sequences = FastaHolder(consensus_file).get_dict_fasta()

            out = True
            for record_id, sequence in consensus_sequences.items():
                seq_length = len(sequence)
                total_coverage = 0
                for pileupcolumn in samfile.pileup(record_id, 0, seq_length):
                    total_coverage += pileupcolumn.nsegments
                    if pileupcolumn.nsegments < self.params.min_consensus_depth:
                        self.logging.debug(
                            f"Low coverage (<{self.params.min_consensus_depth})  for {record_id}: {total_coverage} reads"
                        )
                        # break
                average_depth = total_coverage / seq_length if seq_length > 0 else 0
                if average_depth < self.params.min_consensus_depth:
                    out = False
                    self.logging.warning(f"Depth insufficient for {self.barcode}")
                self.logging.info(f"Mean coverage for {record_id}: {average_depth:.2f}")
            if not out:
                self.logging.info(
                    f"Depth  insufficient  for {self.barcode} (min={self.params.min_consensus_depth} for each allele)"
                )
            samfile.close()
            return out
        except asyncio.CancelledError:
            self.logging.error(
                msg=f"Depth coverage verification for barcode {self.barcode} was canceled."
            )
            return False
        except FileNotFoundError as fnf_error:
            self.logging.error(f"File not found: {fnf_error}")
            return False
        except Exception as e:
            self.logging.error(f"Error verifying depth coverage: {e}")
            return False

    async def query_db_print_result(self):
        new_var = self.postprocessing_dir / "consensus.fasta"
        encoded_consensus = await read_encode_fasta_file(new_var)
        try:
            self.logging.info("Printing query results to file")
            # st_dict = await fetch_consensus_data_new(self.params)
            st_dict = await query_st_database(
                encoded_sequence=encoded_consensus,
                scheme_url=self.params.url,
            )
            self.logging.debug(f"Current results are {st_dict}")

            st_df = pd.DataFrame.from_dict(st_dict, orient="index")
            st_df.to_csv(path_or_buf=self.output_dir / f"st_results_{self.barcode}.csv")
            self.logging.info(f"ST results for {self.barcode} written to CSV")
            if st_dict.get("fields"):
                st = st_dict.get("fields")["ST"]
                st_df.to_csv(
                    path_or_buf=self.output_dir / f"st_results_{self.barcode}.csv"
                )
                self.logging.info(f"{self.barcode}: ST solved {st}")
                self.completed = True
                self.st = st
            elif not st_dict.get("fields") or st_dict.get("fields")["ST"] == "None":
                self.logging.info(f"{self.barcode}: ST not yet solved, continuing")
        except Exception as r:
            self.logging.error(f"Problem with the query of the consensus sequence {r}")
            return


class TerminateTaskGroup(Exception):
    """Exception raised to terminate a task group."""


class ExtensionFinished(Exception):
    """Extension Phase is finished successfully."""


async def force_terminate_task_group():
    """Used to force termination of a task group."""
    raise TerminateTaskGroup()
