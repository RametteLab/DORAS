#!/usr/bin/env python

import argparse
import asyncio
import logging
from pathlib import Path
from src.doras_extension import Phase
from src.doras_manager import DorasManager
from src.config import load_params_from_toml


# from datasets.strandness.RandSeqGen import RandSeq

# Configure logging for asmlst_asm.p
#
# New function for BLAST verification using the provided command
parser = argparse.ArgumentParser()
parser.add_argument("--toml", type=str, help="Valid TOML file", required=True)
# parser.add_argument("--test_mode",store, help="Valid TOML file",required=True)
parser.add_argument(
    "--phase",
    type=Phase,
    choices=list(Phase),
    default=Phase.EXTENSION,
    help="Start the query phase",
)
parser.add_argument(
    "--size",
    type=int,
    default=0,
    help="Desired size of the extension, default is 0 and will be dynamically computed. Meant to be used for debugging",
)
parser.add_argument(
    "--overwrite_previous_ref",
    action="store_true",
    help="Whether to overwrite the previous reference when restarting the script",
)
parser.add_argument(
    "--force_finalize",
    action="store_true",
    help="Finish the trimming process without completing the extension",
)
parser.add_argument(
    "--no_clean_up", action="store_true", help="Do not remove the concatenanted fastq file"
)
parser.add_argument("--debug", action="store_true", help="Enable debug logging")

args = parser.parse_args()

params = load_params_from_toml(args.toml)
logger = logging.getLogger("DORAS")
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(name)s - %(message)s')
# Log to file
log_file = f'{params.experiment_name}.log'
file_handler = logging.FileHandler(log_file)
file_handler.setFormatter(formatter)
file_handler.setLevel(logging.DEBUG)
## Stream printed to console
stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.INFO)
stream_handler.setFormatter(formatter)
if args.debug:
    stream_handler.setLevel(logging.DEBUG)

logger.setLevel(logging.DEBUG)
logger.addHandler(file_handler)
logger.addHandler(stream_handler)
logger.propagate = False

async def main():
    mlst_genes_path = params.mlst_genes_path
    test_folder = params.output_dir_path
    if params.test_mode:
        test_folder = Path(params.output_dir_path)

    barcodes = params.sample_names
    params.force_finalize = args.force_finalize

    manager = DorasManager(
        mlst_ref_genes_path=mlst_genes_path,
        barcodes=barcodes,
        folder_path=test_folder,
        mapq=params.min_map_quality,
        overwrite=args.overwrite_previous_ref,
        target_extension_size=args.size,
        params=params,
        clean_up=args.no_clean_up,
        test_mode=params.test_mode,  # TODO add params for TEST_MODE
        phase=args.phase,
        logger=logger,
    )

    processor_task = asyncio.create_task(manager.main())

    await processor_task


if __name__ == "__main__":
    asyncio.run(main())
