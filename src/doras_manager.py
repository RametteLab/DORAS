from datetime import datetime

import logging
import asyncio

from typing import List

import pandas as pd
from pathlib import Path

from .bigsdb_tools import read_encode_fasta_file, query_st_database
from .doras_extension import ExtensionQuery, TerminateTaskGroup
from .config import DorasParams, Phase
from .fastq_time_extractor import FastQTimeExtractor


class DorasManager:
    def __init__(
        self,
        barcodes: List[str],
        folder_path: str,
        mlst_ref_genes_path: str,
        mapq: int,
        test_mode: bool,
        overwrite: bool,
        clean_up: bool,
        logger: logging.Logger,
        target_extension_size: int,
        params: DorasParams,  # TODO add params from toml for the query
        phase=Phase.EXTENSION,
    ):
        self.barcodes = barcodes  # TODO use the params for that
        self.params = params
        self.total_barcodes = len(barcodes)
        self.finished = False
        self.output_dir = Path(folder_path) / params.experiment_name
        self.mlst_ref_genes_path = mlst_ref_genes_path
        self.processors = []
        self.target_extension_size = target_extension_size
        self.completed_barcodes = 0
        self.total_files_processed = 0
        self.completed_extension_phase = False
        self.mapq = mapq  # Mapping quality to filter the mapping reads
        self.overwrite = overwrite
        self.phase = phase
        self.sts = []
        self.clean_up = clean_up
        self.phase_switching_event = asyncio.Event()
        self.done_event = asyncio.Event()
        self.logging = logger
        self.logging.info(f"Starting DORAS with parameters: {self.params.model_dump_json(indent=4)}")
        self.mock_files_creators = []
        self.test_mode = test_mode
        self._setup_processors()

    def _setup_processors(self):
        self.output_dir.mkdir(parents=True, exist_ok=True)
        for barcode in self.barcodes:
            processor = ExtensionQuery(
                barcode=barcode,
                folder_path=self.output_dir,
                mlst_ref_genes_path=Path(self.mlst_ref_genes_path),
                raw_fastq_path=Path(self.params.fastq_files_path),
                target_size=self.target_extension_size,
                mapq=self.mapq,
                clean_up=self.clean_up,
                overwrite=self.overwrite,
                pause_event=self.phase_switching_event,
                parent_logger=self.logging,
                test_mode=self.test_mode,
                starting_phase=self.phase,
                params=self.params,
            )
            processor.initialize_extended_ref()
            self.processors.append(processor)
            self.logging.info(f"Initialized processor for barcode: {barcode}")

    # For test mode with pre computed fastqs
    async def _setup_mock_files(self):

        for i, barcode in enumerate(self.barcodes):
            # async def run_single_extractor():
            self.logging.info(
                f"Setting up mock files for barcode {barcode} with sample {self.params.test_samples[i]} from {self.params.test_start_time} to {self.params.test_end_time}"
            )
            extractor = FastQTimeExtractor(
                self.params.test_samples[i], end_point=self.params.test_end_time
            )
            tasks = asyncio.create_task(
                extractor.create_mock_files(
                    Path(self.params.fastq_files_path) / barcode,
                    interval_hours=self.params.test_interval,
                    stop_event=self.done_event,
                )
            )
            self.mock_files_creators.append(tasks)

    async def fetch_consensus_for_sample(self, params, bc):
        """Process a single sample - used for concurrent execution."""
        consensus_path = Path(
            Path(self.output_dir)
            / bc
            / f"{bc}_extended_ref_postprocessing"
            / "consensus.fasta"
        )

        try:
            encoded_sequence = await read_encode_fasta_file(consensus_path)

        except ValueError as e:
            self.logging.error(f"Error reading consensus data for barcode {bc}: {e}")
            raise FileNotFoundError(f"Consensus data not found for barcode {bc}")

        try:
            matches = await query_st_database(encoded_sequence, params.url)
            self.logging.debug(f"ST status for {bc}")
            self.logging.debug(f"{matches}")

            st_value = "Unknown"
            final_time = datetime.now()

            if matches.get("fields"):
                st_value = matches.get("fields", {}).get("ST", st_value)
                final_time = datetime.now()
            if st_value != "Unknown":
                self.logging.info(f"ST profile solved for barcode {bc}: ST{st_value}")
            else:
                self.logging.warning(f"ST profile not solved for barcode {bc}")
            if matches.get("exact_matches"):
                self.logging.debug(
                    f"Number of fully identified alleles: {len(matches['exact_matches'])}"
                )
            return {"Isolate": bc, "ST": st_value, "final_timepoint": final_time}
        except asyncio.CancelledError:
            self.logging.error(f"Request for sample {bc} was cancelled.")
            return {"Isolate": bc, "ST": st_value, "final_timepoint": final_time}
        except Exception as e:
            self.logging.error(f"Request problem for sample {bc}: {e}")
            return {"Isolate": bc, "ST": st_value, "final_timepoint": final_time}

    async def fetch_consensus_data_new(self, max_concurrent: int = 1):
        """
        Continuously fetches consensus data and queries the ST database for each sample as they become ready.

        :param params: Parameters defined by the FetchConsensusDataParams Pydantic class.
        :param max_concurrent: Maximum number of concurrent requests (default: 10).
        :return: Dictionary containing isolate names and their corresponding ST statuses.
        """
        semaphore = asyncio.Semaphore(max_concurrent)
        processed_samples = set()  # Track which samples have been processed

        self.logging.info(
            f"Starting continuous fetch_consensus_data_new with max_concurrent={max_concurrent}"
        )

        while True:
            # Collect samples that are ready for ST querying
            ready_samples = []
            for p in self.processors:
                if (
                    p.query_event.is_set()
                    and p.st["ST"] == "Unknown"
                    and p.barcode not in processed_samples
                ):
                    ready_samples.append(p.barcode)
                    self.logging.info(f"Sample {p.barcode} is ready for ST query.")

            # Check if we should exit (no more samples to process)
            if not ready_samples:
                all_done = all(p.st["ST"] != "Unknown" for p in self.processors)
                if all_done:
                    self.logging.info(
                        "All samples have been processed. Exiting ST query loop."
                    )
                    break
                else:
                    self.logging.debug(
                        "No samples ready for ST querying. Waiting 5 seconds..."
                    )
                    await asyncio.sleep(5)
                    continue

            # Create and run tasks for ready samples
            async def process_with_semaphore(bc):
                async with semaphore:
                    return await self.fetch_consensus_for_sample(self.params, bc)

            tasks = [process_with_semaphore(bc) for bc in ready_samples]

            ## Verify that all tasks are completed and ST not unknown
            all_completed = 0
            try:
                results = await asyncio.gather(*tasks)

                # Process results immediately after they arrive
                for bc, result in zip(ready_samples, results):
                    processed_samples.add(bc)
                    # Update the processor state
                    for p in self.processors:
                        if p.barcode == bc:
                            p.st["ST"] = result["ST"] or "Unknown"
                            p.completed = result["ST"] != "Unknown"
                            if p.completed:
                                all_completed += 1
                            p.st["final_timepoint"] = result["final_timepoint"]
                            if result["ST"] == "Unknown":
                                self.logging.warning(
                                    f"ST for sample {p.barcode} could not be determined."
                                )
                            else:
                                self.logging.info(
                                    f"ST for sample {p.barcode} is {p.st}, marking as completed."
                                )
                            break

                self.logging.info(
                    f"Processed {len(results)} samples. Remaining: {len(self.processors) - len(processed_samples)}"
                )

                # Check if all samples are done
                if all_completed == len(self.processors):
                    self.logging.info("All samples have been processed.")
                    break

            except asyncio.CancelledError:
                self.logging.error("The fetch_consensus_data_new task was cancelled.")
                raise
            except Exception as r:
                self.logging.error(f"Error during concurrent ST queries: {r}")

            # Wait a bit before next check
            await asyncio.sleep(1)

        # Return final results
        st_status = {"Isolate": [], "ST": [], "final_timepoint": []}
        for p in self.processors:
            st_status["Isolate"].append(p.barcode)
            st_status["ST"].append(p.st.get("ST", "Unknown"))
            st_status["final_timepoint"].append((p.st.get("final_timepoint", None)))

        self.logging.info(f"Final ST status: {st_status}")
        return st_status

    async def _track_progress_extension(self):
        while True:
            completed_count = sum(1 for p in self.processors if p.extension_complete)
            if self.phase == Phase.QUERY_ST:
                break
            self.completed_barcodes = completed_count
            self.logging.info(
                f"Progress of the extension phase : {completed_count}/{len(self.processors)} barcodes completed"
            )
            if completed_count == len(self.processors):
                await asyncio.sleep(2)
                self.completed_extension_phase = True
                self.done_event.set()
                self.logging.warning(
                    "All barcodes have been completed. Exiting progress tracking. The user is asked to concatenate all samples final_sequences.fasta into a single fasta file and proceed with Phase 2"
                )
                # raise TerminateTaskGroup()
                break  # Exit the loop when all barcodes are completed
            elif self.phase != Phase.EXTENSION:
                self.logging.warning("Not tracking extension anymore")

            await asyncio.sleep(5)  # Poll every second

    async def _track_progress_overall(self):
        while True:
            if self.phase == Phase.EXTENSION:
                break
            completed_count = sum(1 for p in self.processors if p.completed)
            self.sts = [p.st for p in self.processors]
            self.completed_barcodes = completed_count
            self.logging.info(
                f"Progress of the query phase (ST solved): {completed_count}/{len(self.processors)} barcodes completed"
            )
            if completed_count == len(self.processors):
                self.logging.warning(
                    "ST of all barcodes have been solved the output can be found in the st.csv for each barcode"
                )
                self.done_event.set()
                await self.genereate_report_st()
                break  # Exit the loop when all barcodes are completed
            await asyncio.sleep(10)  # Poll every 10 seconds

    async def genereate_report_st(self):
        s = await self.fetch_consensus_data_new()
        try:
            # s["timepoint"] = time.time
            st_df = pd.DataFrame.from_dict(s)

            self.logging.info(
                f"Printing the final results of the ST query {self.output_dir}/st_results.csv"
            )

            st_df.to_csv(path_or_buf=self.output_dir / "st_results.csv")

        except asyncio.CancelledError:
            self.logging.error("The ST report generation task was cancelled.")
            raise
        except Exception as r:
            self.logging.error(f"Error occured when write ST results to file : {r}")

    async def main(self):
        try:
            # Wait until at least one read file is detected
            if self.test_mode:
                self.logging.warning("Running Test mode, initiating mock files")
                await self._setup_mock_files()
                try:
                    # Start progress tasks
                    # Start processor tasks
                    tasks = [asyncio.create_task(p.main()) for p in self.processors]
                    progress_task_extension = asyncio.create_task(
                        self._track_progress_extension()
                    )
                    progress_task = asyncio.create_task(self._track_progress_overall())
                    if self.phase == Phase.QUERY_ST:
                        self.logging.info(
                            "Starting consensus data fetching and ST querying for all samples."
                        )
                        query_task = asyncio.create_task(
                            self.fetch_consensus_data_new()
                        )
                        results = await asyncio.gather(
                            progress_task,
                            *tasks,
                            *self.mock_files_creators,
                            query_task,
                            return_exceptions=True,
                        )
                    else:
                        results = await asyncio.gather(
                            progress_task_extension,
                            progress_task,
                            *tasks,
                            *self.mock_files_creators,
                            return_exceptions=True,
                        )
                    # Start mock file tasks
                    # results = await asyncio.gather(progress_task_extension,progress_task,*processor_tasks,*self.mock_files_creators,return_exceptions=True)

                    # Cancel mock file tasks
                    if self.completed_extension_phase:
                        self.logging.warning("Stopping the mockfiles creation")
                    # results = await asyncio.gather(progress_task_extension, progress_task, *processor_tasks, *self.mock_files_creators, return_exceptions=True)
                    for p in self.mock_files_creators:
                        if isinstance(p.result(), Exception):
                            self.logging.error(
                                f"An error occurred in mock file creation: {p.result()}"
                            )
                        else:
                            self.logging.info(
                                f"Mock file creation completed: {p.result()}"
                            )
                    for p in results:
                        if isinstance(p, Exception):
                            self.logging.error(f"An error occurred: {p}")
                            raise p
                        else:
                            success, data = p
                            if success:
                                self.logging.info(
                                    f"Processor completed successfully: {data}"
                                )

                            self.logging.info(f"Processor completed: {p}")

                except asyncio.CancelledError:
                    self.logging.info("The run was cancelled, all tasks stopping.")
                except TerminateTaskGroup:
                    self.logging.info(
                        "Terminating task group as all processors have completed the extension phase."
                    )
                except Exception as e:
                    self.logging.error(f"Error when running test mode: {e}")
                # progress_task = asyncio.create_task(self._track_progress_overall())
            else:
                while (
                    not list(Path(self.params.fastq_files_path).rglob("*.fastq.gz"))
                    and not self.test_mode
                ):
                    self.logging.info("Waiting for reads to be created (every 30s)...")
                    await asyncio.sleep(delay=30)  # Check every 30 seconds
                # Initialize the progress tracking task

                progress_task = asyncio.create_task(self._track_progress_overall())
                progress_task_extension = asyncio.create_task(
                    self._track_progress_extension()
                )
                # Create tasks for all processors
                tasks = [asyncio.create_task(p.main()) for p in self.processors]
                if self.phase == Phase.QUERY_ST:
                    self.logging.info(
                        "Starting consensus data fetching and ST querying for all samples."
                    )
                    query_task = asyncio.create_task(self.fetch_consensus_data_new())
                    results = await asyncio.gather(
                        progress_task_extension,
                        progress_task,
                        query_task,
                        *tasks,
                        return_exceptions=True,
                    )
                else:
                    results = await asyncio.gather(
                        progress_task_extension,
                        progress_task,
                        *tasks,
                        return_exceptions=True,
                    )

                # if progress_task_extension.done:
                #     self.logging.info("Progress tracking for extension phase completed.")
                #     await self.fetch_consensus_data_new(self.params)
                for result in results:
                    if isinstance(result, Exception):
                        self.logging.error(f"An error occurred: {result}")
                        raise result

        except asyncio.CancelledError:
            self.logging.error("The run was cancelled, all tasks stopping.")
        except Exception as e:
            self.logging.error(f"Error managing processors: {e}")
        finally:
            # Cancel the progress task and wait for it to finish
            if progress_task:
                progress_task.cancel()
                await asyncio.sleep(2)

            if progress_task_extension:
                progress_task_extension.cancel()
                await asyncio.sleep(delay=2)
