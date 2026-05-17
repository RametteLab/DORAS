import asyncio
import datetime
import pandas as pd
from pathlib import Path
import aiohttp
import logging
import requests
import base64
import json
import marimo as mo
from Bio import SeqIO

from .config import DorasParams

# Set up logging with debug mode included
logger = logging.getLogger("QueryDB")
logging.basicConfig(
    level=logging.INFO,  # Changed to DEBUG for detailed logs
    format="%(asctime)s - %(levelname)s - %(message)s",
)

BASE_API = {
    "PubMLST": "https://rest.pubmlst.org",
    "Pasteur": "https://bigsdb.pasteur.fr/api",
}


def generate_dataframe(response_list):
    data = []
    for response in response_list:
        exact_matches = response.get("exact_matches", [])
        best_match = response.get("best_match", {})
        if exact_matches:
            for match in exact_matches:
                data.append(
                    {
                        "Isolate": response.get("Isolate", "None"),
                        "Loci": response.get("Loci", "None"),
                        "complete": True,
                        "allele_id": match.get("allele_id"),
                        "href": match.get("href"),
                    }
                )
        elif best_match:
            data.append(
                {
                    "Isolate": response.get("Isolate", "None"),
                    "Loci": response.get("Loci", "None"),
                    "complete": False,
                    "allele_id": best_match.get("allele_id"),
                    "identity": best_match.get("identity"),
                    "alignment": best_match.get("alignment"),
                }
            )
    df = pd.DataFrame(data)
    return df


def get_loci_status(params: DorasParams):
    """
    params are defined by the DorasParams Pydantic class which is defined by the toml file doras.toml
    """
    response_list = []
    for bc in mo.status.progress_bar(params.sample_names):
        logger.info(f"Fetching alignment details from {params.scheme}")
        logger.debug(bc)  # Changed to debug
        consensus = (
            Path(params.path)
            / str(bc)
            / f"{bc}_extended_ref_postprocessing"
            / "consensus.fasta"
        )
        logger.debug(f"Fetching data for {consensus}")
        fa = SeqIO.parse(consensus, "fasta")

        for seq in fa:
            _sequence_data = {"sequence": f"{seq.seq}"}
            loci = seq.id.split("_")[0]
            _api_url = f"{params.selected_db}/loci/{loci}/sequence"
            logger.debug(f"Fetching alignement {_api_url} ")
            try:
                _response = requests.post(
                    _api_url,
                    data=json.dumps(_sequence_data),
                    params=json.dumps(params.request_params),
                )
                _text = _response.json()
                logger.debug(_text)  # Changed to debug
                _text["Isolate"] = bc
                _text["Loci"] = loci
                response_list.append(_text)
            except Exception:
                logging.error(
                    f"Alignment details failed for {bc}"
                )  # Fixed variable name

    return response_list


def get_loci_status_single(params: DorasParams):
    """
    params are defined by the FetchConsensusDataParams Pydantic class which is defined by the toml file doras.toml
    """
    fa = SeqIO.parse(consensus, "fasta")
    params = {"details": "true", "base64": "true"}
    list_loci = []
    list_isolates = []
    for seq in fa:
        list_isolates.append(sample_name)
        _sequence_data = {"sequence": f"{seq.seq}"}
        loci = seq.id.split("_")[0]
        list_loci.append(loci)
        _api_url = f"{params.selected_db}/loci/{loci}/sequence"

        _response = requests.post(
            _api_url,
            data=json.dumps(_sequence_data),
            params=params,  # Corrected to use params
        )
        _text = _response.json()
        logger.debug(_text)  # Changed to debug
        logger.info(f"Processed loci: {_loci}")

    return list_loci, list_isolates


async def read_encode_fasta_file(fasta_path: Path):
    """
    Reads a FASTA file and encodes its content in base64.

    :param fasta_path: Path to the FASTA file.
    :return: Dictionary containing the base64-encoded sequence.
    """
    if not fasta_path.exists():
        logger.error(f"File {fasta_path} does not exist.")
        raise FileNotFoundError(f"File {fasta_path} does not exist.")

    with open(fasta_path, "r") as fasta_file:
        fasta_content = fasta_file.read()

    encoded_sequence = {
        "base64": "true",
        "sequence": base64.b64encode(fasta_content.encode()).decode(),
    }
    return encoded_sequence


async def query_st_database(encoded_sequence: dict, scheme_url: str) -> dict:
    """
    Queries the ST database with an encoded FASTA sequence.

    :param encoded_sequence: Dictionary containing the base64-encoded sequence.
    :param scheme_url: URL of the ST scheme database endpoint.
    :return: JSON response from the database query.
    """
    try:
        logger.info(f"Querying {scheme_url}")
        async with aiohttp.ClientSession() as session:
            async with session.post(
                scheme_url,
                data=json.dumps(encoded_sequence),
            ) as response:
                return await response.json(content_type=None)
    except Exception as e:
        logger.error(f"Request problem: {e}")
        return {"error": str(e)}


async def fetch_consensus_for_sample(params, bc):
    """Process a single sample - used for concurrent execution."""
    consensus_path = Path(
        Path(params.path) / bc / f"{bc}_extended_ref_postprocessing" / "consensus.fasta"
    )

    encoded_sequence = await read_encode_fasta_file(consensus_path)
    if not encoded_sequence:
        return {"Isolate": bc, "ST": None, "final_timepoint": None}

    try:
        matches = await query_st_database(encoded_sequence, params.url)
        logger.debug(f"ST status for {bc}")
        logger.debug(f"{matches}")

        st_value = "Unknown"
        final_time = None

        if matches.get("fields"):
            st_value = matches.get("fields", {}).get("ST", "None")
            final_time = datetime.datetime.now()
            logger.info(f"ST profile solved for barcode {bc}: ST{st_value}")

        if matches.get("exact_matches"):
            logger.debug(
                f"Number of fully identified alleles: {len(matches['exact_matches'])}"
            )

        return {"Isolate": bc, "ST": st_value, "final_timepoint": final_time}

    except Exception as e:
        logger.error(f"Request problem for sample {bc}: {e}")
        return {"Isolate": bc, "ST": None, "final_timepoint": None}


async def fetch_consensus_data_new(params, max_concurrent: int = 10):
    """
    Fetches consensus data and queries the ST database for each sample concurrently.

    :param params: Parameters defined by the FetchConsensusDataParams Pydantic class.
    :param max_concurrent: Maximum number of concurrent requests (default: 10).
    :return: Dictionary containing isolate names and their corresponding ST statuses.
    """
    # Semaphore to limit concurrent requests and prevent host timeouts
    semaphore = asyncio.Semaphore(max_concurrent)

    async def process_with_semaphore(bc):
        async with semaphore:
            return await fetch_consensus_for_sample(params, bc)

    # Create tasks for all samples to run concurrently (limited by semaphore)
    tasks = [process_with_semaphore(bc) for bc in params.sample_names]

    # Run all queries in parallel
    results = await asyncio.gather(*tasks)

    # Aggregate results
    st_status = {"Isolate": [], "ST": [], "final_timepoint": []}
    for result in results:
        st_status["Isolate"].append(result["Isolate"])
        st_status["ST"].append(result["ST"] or "Unknown")
        st_status["final_timepoint"].append(result["final_timepoint"])

    return st_status


async def query_consensus(params, details=True):
    # params of type FetchConsensusDataParams
    # details to get alignment details

    st_infos = await fetch_consensus_data_new(
        params=params
    )  # TODO FInish the params section
    df_st_infos = pd.DataFrame(st_infos)
    if details:
        list_details = get_loci_status(params)
        df_details = generate_dataframe(list_details)
        return {"details": df_details, "st": df_st_infos}
    else:
        return {"st": df_st_infos}


def verify_connection(BASEURL):
    ok = requests.get(url=BASEURL)
    if ok.status_code == 200:
        return ok.json()
    else:
        return False


def get_available_dbs(params):
    "Return a list of all available DBs (e.g E.Coli, Dipht, Clostridium...) available"
    responses = verify_connection(params.base_api)
    if responses:
        try:
            return {
                response["name"]: response["databases"][1]["href"]
                for response in responses
                if len(response["databases"]) > 1
            }
        except Exception as e:
            logging.error(f"Some field missing in reponse {e}")


def is_db_valid(BASEURL):
    raise NotImplementedError
