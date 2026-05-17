import tomlkit  # noqa: E402
from pydantic import BaseModel, Field
from typing import List, Dict, Optional
from enum import Enum


class Phase(Enum):
    EXTENSION = "extension"
    QUERY_ST = "query"

    def __str__(self) -> str:
        return self.value


def load_params_from_toml(file_path):
    """Load DorasParams from TOML configuration file"""
    with open(file_path, "r") as file:
        config = tomlkit.load(file)

    # Extract parameters from the TOML file
    params = {
        "experiment_name": config["experiment_name"],
        "selected_db": config["bigsdb"]["db_selected"],
        "url": config["bigsdb"]["url"],
        "base_api": config["bigsdb"]["base_api"],
        "scheme": config["bigsdb"]["scheme"],
        "fastq_files_path": config["paths"]["fastq_files_path"],
        "quantile": config["run_params"]["quantile"],
        "min_consensus_depth": config["run_params"]["min_consensus_depth"],
        "sample_names": config["sample_names"]["list"],
        "mlst_genes_path": config["mlst_genes_path"]["value"],
        "genome_size": config.get(
            "genome_size", 5000000
        ),  # Default 5MB if not specified
        "test_mode": config["test_mode"]["value"],
        "test_start_time": config["test_mode"].get("test_start_time"),
        "test_end_time": config["test_mode"].get("test_end_time"),
        "test_samples": config["test_mode"].get("test_samples"),
        "min_quantile": config["run_params"].get("min_quantile", 0.5),
        "min_map_quality": config["run_params"].get("min_map_quality", 10),
    }

    # Create an instance of DorasParams
    fetch_params = DorasParams(**params)
    return fetch_params


def save_params_to_toml(params, file_path: str):
    """Save DorasParams to a TOML file"""
    config = tomlkit.document()
    # Experiment Name
    config["experiment_name"] = params.experiment_name
    # Genome size at top level
    config["genome_size"] = params.genome_size

    # BigsDB section
    bigsdb = tomlkit.table()
    bigsdb["base_api"] = params.base_api
    bigsdb["url"] = params.url
    bigsdb["scheme"] = params.scheme
    bigsdb["db_selected"] = params.selected_db
    config["bigsdb"] = bigsdb

    # MLST genes path section
    mlst_path = tomlkit.table()
    mlst_path["value"] = params.mlst_genes_path
    config["mlst_genes_path"] = mlst_path

    # Sequencing directory path section
    paths = tomlkit.table()
    paths["fastq_files_path"] = params.fastq_files_path
    paths["output_dir"] = params.output_dir_path
    config["paths"] = paths

    # Sample names section
    sample_names = tomlkit.table()
    sample_names["list"] = tomlkit.array(params.sample_names)
    config["sample_names"] = sample_names

    # Run parameters section
    run_params = tomlkit.table()
    run_params["quantile"] = params.quantile
    run_params["min_quantile"] = params.min_quantile
    run_params["min_map_quality"] = getattr(params, "min_map_quality", 10)
    run_params["min_consensus_depth"] = params.min_consensus_depth
    config["run_params"] = run_params

    # Test mode section
    test_mode = tomlkit.table()
    test_mode["value"] = params.test_mode or False

    if hasattr(params, "fastq_files_path") and params.fastq_files_path:
        test_mode["fastq_files_path"] = params.fastq_files_path

    if hasattr(params, "test_samples") and params.test_samples:
        test_mode["test_samples"] = tomlkit.array(params.test_samples)

    test_mode["test_start_time"] = params.test_start_time or 0
    test_mode["test_end_time"] = params.test_end_time or 1

    if hasattr(params, "test_start_time_query"):
        test_mode["test_start_time_query"] = params.test_start_time_query
    if hasattr(params, "test_end_time_query"):
        test_mode["test_end_time_query"] = params.test_end_time_query

    config["test_mode"] = test_mode

    # Write to file
    with open(file_path, "w") as file:
        tomlkit.dump(config, file)
    return True


class DorasParams(BaseModel):
    genome_size: int
    experiment_name: str = Field(default="DORAS")
    selected_db: str
    base_api: str = Field(default="https://rest.pubmlst.org")
    scheme: str = Field(default="schemes/1/")
    url: str = Field(
        default="https://rest.pubmlst.org/db/pubmlst_escherichia_seqdef/schemes/1/sequence"
    )
    request_params: Dict = Field(default={"details": "true", "base64": "true"})
    output_dir_path: str = Field(
        default=".", title="Output path of the experiment"
    )  # Not implemented yet
    sample_names: List[str]
    force_finalize: Optional[bool] = Field(default=False)
    mlst_genes_path: str
    fastq_files_path: str = Field(
        title="Location of the fastq_pass e.g /my/path/fastq_pass/"
    )
    quantile: float = Field(default=0.90)
    min_consensus_depth: int = Field(default=20)
    min_quantile: float = Field(default=0.50)
    min_map_quality: int = Field(default=10)
    test_mode: Optional[bool] = Field(default=False)
    test_interval: Optional[int] = Field(default=1)
    test_samples: Optional[List[str]] = Field(default=None)
    test_start_time: Optional[int] = Field(default=0)
    test_end_time: Optional[int] = Field(default=1)
    test_start_time_query: Optional[int] = Field(default=0)
    test_end_time_query: Optional[int] = Field(default=1)

    class Config:
        validate_assignment = True
