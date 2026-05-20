# DORAS (Dynamically Optimized Reference for Adaptive Sampling)

DORAS is a tool designed to optimize enrichment of adaptive sampling (e.g., Oxford Nanopore Technologies' ReadUntil) by dynamically optimizing the length of the reference based on the read length distribution. 

This tool circumvents the limitations of Multi-Locus Sequence Typing (MLST) enrichment where the genomic context of the *loci* is unknown for clinical isolates. It works in two main phases: extending base MLST loci to capture context and then using the extended references for high-accuracy sequence typing.

## Table of Contents
- [Supported OS](#supported-os)
- [Requirements](#requirements)
- [Installation](#installation)
- [Project Structure](#project-structure)
- [Configuration](#configuration)
- [Usage](#usage)
- [Testing & Development](#testing--development)
- [Scripts](#scripts)
- [Citation](#citation)
- [License](#license)

## Supported OS
- **Linux** (Primary support)
- **macOS** (Tested on Apple Silicon)

## Requirements
- **Python**: >= 3.12.3
- **Tools**:
    - `minimap2` (>= 2.11)
    - `samtools`
    - `bcftools`
    - `seqkit`
    - `medaka` (for polishing)
    - `pyabpoa`

## Installation

### Pixi (Preferred)
[Pixi](https://pixi.sh/) is the recommended package manager. It handles both Python and system dependencies automatically.

1. Install Pixi:
   ```bash
   curl -fsSL https://pixi.sh/install.sh | sh
   ```
2. Clone and install:
   ```bash
   git clone <REPO_URL>
   cd doras
   pixi install
   ```

### Conda 
Alternatively, use the provided environment file:
```bash
git clone <REPO_URL>
cd doras
conda env create -f environment.yaml
conda activate doras
```

## Project Structure
```text
doras/
├── src/                # Core source code
│   ├── doras_manager.py   # Main orchestration logic
│   ├── doras_extension.py # Extension and trimming logic
│   ├── bigsdb_tools.py    # Integration with PubMLST/BigsDB
│   ├── config.py          # Configuration and Pydantic models
│   └── ...
├── base_refs/          # Initial MLST loci FASTA files
├── README.md            
├── tests/              # Unit and integration tests
├── main.py             # Primary CLI entry point
└── pixi.toml           # Pixi configuration and tasks
```

## Configuration
DORAS is configured via TOML files. You can generate or edit these manually or use the Marimo UI.

### Configuration Parameters

#### Root Level
- `experiment_name`: (String) A unique identifier for your experiment. This name is used for log files and output naming.
- `genome_size`: (Integer) The approximate genome size of the target organism (e.g., `5000000` for 5MB). Used to estimate required coverage and extension limits.

#### `[bigsdb]` - PubMLST/BigsDB Integration
- `base_api`: (URL) The base REST API endpoint for PubMLST (e.g., `https://rest.pubmlst.org`).
- `url`: (URL) The specific sequence definition URL for the target scheme.
- `scheme`: (String) The scheme path on PubMLST (e.g., `schemes/1/`).
- `db_selected`: (String) The database name on PubMLST (e.g., `pubmlst_escherichia_seqdef`).

#### `[mlst_genes_path]`
- `value`: (Path) Path to the local FASTA file containing the initial MLST loci sequences used as seeds for extension.

#### `[paths]`
- `fastq_files_path`: (Path) Directory where the sequencer/basecaller deposits FastQ files (e.g., `path/to/fastq_pass`).
- `output_dir`: (Path) Directory where results, temporary files, and logs will be stored.

#### `[sample_names]`
- `list`: (Array of Strings) List of barcode or sample names to process (e.g., `["barcode01", "barcode02"]`).

#### `[run_params]` - Core Algorithm Settings
- `quantile`: (Float, 0.0-1.0) The target read length quantile to use for reference extension. A higher value (e.g., `0.95`) means longer extensions based on the longer reads in the distribution.
- `min_quantile`: (Float, 0.0-1.0) The minimum quantile allowed during dynamic adjustment.
- `min_map_quality`: (Integer) Minimum mapping quality (MAPQ) score for a read to be considered for extension or typing.
- `min_consensus_depth`: (Integer) Minimum read depth required at a position to confidently call a consensus base during extension.

### Test Mode
Enable `test_mode` in the TOML to use local FastQ files instead of monitoring a directory for live sequencing. This is useful for simulating runs or re-running analysis on existing data.

```toml
[test_mode]
value = true
test_samples = ["path/to/test.fastq.gz"]
test_start_time = 0
test_end_time = 1
test_start_time_query = 0
test_end_time_query = 1
```
- `value`: (Boolean) Set to `true` to enable test mode.
- `test_samples`: (Array of Paths) List of paths to FastQ files to use for simulation.
- `test_start_time` / `test_end_time`: (Integer) The start and end time (in hours) relative to the beginning of the "virtual" run for the **Extension Phase**.
- `test_start_time_query` / `test_end_time_query`: (Integer) The start and end time (in hours) for the **Query Phase** simulation.

## Usage
DORAS is designed to run alongside a live sequencing run with basecalling enabled.

### Phase 1: Extension
Constructs an organism-specific extended reference.
```bash
# Using pixi
pixi run python main.py --toml config.toml --phase extension

# Or direct python
python main.py --toml config.toml --phase extension
```
Once finished, you will be prompted to concatenate the generated `final_ref.fasta` files for Phase 2.

### Phase 2: Query
Performs sequence typing against the extended reference.
```bash
python main.py --toml config.toml --phase query
```

## Testing & Development

### Automated Tasks
Use `pixi run` to execute predefined tasks:
- `test_query`: Runs a full simulation from extension to query.
- `test_ext`: Runs the extension phase on test data.
- `test_long`: Runs an extended test suite.
- `clean_up`: Removes temporary test directories.


## Citation
TODO: Add citation information when available.

## License
TODO: Add license information (e.g., MIT, GPL).RAS