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
   git clone git@github.com:RametteLab/DORAS.git
   cd doras
   pixi install
   ```

### Conda (Coming soon)


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
- `url`: (URL) The specific sequence definition URL for the target scheme. this is the URL that you would use to access the scheme in the browser.
- `scheme`: (String) The scheme path on PubMLST (e.g., `schemes/1/`).
- `db_selected`: (String) The database name on PubMLST (e.g., `pubmlst_escherichia_seqdef`).

#### `[mlst_genes_path]`
- `value`: (Path) Path to the local FASTA file containing the initial MLST loci sequences used as seeds for extension.

These loci are used to construct the initial reference and are not modified during extension. They should be created by the user based on the MLST scheme used in the experiment.
see examples [here](https://github.com/RametteLab/DORAS/blob/main/base_refs) for more details.

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

#### Understanding Test Mode & Simulation
Test mode allows you to simulate the dynamic behavior of DORAS using existing data. It is particularly useful for optimizing parameters or validating results without needing a live sequencing run.

**How it works:**
1. **Historical Data:** It takes existing FastQ files (defined in `test_samples`) that contain reads with ONT-standard timestamps in their headers.
2. **Read Extraction by Time:** DORAS parses these timestamps to determine the relative start time of the experiment. It then extracts reads into "virtual" time windows (controlled by `test_interval`).
3. **Simulation Phases:**
   - **Extension Simulation:** Reads from `test_start_time` to `test_end_time` are progressively "released" into the `fastq_files_path`. This simulates the first few hours of a run where the reference is being extended.
   - **Query Simulation:** Once the reference is ready (or if running in query phase), it simulates the typing process using reads from `test_start_time_query` to `test_end_time_query`.
4. **Mock File Creation:** In the background, a `FastQTimeExtractor` creates mock FastQ files in the input directory at regular intervals, mimicking how MinKNOW/Guppy deposits files during a real run.

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