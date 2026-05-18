# DO# DORAS (Dynamically Optimized Reference for Adaptive Sampling)

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
    - `bedtools`
    - `blast`
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
├── scripts/            # Auxiliary analysis scripts
├── base_refs/          # Initial MLST loci FASTA files
├── datasets/           # Sample datasets for verification
├── tests/              # Unit and integration tests
├── main.py             # Primary CLI entry point
├── marimo_ui.py        # Interactive UI for config generation
├── pixi.toml           # Pixi configuration and tasks
└── environment.yaml    # Conda environment definition
```

## Configuration
DORAS is configured via TOML files. You can generate or edit these manually or use the Marimo UI.

### Interactive Configuration
Run the Marimo notebook to use a GUI for generating your TOML configuration:
```bash
pixi run marimo edit marimo_ui.py
```

### TOML Key Sections
- `genome_size`: Approximate genome size of the target organism.
- `[bigsdb]`: API endpoints for PubMLST/BigsDB.
- `[mlst_genes_path]`: Path to the initial FASTA containing MLST loci.
- `[paths]`: Input FastQ and output directory locations.
- `[sample_names]`: List of barcodes/samples to process.
- `[run_params]`: Parameters for mapping quality, depth, and quantiles.

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

### Test Mode
Enable `test_mode` in the TOML to use local FastQ files instead of monitoring a directory for live sequencing:
```toml
[test_mode]
value = true
test_samples = ["path/to/test.fastq.gz"]
```

## Scripts
Additional tools are available in the `scripts/` directory:
- `query_st.py`: Query Sequence Type from PubMLST.
- `analysis_assemblies.py`: Analyze assembly results.
- `mapping_verif_pipeline_script.py`: Verification pipeline for mapping.

## Citation
TODO: Add citation information when available.

## License
TODO: Add license information (e.g., MIT, GPL).RAS