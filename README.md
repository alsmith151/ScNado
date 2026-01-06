# scnado

scnado is a high-performance pipeline for processing single-cell CUT&TAG (CAT) and RNA-seq data. It combines a fast Rust-based barcode detection engine with a Python-based analysis suite and a containerized Snakemake workflow.

## Features

- **Fast Barcode Detection**: Rust-based engine for high-throughput barcode matching and UMI extraction.
- **Mixed-Language Package**: Seamless integration of Rust performance and Python's single-cell ecosystem (Scanpy, SnapATAC2, Muon).
- **Containerized Workflow**: Snakemake pipeline with full support for Docker/Singularity and GitHub Container Registry (GHCR).
- **Multi-modal Integration**: Built-in support for linking CAT and RNA data using `muon` (MuData).

## Project Structure

| Directory | Description |
|-----------|-------------|
| [python/scnado/](python/scnado/) | Python package for single-cell analysis |
| [src/](src/) | Rust source code for barcode processing |
| [workflow/](workflow/) | Snakemake workflow and scripts |
| [config/](config/) | Configuration files and sample sheets |

## Installation

### From Source (Recommended)

The package is built using `maturin`. To install locally:

```bash
pip install .
```

This installs:
- **`scnado`** - Python CLI for the complete analysis pipeline
- **`scnado`** Python module - for programmatic access

The Rust functions (barcode detection, fragment generation) are accessible through the Python CLI and module:

```python
import scnado
barcode_path = scnado.get_barcode_csv()
```

### Docker

The project is containerized and available via GitHub Container Registry:

```bash
docker pull ghcr.io/alsmith151/scnado:latest
```

## Development and Testing

### Running Tests

The project includes both Rust and Python test suites.

**Python Tests**:
```bash
pytest
```

**Rust Tests**:
```bash
cargo test
```

Note: Rust tests are run without the `extension-module` feature by default to avoid linker issues with Python symbols.

### Building the Extension

To build the Python extension module for development:

```bash
maturin develop
```

## Usage

### Quick Start

1. **Initialize a new project**:
   ```bash
   scnado workflow init --outdir ./my_analysis
   cd ./my_analysis
   ```

2. **Configure your analysis**:
   - Edit `config/config.yaml` to set reference paths and container image
   - Edit `config/samples.tsv` with your sample information

3. **Run the pipeline**:
   ```bash
   scnado workflow run --cores 8
   ```

### Workflow Management

The `scnado workflow` commands manage the Snakemake pipeline:

1. **Initialize a new project**:
   ```bash
   scnado workflow init --outdir ./my_project
   ```

2. **Run the pipeline**:
   ```bash
   # Run with default settings (8 cores)
   scnado workflow run
   
   # Run with custom core count
   scnado workflow run --cores 16
   
   # Run with Docker containers
   scnado workflow run --use-docker --cores 8
   
   # Dry run to see what would be executed
   scnado workflow run --dry-run
   
   # Use conda for dependencies
   scnado workflow run --use-conda
   ```

### Core Commands

scnado provides direct access to the Rust-based processing tools:

**Find barcodes in FASTQ files**:
```bash
scnado find-barcodes \
  --r1 input_R1.fastq.gz \
  --r2 input_R2.fastq.gz \
  --barcodes barcodes.csv \
  --output-prefix results/barcoded/sample1 \
  --n-missmatches 2 \
  --enable-n-to-match
```

**Extract fragments from BAM file**:
```bash
scnado fragments \
  --bam aligned.bam \
  --output fragments.tsv.gz \
  --barcode-regex '^[^|]+\|(?P<barcode>[ACGT-]+)$' \
  --shift-plus 4 \
  --shift-minus -5
```

### Docker Usage

To run the entire pipeline from within the Docker container:

```bash
docker run -v $(pwd):/work -w /work ghcr.io/alsmith151/scnado:latest \
  scnado workflow init --outdir ./my_analysis

cd ./my_analysis

docker run -v $(pwd):/work -w /work ghcr.io/alsmith151/scnado:latest \
  scnado workflow run --use-docker --cores 8
```

Or with Singularity:

```bash
singularity exec docker://ghcr.io/alsmith151/scnado:latest \
  scnado workflow run --use-singularity --cores 8
```

## Development

### Prerequisites

- Rust 1.83+
- Python 3.11+
- Maturin

### Building from Source

```bash
maturin develop
```

## License

[MIT](LICENSE)
