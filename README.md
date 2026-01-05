# ScNado

ScNado is a high-performance pipeline for processing single-cell CUT&TAG (CAT) and RNA-seq data. It combines a fast Rust-based barcode detection engine with a Python-based analysis suite and a containerized Snakemake workflow.

## Features

- **Fast Barcode Detection**: Rust-based engine for high-throughput barcode matching and UMI extraction.
- **Mixed-Language Package**: Seamless integration of Rust performance and Python's single-cell ecosystem (Scanpy, SnapATAC2, Muon).
- **Containerized Workflow**: Snakemake pipeline with full support for Docker/Singularity and GitHub Container Registry (GHCR).
- **Multi-modal Integration**: Built-in support for linking CAT and RNA data using `muon` (MuData).

## Project Structure

- [python/scnado/](python/scnado/): Python package for single-cell analysis.
- [src/](src/): Rust source code for the `scnado` CLI.
- [workflow/](workflow/): Snakemake workflow and scripts.
- [config/](config/): Configuration files and sample sheets.
- [initial_workflow/](initial_workflow/): Original prototype notebooks and scripts.

## Installation

### Python Package (Mixed Rust/Python)

The package is built using `maturin`. To install it locally:

```bash
pip install .
```

This will compile the Rust components and install the `scnado` CLI and `scnado` Python module. The package includes default barcode reference files which can be accessed via:

```python
import scnado
barcode_path = scnado.get_barcode_csv()
```

### Docker

The project is containerized and available via GitHub Container Registry:

```bash
docker pull ghcr.io/<username>/scnado:latest
```

## Usage

### Snakemake Pipeline

The primary way to run the pipeline is via Snakemake.

1. **Configure the pipeline**:
   Edit [config/config.yaml](config/config.yaml) to set your reference paths and parameters.
   Update [config/samples.tsv](config/samples.tsv) with your sample metadata and FASTQ paths.

2. **Run the pipeline**:
   ```bash
   snakemake --use-conda --use-singularity --cores 8
   ```

### Rust CLI

The `scnado` binary provides low-level tools for barcode processing:

```bash
# Find barcodes and extract UMIs
scnado find-barcodes \
  --r1 input_R1.fastq.gz \
  --r2 input_R2.fastq.gz \
  --barcodes barcodes.csv \
  --output-prefix results/barcoded/sample1

# Generate fragments from BAM
scnado fragments \
  --bam aligned.bam \
  --output fragments.tsv.gz \
  --barcode-regex '^[^|]+\|(?P<barcode>[ACGT-]+)$'
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
