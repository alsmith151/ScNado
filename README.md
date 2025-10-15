# ScNado

ScNado is a high-performance tool for barcode detection and processing in single-cell sequencing data.

## Features

- Fast barcode detection in FASTQ files
- Support for mismatches and slack regions
- Multi-threaded processing with Rayon
- Gzip compression support

## Installation

### Pre-built Binaries

Download the latest pre-built binaries from the [Releases](https://github.com/alsmith151/ScNado/releases) page:

- **macOS**: `scnado-macos-latest`
- **Ubuntu/Linux**: `scnado-ubuntu-latest`

After downloading, make the binary executable and optionally move it to your PATH:

```bash
chmod +x scnado-*
sudo mv scnado-* /usr/local/bin/scnado
```

### Building from Source

#### Prerequisites

- Rust 1.70 or later (install from [rustup.rs](https://rustup.rs/))

#### Build Instructions

1. Clone the repository:
```bash
git clone https://github.com/alsmith151/ScNado.git
cd ScNado
```

2. Build the project:
```bash
cargo build --release
```

3. The binary will be available at `target/release/scnado`

4. Optionally, install it to your PATH:
```bash
cargo install --path .
```

## Usage

### Find Barcodes

Detect barcodes in FASTQ files:

```bash
scnado find-barcodes \
  -i input.fastq.gz \
  -b barcodes.csv \
  -o output.fastq.gz \
  --n-missmatches 1
```

#### Options

- `-i, --input <FILE>`: Input FASTQ file (supports gzip compression)
- `-b, --barcodes <FILE>`: CSV file with columns 'barcode_type' and 'barcode_sequence'
- `-o, --output <FILE>`: Output FASTQ file
- `--slack-left <NUM>`: Number of bases to allow slack on the left (default: 0)
- `--slack-right <NUM>`: Number of bases to allow slack on the right (default: 0)
- `-n, --n-missmatches <NUM>`: Number of allowed mismatches (default: 0)
- `-v, --verbose <LEVEL>`: Verbosity level (default: 2)

#### Barcode File Format

The barcode file should be a CSV file with the following columns:

```csv
barcode_type,barcode_sequence
barcode1,ATCGATCG
barcode2,GCTAGCTA
```

## Examples

### Basic barcode detection
```bash
scnado find-barcodes -i input.fastq.gz -b barcodes.csv -o output.fastq.gz
```

### Allow up to 2 mismatches
```bash
scnado find-barcodes -i input.fastq.gz -b barcodes.csv -o output.fastq.gz --n-missmatches 2
```

### With slack regions
```bash
scnado find-barcodes -i input.fastq.gz -b barcodes.csv -o output.fastq.gz \
  --slack-left 2 --slack-right 2 --n-missmatches 1
```

## Development

### Running Tests

```bash
cargo test
```

### Building for Release

```bash
cargo build --release
```

## License

[Add your license here]

## Citation

[Add citation information if applicable]

## Contact

For questions or issues, please open an issue on the [GitHub repository](https://github.com/alsmith151/ScNado/issues).
