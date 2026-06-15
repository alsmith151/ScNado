# scnado

**scnado** is a pipeline for processing single-cell combinatorial indexing CUT&TAG (scCAT) and matched RNA-seq data. It combines a Rust-based barcode detection engine with a Python analysis suite and a containerised Snakemake workflow.

## Overview

scnado handles the full journey from raw FASTQ files to analysed AnnData/MuData objects:

1. **Barcode finding** — Rust engine identifies combinatorial barcodes in CAT reads and tags them into the read name
2. **Alignment** — Bowtie2 (CAT) and STARsolo (RNA) via containerised rules
3. **Fragment generation** — Rust engine converts aligned BAM to fragment files with Tn5 shift correction
4. **Downstream analysis** — SnapATAC2 (CAT) and Scanpy (RNA) processing scripts produce per-sample `.h5ad` files
5. **Multiome integration** — Optional Muon-based joint analysis producing `.h5mu` files

## Project Structure

| Path | Description |
|------|-------------|
| [python/scnado/](python/scnado/) | Python package and CLI |
| [src/](src/) | Rust source (barcode finding, fragment extraction) |
| [workflow/](workflow/) | Snakemake workflow |
| [workflow/rules/](workflow/rules/) | CAT and RNA Snakemake rule files |
| [workflow/scripts/](workflow/scripts/) | Python analysis scripts called by rules |

## Installation

Requires Python ≥ 3.11 and Rust (for building the extension module).

```bash
pip install .
```

This builds the Rust extension (`scnado._scnado`) via `maturin` and installs the `scnado` CLI entry point.

### Docker / Singularity

A pre-built image is available from GHCR:

```bash
docker pull ghcr.io/alsmith151/scnado:latest
# or for Singularity/Apptainer:
singularity pull oras://ghcr.io/alsmith151/scnado:latest
```

---

## Quick Start

### 1. Configure the pipeline

Run the interactive configuration wizard. It reads genome paths from your seqnado genome registry (`~/.config/seqnado/genome_config.json`) — run `seqnado init` first if you haven't already.

```bash
scnado workflow config --outdir ./my_project
```

This creates `my_project/config/config.yaml`.

### 2. Create the sample sheet

Point `design` at the directory containing your raw FASTQ files. It auto-classifies files as CAT or RNA by filename pattern and pairs sublibraries:

```bash
scnado workflow design /path/to/fastqs --output my_project/config/metadata.csv
```

The output is a CSV with columns `sample, sublibrary, cat_r1, cat_r2, rna_r1, rna_r2`.

### 3. Run the pipeline

```bash
cd my_project
scnado workflow run --cores 16
```

Use `--dry-run` to preview jobs without executing them. Use `--profile` to pass a seqnado Snakemake profile (e.g. `ss` for Slurm + Singularity):

```bash
scnado workflow run --profile ss --cores 32
```

---

## CLI Reference

### `scnado workflow config`

Interactively generate `config/config.yaml`.

```
Options:
  --outdir / -o   Project root directory [default: .]
```

Prompts for genome, alignment index paths, barcode settings, and which modalities (CAT / RNA / both) to enable.

---

### `scnado workflow design`

Scan a directory for FASTQ files and produce `config/metadata.csv`.

```
Arguments:
  FASTQ_DIR   Directory containing *.fastq.gz files [default: .]

Options:
  --output / -o         Output CSV path [default: config/metadata.csv]
  --cat-pattern         Regex to identify CAT files by name [default: (?i)cat]
  --rna-pattern         Regex to identify RNA files by name [default: (?i)rna]
  --pairing-pattern     Regex to extract sublibrary key from filename [default: \d+]
  --sample / -s         Assign all sublibraries to this sample name (skips prompt)
```

**Example:** If you have files like `scCAT1_R1.fastq.gz`, `scRNA1_R1.fastq.gz`, the default patterns will classify and pair them automatically.

---

### `scnado workflow run`

Execute the Snakemake pipeline.

```
Options:
  --workflow-dir / -d   Project root [default: current directory]
  --cores / -c          Threads [default: 8]
  --profile / -p        seqnado profile shortcode (e.g. ls, ss)
  --dry-run / -n        Show jobs without running
  --quiet / -q          Suppress output
  --verbose / -v        Verbose output
  --snakemake-args      Extra arguments forwarded to snakemake (space-separated)
```

---

### `scnado find-barcodes`

Low-level command (also called by the workflow). Finds combinatorial barcodes in a CAT FASTQ pair and writes barcoded output FASTQs.

```
Options:
  --r1                  Input R1 FASTQ
  --r2                  Input R2 FASTQ
  --barcodes            Barcode reference CSV
  --output-prefix       Output file prefix
  --n-missmatches       Allowed mismatches per barcode [default: 0]
  --enable-n-to-match   Treat N bases as matching any base
```

The barcode sequence is appended to the read name in the output FASTQ using a `|barcode|UMI` suffix that downstream tools (fragment extraction, STARsolo) parse via regex.

---

### `scnado fragments`

Extract a fragment file from an aligned BAM. Barcode and UMI are parsed from read names using Python regex named groups.

```
Options:
  --bam                       Input BAM file
  --output                    Output fragments.tsv.gz
  --barcode-regex             Regex with named group (?P<barcode>...)
  --umi-regex                 Regex with named group (?P<UMI>...)
  --shift-plus                Tn5 shift on plus strand [default: 0]
  --shift-minus               Tn5 shift on minus strand [default: 0]
  --fragment-length-extension Optional extension applied to fragment length
```

Default Tn5 shifts used by the workflow: `--shift-plus 4 --shift-minus -5`.

---

## Configuration File

`config/config.yaml` (generated by `scnado workflow config`):

```yaml
metadata: "config/metadata.csv"
barcode_mismatches: 2
# barcode_csv: "/path/to/custom_barcodes.csv"  # omit to use built-in

genome:
  name: hg38
  index:
    type: Bowtie2
    prefix: "/path/to/bowtie2/index"
  gtf: "/path/to/genes.gtf"
  fasta: "/path/to/genome.fa"
  chromosome_sizes: "/path/to/chrom.sizes"
  blacklist: "/path/to/blacklist.bed"

star_index: "/path/to/star/index"   # required when enable_rna: true

cat:
  bin_size: 5000
  n_features: 50000

rna:
  n_top_genes: 2000

enable_cat: true
enable_rna: true
enable_integration: true            # auto-enabled when both modalities active

scnado_container: "oras://ghcr.io/alsmith151/scnado:latest"
```

---

## Sample Sheet Format

`config/metadata.csv` — one row per sublibrary:

| sample | sublibrary | cat_r1 | cat_r2 | rna_r1 | rna_r2 |
|--------|-----------|--------|--------|--------|--------|
| treated | 1 | /data/scCAT1_R1.fastq.gz | /data/scCAT1_R2.fastq.gz | /data/scRNA1_R1.fastq.gz | /data/scRNA1_R2.fastq.gz |
| treated | 2 | /data/scCAT2_R1.fastq.gz | /data/scCAT2_R2.fastq.gz | /data/scRNA2_R1.fastq.gz | /data/scRNA2_R2.fastq.gz |

Multiple sublibraries belonging to the same `sample` are merged before the downstream analysis steps.

---

## Workflow Outputs

| Path | Description |
|------|-------------|
| `results/fragments/{sample}_{sublibrary}_fragments.tsv.gz` | Per-sublibrary fragment files |
| `results/processed_cat/{sample}.h5ad` | SnapATAC2 AnnData (CAT) |
| `results/coverage/{sample}/` | Per-cluster bigWig coverage tracks |
| `results/processed_rna/{sample}.h5ad` | Scanpy AnnData (RNA) |
| `results/integrated/{sample}.h5mu` | Muon MuData (joint CAT + RNA) |

---

## Development

### Prerequisites

- Python ≥ 3.11
- Rust ≥ 1.83
- `maturin`

### Build the extension

```bash
maturin develop
```

### Tests

```bash
pytest          # Python tests
cargo test      # Rust tests
```

Rust tests require omitting the `extension-module` feature to avoid Python linker symbols:

```bash
cargo test --no-default-features
```

---

## License

[MIT](LICENSE)
