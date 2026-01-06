import subprocess
import sys
from pathlib import Path
from typing import Optional

import typer

from .cat import process_cat, analyze_cat_dataset
from .multiome import integrate_cat_rna
from .rna import process_rna

app = typer.Typer(help="pyscnado: Single-cell CUT&TAG and RNA-seq analysis pipeline")
workflow_app = typer.Typer(help="Workflow management")
app.add_typer(workflow_app, name="workflow")


@workflow_app.command()
def init(
    outdir: str = typer.Option(
        "scnado_project",
        "--outdir",
        "-o",
        help="Output directory for project files",
    ),
) -> None:
    """Initialize a new scnado project with configuration templates."""
    config_dir = Path(outdir) / "config"
    config_dir.mkdir(parents=True, exist_ok=True)

    # Create config.yaml template
    config_template = """samples: "config/samples.tsv"
# barcode_csv: "path/to/custom_barcodes.csv" # Optional: defaults to packaged barcodes
barcode_mismatches: 2

# Analysis flags - control which analyses to run
enable_cat: true          # Run CUT&TAG analysis
enable_rna: true          # Run RNA analysis
# enable_integration: true # Optional: defaults to true if both CAT and RNA enabled

# Reference paths
bowtie2_index: "reference/hg38"
minimap2_index: "reference/hg38.fa"
gtf_file: "reference/hg38.gtf"

# Containers
scnado_container: "docker://ghcr.io/alsmith151/scnado:latest"

# Tool parameters
cat_bin_size: 5000
cat_n_features: 50000
rna_n_top_genes: 2000
"""

    # Create samples.tsv template
    samples_template = """sample\tsublibrary\tfastq_r1\tfastq_r2\tmodality
sample1\tsublib1\tfastq_files/sample1_sublib1_R1.fastq.gz\tfastq_files/sample1_sublib1_R2.fastq.gz\tmixed
sample1\tsublib2\tfastq_files/sample1_sublib2_R1.fastq.gz\tfastq_files/sample1_sublib2_R2.fastq.gz\tmixed
"""

    config_file = config_dir / "config.yaml"
    samples_file = config_dir / "samples.tsv"

    with open(config_file, "w") as f:
        f.write(config_template)

    with open(samples_file, "w") as f:
        f.write(samples_template)

    typer.echo(f"✓ Created config template: {config_file}")
    typer.echo(f"✓ Created samples template: {samples_file}")
    typer.echo(f"\nNext steps:")
    typer.echo(f"1. Edit {config_file} with your reference paths and container")
    typer.echo(f"2. Edit {samples_file} with your sample information")
    typer.echo(f"3. Run 'pyscnado workflow run' to execute the pipeline")


@workflow_app.command()
def run(
    workflow_dir: Optional[str] = typer.Option(
        None,
        "--workflow-dir",
        "-d",
        help="Directory containing the workflow (default: current directory)",
    ),
    cores: int = typer.Option(
        8,
        "--cores",
        "-c",
        help="Number of cores to use",
    ),
    use_conda: bool = typer.Option(
        False,
        "--use-conda",
        help="Use conda for dependencies",
    ),
    use_singularity: bool = typer.Option(
        False,
        "--use-singularity",
        help="Use Singularity for containers",
    ),
    use_docker: bool = typer.Option(
        False,
        "--use-docker",
        help="Use Docker for containers",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help="Do a dry run (show what would be executed)",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet",
        "-q",
        help="Quiet output",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="Verbose output",
    ),
    snakemake_args: Optional[str] = typer.Option(
        None,
        "--snakemake-args",
        help="Additional arguments to pass to snakemake (space-separated)",
    ),
) -> None:
    """Run the Snakemake pipeline."""
    # Find the workflow directory
    base_dir = Path(workflow_dir) if workflow_dir else Path.cwd()

    # Check if Snakefile exists
    snakefile = base_dir / "workflow" / "Snakefile"
    if not snakefile.exists():
        typer.echo(
            f"Error: Snakefile not found at {snakefile}",
            err=True,
        )
        typer.echo(
            "Make sure you're in the project root directory or specify --workflow-dir",
            err=True,
        )
        raise typer.Exit(code=1)

    # Check if config exists
    config_file = base_dir / "config" / "config.yaml"
    if not config_file.exists():
        typer.echo(
            f"Error: config.yaml not found at {config_file}",
            err=True,
        )
        typer.echo(
            "Run 'pyscnado workflow init' to create configuration templates",
            err=True,
        )
        raise typer.Exit(code=1)

    # Build snakemake command
    cmd = [
        "snakemake",
        "-s",
        str(snakefile),
        "--configfile",
        str(config_file),
        "--cores",
        str(cores),
    ]

    if use_conda:
        cmd.append("--use-conda")

    if use_singularity:
        cmd.append("--use-singularity")

    if use_docker:
        cmd.append("--use-docker")

    if dry_run:
        cmd.append("-n")

    if quiet:
        cmd.append("-q")

    if verbose:
        cmd.append("-v")

    # Add any additional snakemake arguments
    if snakemake_args:
        cmd.extend(snakemake_args.split())

    typer.echo(f"Running: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, cwd=str(base_dir))
        raise typer.Exit(code=result.returncode)
    except FileNotFoundError:
        typer.echo(
            "Error: snakemake not found. Install it with: pip install snakemake",
            err=True,
        )
        raise typer.Exit(code=1)


def main() -> None:
    """Entry point for the CLI."""
    app()


if __name__ == "__main__":
    main()
