import importlib.resources as pkg_resources
import re
import subprocess
import sys
from pathlib import Path
from typing import Optional

import jinja2
import typer
from loguru import logger

from . import data as _data_pkg
from .cat import process_cat, analyze_cat_dataset
from .multiome import integrate_cat_rna
from .rna import process_rna

app = typer.Typer(help="scnado: Single-cell CUT&TAG and RNA-seq analysis pipeline")
workflow_app = typer.Typer(help="Workflow management")
app.add_typer(workflow_app, name="workflow")


@app.command()
def find_barcodes(
    r1: str = typer.Option(..., "--r1", help="Input R1 FASTQ file"),
    r2: str = typer.Option(..., "--r2", help="Input R2 FASTQ file"),
    barcodes: str = typer.Option(..., "--barcodes", help="Barcode CSV file"),
    output_prefix: str = typer.Option(..., "--output-prefix", help="Output file prefix"),
    n_missmatches: int = typer.Option(0, "--n-missmatches", help="Number of allowed mismatches"),
    enable_n_to_match: bool = typer.Option(False, "--enable-n-to-match", help="Allow N in barcode matching"),
) -> None:
    """Find barcodes in FASTQ files and extract UMIs."""
    try:
        from ._scnado import find_barcodes as rust_find_barcodes
        rust_find_barcodes(r1, r2, barcodes, output_prefix, n_missmatches, enable_n_to_match)
        logger.success("Barcode finding complete")
    except ImportError:
        logger.error("scnado Rust module not installed. Run: pip install .")
        raise typer.Exit(code=1)
    except Exception as e:
        logger.error(f"{e}")
        raise typer.Exit(code=1)


@app.command()
def fragments(
    bam: str = typer.Option(..., "--bam", help="Input BAM file"),
    output: str = typer.Option(..., "--output", help="Output fragments file"),
    barcode_regex: Optional[str] = typer.Option(None, "--barcode-regex", help="Regex to extract barcode from read name"),
    umi_regex: Optional[str] = typer.Option(None, "--umi-regex", help="Regex to extract UMI from read name"),
    shift_plus: int = typer.Option(0, "--shift-plus", help="Shift for plus strand"),
    shift_minus: int = typer.Option(0, "--shift-minus", help="Shift for minus strand"),
    fragment_length_extension: Optional[int] = typer.Option(None, "--fragment-length-extension", help="Fragment length extension"),
) -> None:
    """Extract fragments from BAM file."""
    if barcode_regex:
        try:
            re.compile(barcode_regex)
        except re.error as e:
            logger.error(f"Invalid barcode regex pattern: {e}")
            raise typer.Exit(code=1)

    if umi_regex:
        try:
            re.compile(umi_regex)
        except re.error as e:
            logger.error(f"Invalid umi regex pattern: {e}")
            raise typer.Exit(code=1)

    try:
        from ._scnado import fragments as rust_fragments
        rust_fragments(bam, output, barcode_regex, umi_regex, shift_plus, shift_minus, fragment_length_extension)
        logger.success("Fragment generation complete")
    except ImportError:
        logger.error("scnado Rust module not installed. Run: pip install .")
        raise typer.Exit(code=1)
    except Exception as e:
        logger.error(f"{e}")
        raise typer.Exit(code=1)


def _load_seqnado_genomes() -> dict:
    """Load genome configs from the seqnado registry (~/.config/seqnado/genome_config.json)."""
    import json
    registry = Path.home() / ".config" / "seqnado" / "genome_config.json"
    if not registry.exists():
        return {}
    with open(registry) as f:
        return json.load(f)


def _render_config(template_vars: dict, outfile: Path) -> None:
    template_path = str(pkg_resources.files(_data_pkg).joinpath("config_template.yaml.j2"))
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(str(Path(template_path).parent)))
    tmpl = env.get_template("config_template.yaml.j2")
    outfile.parent.mkdir(parents=True, exist_ok=True)
    outfile.write_text(tmpl.render(**template_vars))


@workflow_app.command()
def init(
    outdir: str = typer.Option(
        ".",
        "--outdir",
        "-o",
        help="Project root directory (default: current directory)",
    ),
) -> None:
    """Interactively create config/config.yaml and config/samples.tsv."""
    from seqnado.config.user_input import get_user_input

    typer.echo("=== scnado project initialisation ===\n")

    # Genome selection — try the seqnado registry first
    genomes = _load_seqnado_genomes()
    if not genomes:
        typer.echo(
            "No genomes found in ~/.config/seqnado/genome_config.json.\n"
            "Run 'seqnado init' first to register your genomes, then re-run this command.\n",
            err=True,
        )
        raise typer.Exit(code=1)

    available = list(genomes.keys())
    genome_name = get_user_input(
        "Genome?",
        choices=available,
        default=available[0],
    )
    gdata = genomes[genome_name]

    bowtie2_prefix = gdata.get("bt2_index", "")
    gtf = gdata.get("gtf", None)
    fasta = gdata.get("fasta", None)
    chromosome_sizes = gdata.get("chromosome_sizes", None)
    blacklist = gdata.get("blacklist", None)

    samples_path = get_user_input("Samples TSV path?", default="config/samples.tsv")
    barcode_mismatches = int(get_user_input("Barcode mismatches?", default="2"))
    barcode_csv = get_user_input("Custom barcode CSV? (leave blank for default)", required=False) or None

    enable_cat = get_user_input("Enable CAT analysis?", default="yes", is_boolean=True)
    enable_rna = get_user_input("Enable RNA analysis?", default="no", is_boolean=True)

    cat_bin_size = int(get_user_input("CAT bin size?", default="5000"))
    cat_n_features = int(get_user_input("CAT n features?", default="50000"))
    rna_n_top_genes = int(get_user_input("RNA n top genes?", default="2000"))

    scnado_container = get_user_input(
        "Container image?",
        default="oras://ghcr.io/alsmith151/scnado:latest",
    )

    template_vars = dict(
        samples=samples_path,
        barcode_mismatches=barcode_mismatches,
        barcode_csv=barcode_csv,
        genome_name=genome_name,
        bowtie2_prefix=bowtie2_prefix,
        gtf=gtf,
        fasta=fasta,
        chromosome_sizes=chromosome_sizes,
        blacklist=blacklist,
        cat_bin_size=cat_bin_size,
        cat_n_features=cat_n_features,
        rna_n_top_genes=rna_n_top_genes,
        enable_cat=enable_cat,
        enable_rna=enable_rna,
        enable_integration=None,
        scnado_container=scnado_container,
    )

    config_file = Path(outdir) / "config" / "config.yaml"
    _render_config(template_vars, config_file)
    logger.success(f"Created {config_file}")

    samples_file = Path(outdir) / "config" / "samples.tsv"
    if not samples_file.exists():
        samples_file.write_text(
            "sample\tsublibrary\tfastq_r1\tfastq_r2\tmodality\n"
            "sample1\tsublib1\tfastq_files/sample1_sublib1_R1.fastq.gz\t"
            "fastq_files/sample1_sublib1_R2.fastq.gz\tmixed\n"
        )
        logger.success(f"Created {samples_file}")

    typer.echo(f"\nNext steps:")
    typer.echo(f"  1. Edit {samples_file} with your sample information")
    typer.echo(f"  2. Run 'scnado workflow run --profile <profile>' to execute the pipeline")
    typer.echo(f"     Available profiles (from seqnado init): ls, lc, ss, sc ...")


@workflow_app.command()
def run(
    workflow_dir: Optional[str] = typer.Option(
        None,
        "--workflow-dir",
        "-d",
        help="Project root directory (default: current directory)",
    ),
    cores: int = typer.Option(
        8,
        "--cores",
        "-c",
        help="Number of cores to use",
    ),
    profile: Optional[str] = typer.Option(
        None,
        "--profile",
        "-p",
        help="SeqNado Snakemake profile shortcode, e.g. ls (local+singularity), ss (slurm+singularity)",
    ),
    dry_run: bool = typer.Option(False, "--dry-run", "-n", help="Show what would be executed"),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Quiet output"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Verbose output"),
    snakemake_args: Optional[str] = typer.Option(
        None,
        "--snakemake-args",
        help="Additional arguments to pass to snakemake (space-separated)",
    ),
) -> None:
    """Run the Snakemake pipeline."""
    base_dir = Path(workflow_dir) if workflow_dir else Path.cwd()

    snakefile = base_dir / "workflow" / "Snakefile"
    if not snakefile.exists():
        logger.error(f"Snakefile not found at {snakefile}. Make sure you're in the project root or specify --workflow-dir")
        raise typer.Exit(code=1)

    config_file = base_dir / "config" / "config.yaml"
    if not config_file.exists():
        logger.error(f"config.yaml not found at {config_file}. Run 'scnado workflow init' first.")
        raise typer.Exit(code=1)

    cmd = [
        "snakemake",
        "-s", str(snakefile),
        "--configfile", str(config_file),
        "--cores", str(cores),
    ]

    if profile:
        try:
            from seqnado.utils import get_preset_profiles, resolve_profile_path
            import importlib.resources as _pkg_resources
            import seqnado as _seqnado_pkg

            pkg_trav = _pkg_resources.files(_seqnado_pkg)
            profile_path = resolve_profile_path(profile, pkg_trav)
            if profile_path is None:
                available = list(get_preset_profiles().keys())
                logger.error(f"Unknown profile '{profile}'. Available: {', '.join(available)}")
                raise typer.Exit(code=1)
            cmd += ["--profile", str(profile_path)]
        except ImportError:
            logger.warning("seqnado not importable; passing --profile value directly to snakemake")
            cmd += ["--profile", profile]

    if dry_run:
        cmd.append("-n")
    if quiet:
        cmd.append("-q")
    if verbose:
        cmd.append("-v")
    if snakemake_args:
        cmd.extend(snakemake_args.split())

    logger.info(f"Running: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, cwd=str(base_dir))
        raise typer.Exit(code=result.returncode)
    except FileNotFoundError:
        logger.error("snakemake not found. Install it with: pip install snakemake")
        raise typer.Exit(code=1)


def main() -> None:
    """Entry point for the CLI."""
    app()


if __name__ == "__main__":
    main()
