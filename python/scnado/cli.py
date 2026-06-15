import importlib.resources as pkg_resources
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Optional

import typer
from loguru import logger

from . import data as _data_pkg
from . import workflow as _workflow_pkg

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
    import jinja2
    template_path = str(pkg_resources.files(_data_pkg).joinpath("config_template.yaml.j2"))
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(str(Path(template_path).parent)))
    tmpl = env.get_template("config_template.yaml.j2")
    outfile.parent.mkdir(parents=True, exist_ok=True)
    outfile.write_text(tmpl.render(**template_vars))


def _prompt_config(outdir: str) -> None:
    """Interactively create config/config.yaml."""
    from seqnado.config.user_input import get_user_input

    genomes = _load_seqnado_genomes()
    if not genomes:
        typer.echo(
            "No genomes found in ~/.config/seqnado/genome_config.json.\n"
            "Run 'seqnado init' first to register your genomes, then re-run this command.\n",
            err=True,
        )
        raise typer.Exit(code=1)

    available = list(genomes.keys())
    genome_name = get_user_input("Genome?", choices=available, default=available[0])
    gdata = genomes[genome_name]

    bowtie2_prefix = gdata.get("bt2_index", "")
    gtf = gdata.get("gtf", None)
    fasta = gdata.get("fasta", None)
    chromosome_sizes = gdata.get("chromosome_sizes", None)
    blacklist = gdata.get("blacklist", None)

    metadata_path = get_user_input("Metadata CSV path?", default="config/metadata.csv")
    barcode_mismatches = int(get_user_input("Barcode mismatches?", default="2"))
    barcode_csv = get_user_input("Custom barcode CSV? (leave blank for default)", required=False) or None

    enable_cat = get_user_input("Enable CAT analysis?", default="yes", is_boolean=True)
    enable_rna = get_user_input("Enable RNA analysis?", default="no", is_boolean=True)

    star_index = None
    if enable_rna:
        star_index = get_user_input("STAR genome directory?", required=True) or None

    cat_bin_size = int(get_user_input("CAT bin size?", default="5000"))
    cat_n_features = int(get_user_input("CAT n features?", default="50000"))
    rna_n_top_genes = int(get_user_input("RNA n top genes?", default="2000"))

    scnado_container = get_user_input(
        "Container image?",
        default="oras://ghcr.io/alsmith151/scnado:latest",
    )

    template_vars = dict(
        metadata=metadata_path,
        barcode_mismatches=barcode_mismatches,
        barcode_csv=barcode_csv,
        genome_name=genome_name,
        bowtie2_prefix=bowtie2_prefix,
        gtf=gtf,
        fasta=fasta,
        chromosome_sizes=chromosome_sizes,
        blacklist=blacklist,
        star_index=star_index,
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


@workflow_app.command()
def config(
    outdir: str = typer.Option(
        ".",
        "--outdir",
        "-o",
        help="Project root directory (default: current directory)",
    ),
) -> None:
    """Interactively create config/config.yaml for the scnado pipeline."""
    typer.echo("=== scnado workflow config ===\n")
    _prompt_config(outdir)
    typer.echo("\nNext steps:")
    typer.echo("  1. Run 'scnado workflow design <fastq_dir>' to create config/metadata.csv")
    typer.echo("  2. Run 'scnado workflow run' to execute the pipeline")


@workflow_app.command()
def design(
    fastq_dir: str = typer.Argument(
        ".",
        help="Directory containing FASTQ files",
    ),
    output: str = typer.Option(
        "config/metadata.csv",
        "--output",
        "-o",
        help="Output metadata CSV path",
    ),
    cat_pattern: str = typer.Option(
        r"(?i)cat",
        "--cat-pattern",
        help="Regex to identify CAT FASTQ files by filename",
    ),
    rna_pattern: str = typer.Option(
        r"(?i)rna",
        "--rna-pattern",
        help="Regex to identify RNA FASTQ files by filename",
    ),
    pairing_pattern: str = typer.Option(
        r"\d+",
        "--pairing-pattern",
        help=(
            "Regex to extract the sublibrary key from each filename stem. "
            "The first capture group (or full match if no groups) is used to pair "
            "CAT and RNA files that belong to the same sublibrary. "
            "Default '\\d+' extracts the first number, e.g. 'scCAT1' and 'scRNA1' both yield '1'."
        ),
    ),
    sample: Optional[str] = typer.Option(
        None,
        "--sample",
        "-s",
        help="Assign all sublibraries to this sample name (skips interactive prompt)",
    ),
) -> None:
    """Generate config/metadata.csv linking CAT and RNA FASTQ files by sublibrary.

    Scans FASTQ_DIR for *.fastq.gz files, classifies them as CAT or RNA using
    --cat-pattern / --rna-pattern, then pairs them by the sublibrary key extracted
    with --pairing-pattern. Produces a wide-format CSV with one row per sublibrary:

        sample, sublibrary, cat_r1, cat_r2, rna_r1, rna_r2
    """
    import pandas as pd

    fastq_path = Path(fastq_dir).resolve()
    if not fastq_path.exists():
        logger.error(f"Directory not found: {fastq_path}")
        raise typer.Exit(code=1)

    fastq_files = sorted(fastq_path.glob("*.fastq.gz"))
    if not fastq_files:
        logger.error(f"No *.fastq.gz files found in {fastq_path}")
        raise typer.Exit(code=1)

    # Validate user-supplied patterns up front
    for name, pat in [("--cat-pattern", cat_pattern), ("--rna-pattern", rna_pattern), ("--pairing-pattern", pairing_pattern)]:
        try:
            re.compile(pat)
        except re.error as e:
            logger.error(f"Invalid regex for {name}: {e}")
            raise typer.Exit(code=1)

    # Classify each file: classified[sublib_key][modality][read_num] = Path
    classified: dict[str, dict[str, dict[int, Path]]] = defaultdict(
        lambda: {"cat": {}, "rna": {}}
    )
    skipped = []

    for f in fastq_files:
        stem = f.name.removesuffix(".gz").removesuffix(".fastq")

        # Extract read number from _R1/_R2 suffix
        read_match = re.search(r"_R([12])(?:_\d+)?$", stem)
        if read_match is None:
            logger.warning(f"Cannot detect read number in '{f.name}', skipping")
            skipped.append(f.name)
            continue
        read_num = int(read_match.group(1))

        # Classify modality
        if re.search(cat_pattern, stem):
            modality = "cat"
        elif re.search(rna_pattern, stem):
            modality = "rna"
        else:
            logger.warning(f"Cannot classify '{f.name}' as CAT or RNA, skipping")
            skipped.append(f.name)
            continue

        # Extract sublibrary key via pairing_pattern
        key_match = re.search(pairing_pattern, stem)
        if key_match is None:
            logger.warning(f"Pairing pattern '{pairing_pattern}' did not match '{f.name}', skipping")
            skipped.append(f.name)
            continue
        sublib_key = key_match.group(1) if key_match.lastindex else key_match.group(0)

        classified[sublib_key][modality][read_num] = f

    if not classified:
        logger.error("No FASTQ files could be classified. Check --cat-pattern, --rna-pattern, and --pairing-pattern.")
        raise typer.Exit(code=1)

    # Report detections
    typer.echo(f"\nScanned {fastq_path}")
    typer.echo(f"Found {len(fastq_files)} FASTQ file(s), {len(classified)} sublibrary key(s) detected:\n")
    for key in sorted(classified.keys()):
        cat = classified[key]["cat"]
        rna = classified[key]["rna"]
        cat_r1 = cat.get(1, Path("MISSING")).name
        cat_r2 = cat.get(2, Path("MISSING")).name
        rna_r1 = rna.get(1, Path("MISSING")).name
        rna_r2 = rna.get(2, Path("MISSING")).name
        typer.echo(f"  [{key}]  CAT: {cat_r1} / {cat_r2}   RNA: {rna_r1} / {rna_r2}")

    if skipped:
        typer.echo(f"\nSkipped {len(skipped)} file(s): {', '.join(skipped)}")

    typer.echo("")

    # Collect sample name assignments
    rows = []
    for key in sorted(classified.keys()):
        cat = classified[key]["cat"]
        rna = classified[key]["rna"]

        if sample:
            sname = sample
        else:
            sname = typer.prompt(f"Sample name for sublibrary '{key}'", default=f"sample_{key}")

        rows.append({
            "sample": sname,
            "sublibrary": key,
            "cat_r1": str(cat[1]) if 1 in cat else None,
            "cat_r2": str(cat[2]) if 2 in cat else None,
            "rna_r1": str(rna[1]) if 1 in rna else None,
            "rna_r2": str(rna[2]) if 2 in rna else None,
        })

    df = pd.DataFrame(rows)
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    logger.success(f"Metadata saved → {output_path}")
    typer.echo(f"\nNext step: run 'scnado workflow run' to execute the pipeline")


@workflow_app.command(hidden=True)
def init(
    outdir: str = typer.Option(
        ".",
        "--outdir",
        "-o",
        help="Project root directory (default: current directory)",
    ),
) -> None:
    """Deprecated: use 'scnado workflow config' instead."""
    logger.warning("'scnado workflow init' is deprecated. Use 'scnado workflow config' instead.")
    typer.echo("=== scnado project initialisation ===\n")
    _prompt_config(outdir)
    typer.echo("\nNext steps:")
    typer.echo("  1. Run 'scnado workflow design <fastq_dir>' to create config/metadata.csv")
    typer.echo("  2. Run 'scnado workflow run' to execute the pipeline")


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
    base_dir = Path(workflow_dir).resolve() if workflow_dir else Path.cwd()

    # Prefer a local workflow/Snakefile; fall back to the packaged one
    local_snakefile = base_dir / "workflow" / "Snakefile"
    if local_snakefile.exists():
        snakefile = local_snakefile
    else:
        snakefile = Path(str(pkg_resources.files(_workflow_pkg).joinpath("Snakefile")))
        if not snakefile.exists():
            logger.error("Snakefile not found locally or in the installed package.")
            raise typer.Exit(code=1)
        logger.info(f"Using packaged Snakefile: {snakefile}")

    config_file = base_dir / "config" / "config.yaml"
    if not config_file.exists():
        logger.error(f"config.yaml not found at {config_file}. Run 'scnado workflow config' first.")
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
