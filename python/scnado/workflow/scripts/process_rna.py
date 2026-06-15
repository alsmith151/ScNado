from pathlib import Path

import scanpy as sc
from scnado.rna import process_rna

# Snakemake variables
# Each entry in snakemake.input.counts is a Solo.out directory from one sublibrary
solo_out_dirs = snakemake.input.counts
output_h5ad = snakemake.output.h5ad
n_top_genes = snakemake.config.get("rna", {}).get("n_top_genes", 2000)

# STARsolo GeneFull_Ex50pAS filtered matrix is at:
# {Solo.out}/GeneFull_Ex50pAS/filtered/
counts_dirs = [
    str(Path(d) / "GeneFull_Ex50pAS" / "filtered")
    for d in solo_out_dirs
]

process_rna(counts_dirs, output_h5ad, n_top_genes=n_top_genes)
