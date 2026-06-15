import snapatac2 as snap
from scnado.cat import export_coverage
import os

# Snakemake variables
h5ad = snakemake.input.h5ad
out_dir = snakemake.output[0]

# Load data
data = snap.read(h5ad)

# Export
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

export_coverage(data, out_dir)
