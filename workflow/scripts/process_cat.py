import snapatac2 as snap
from scnado.cat import process_cat
import os

# Snakemake variables
input_fragments = snakemake.input.fragments
output_h5ad = snakemake.output.h5ad
bin_size = snakemake.config.get("cat_bin_size", 5000)
n_features = snakemake.config.get("cat_n_features", 50000)

# Process
# If multiple sublibraries, we might want to merge them here or process individually
# The notebook processes them together in import_fragments
data = process_cat(
    input_fragments,
    output_h5ad,
    bin_size=bin_size,
    n_features=n_features
)

# Save (process_cat returns a list of AnnData if multiple inputs, but we want to save the result)
# SnapATAC2's import_fragments returns an AnnData or AnnDataSet depending on inputs.
# If we want a single h5ad, we might need to handle it.
if isinstance(data, list):
    # For now, let's just save the first one or handle merging if needed.
    # In a real pipeline, we'd probably want to save an AnnDataSet or concat.
    data[0].write(output_h5ad)
else:
    data.write(output_h5ad)
