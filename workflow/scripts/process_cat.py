from pathlib import Path

from scnado.cat import process_cat
from scnado.config import ScnadoConfig

CONFIG = ScnadoConfig.from_yaml("config/config.yaml")

fragment_files = list(snakemake.input.fragments)
output_h5ad = snakemake.output.h5ad

data = process_cat(
    fragment_files,
    output_h5ad,
    bin_size=CONFIG.cat.bin_size,
    n_features=CONFIG.cat.n_features,
)

data.write(output_h5ad)
