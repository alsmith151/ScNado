from scnado._scnado import find_barcodes

find_barcodes(
    snakemake.input.r1,
    snakemake.input.r2,
    snakemake.input.barcodes,
    snakemake.params.output_prefix,
    snakemake.params.n_mismatches,
    snakemake.params.enable_n_to_match,
    getattr(snakemake.params, "trim_r2", True),
)
