from scnado._scnado import fragments

fragments(
    snakemake.input.bam,
    snakemake.output.fragments,
    snakemake.params.barcode_regex,
    snakemake.params.umi_regex,
    snakemake.params.shift_plus,
    snakemake.params.shift_minus,
    None,
)
