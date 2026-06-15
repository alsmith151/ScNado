from scnado.multiome import integrate_cat_rna

# Snakemake variables
cat_h5ad = snakemake.input.cat
rna_h5ad = snakemake.input.rna
output_h5mu = snakemake.output.h5mu

# Integrate
integrate_cat_rna(cat_h5ad, rna_h5ad, output_h5mu)
