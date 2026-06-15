from seqnado.workflow.helpers.common import define_memory_requested, define_time_requested


rule integrate:
    input:
        cat="results/processed_cat/{sample}.h5ad",
        rna="results/processed_rna/{sample}.h5ad",
    output:
        h5mu="results/integrated/{sample}.h5mu",
    log:
        "logs/integrate/{sample}.log",
    benchmark:
        "benchmarks/integrate/{sample}.tsv"
    message:
        "Integrating CAT and RNA data for {wildcards.sample}"
    container:
        CONFIG.scnado_container
    resources:
        mem=lambda wc, attempt: define_memory_requested(initial_value=32, attempts=attempt),
        runtime=lambda wc, attempt: define_time_requested(initial_value=2, attempts=attempt),
    script:
        "../scripts/integrate.py"
