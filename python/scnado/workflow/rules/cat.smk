from seqnado.workflow.helpers.common import define_memory_requested, define_time_requested


rule find_barcodes_cat:
    input:
        r1=get_cat_r1,
        r2=get_cat_r2,
        barcodes=get_barcode_csv,
    output:
        r1=temp("results/barcoded/{sample}_{sublibrary}_R1.fastq.gz"),
        r2=temp("results/barcoded/{sample}_{sublibrary}_R2.fastq.gz"),
    log:
        "logs/find_barcodes/{sample}_{sublibrary}.log",
    benchmark:
        "benchmarks/find_barcodes/{sample}_{sublibrary}.tsv"
    message:
        "Finding barcodes for {wildcards.sample} / {wildcards.sublibrary}"
    resources:
        mem=lambda wc, attempt: define_memory_requested(initial_value=1, attempts=attempt),
        runtime=lambda wc, attempt: define_time_requested(initial_value=4, attempts=attempt),
    threads: 8
    params:
        output_prefix="results/barcoded_cat/{sample}_{sublibrary}",
        n_mismatches=CONFIG.barcode_mismatches,
        enable_n_to_match=True,
        trim_r2=False,
    script:
        "../scripts/find_barcodes.py"


rule trim_cat:
    input:
        r1="results/barcoded_cat/{sample}_{sublibrary}_R1.fastq.gz",
    output:
        r1=temp("results/trimmed_cat/{sample}_{sublibrary}_trimmed.fq.gz"),
    log:
        "logs/trim_cat/{sample}_{sublibrary}.log",
    benchmark:
        "benchmarks/trim_cat/{sample}_{sublibrary}.tsv"
    message:
        "Trimming CAT reads for {wildcards.sample} / {wildcards.sublibrary}"
    container:
        "docker://ghcr.io/felixkrueger/trimgalore:latest"
    threads: 4
    resources:
        mem=lambda wc, attempt: define_memory_requested(initial_value=4, attempts=attempt),
        runtime=lambda wc, attempt: define_time_requested(initial_value=2, attempts=attempt),
    params:
        outdir="results/trimmed_cat/",
    shell:
        """
        trim_galore --poly_g -j {threads} --2colour 20 --nextera \
          --basename {wildcards.sample}_{wildcards.sublibrary} \
          --output_dir {params.outdir} \
          {input.r1} \
          > {log} 2>&1
        """


rule align_cat:
    input:
        r1="results/trimmed_cat/{sample}_{sublibrary}_trimmed.fq.gz",
    output:
        bam=temp("results/aligned_cat/{sample}_{sublibrary}.bam"),
    log:
        "logs/align_cat/{sample}_{sublibrary}.log",
    benchmark:
        "benchmarks/align_cat/{sample}_{sublibrary}.tsv"
    message:
        "Aligning CAT reads for {wildcards.sample} / {wildcards.sublibrary}"
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    threads: 8
    resources:
        mem=lambda wc, attempt: define_memory_requested(initial_value=16, attempts=attempt),
        runtime=lambda wc, attempt: define_time_requested(initial_value=4, attempts=attempt),
    params:
        index=CONFIG.genome.index.prefix,
    shell:
        """
        bowtie2 -p {threads} -x {params.index} -U {input.r1} 2>> {log} | \
        samtools view -@ {threads} -b - > {output.bam}.unsorted.bam && \
        samtools sort -@ {threads} -o {output.bam} {output.bam}.unsorted.bam 2>> {log} && \
        rm {output.bam}.unsorted.bam
        """


rule make_fragments:
    input:
        bam="results/aligned_cat/{sample}_{sublibrary}.bam",
    output:
        fragments=temp("results/fragments/{sample}_{sublibrary}_fragments.tsv.gz"),
    log:
        "logs/make_fragments/{sample}_{sublibrary}.log",
    benchmark:
        "benchmarks/make_fragments/{sample}_{sublibrary}.tsv"
    message:
        "Extracting fragments for {wildcards.sample} / {wildcards.sublibrary}"
    resources:
        mem=lambda wc, attempt: define_memory_requested(initial_value=6, attempts=attempt),
        runtime=lambda wc, attempt: define_time_requested(initial_value=1, attempts=attempt),
    params:
        barcode_regex=r"^[^|]+\|(?P<barcode>[ACGT-]+)\|",
        umi_regex=r"|(?P<UMI>[ATGC]+)$",
        shift_plus=4,
        shift_minus=-5,
    script:
        "../scripts/make_fragments.py"


rule process_cat:
    input:
        fragments=get_fragments_for_sample,
    output:
        h5ad="results/processed_cat/{sample}.h5ad",
    log:
        "logs/process_cat/{sample}.log",
    benchmark:
        "benchmarks/process_cat/{sample}.tsv"
    message:
        "Processing CAT data for {wildcards.sample}"
    container:
        CONFIG.scnado_container
    resources:
        mem=lambda wc, attempt: define_memory_requested(initial_value=32, attempts=attempt),
        runtime=lambda wc, attempt: define_time_requested(initial_value=4, attempts=attempt),
    script:
        "../scripts/process_cat.py"


rule export_coverage:
    input:
        h5ad="results/processed_cat/{sample}.h5ad",
    output:
        directory("results/coverage/{sample}/"),
    log:
        "logs/export_coverage/{sample}.log",
    benchmark:
        "benchmarks/export_coverage/{sample}.tsv"
    message:
        "Exporting coverage for {wildcards.sample}"
    container:
        CONFIG.scnado_container
    resources:
        mem=lambda wc, attempt: define_memory_requested(initial_value=12, attempts=attempt),
        runtime=lambda wc, attempt: define_time_requested(initial_value=4, attempts=attempt),
    script:
        "../scripts/export_coverage.py"
