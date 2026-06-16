from scnado import get_barcode_whitelists
from seqnado.workflow.helpers.common import define_memory_requested, define_time_requested


rule find_barcodes_rna:
    input:
        r1=get_rna_r1,
        r2=get_rna_r2,
        barcodes=get_barcode_csv,
    output:
        r1=temp("results/barcoded_rna/{sample}_{sublibrary}_R1.fastq.gz"),
        r2=temp("results/barcoded_rna/{sample}_{sublibrary}_R2.fastq.gz"),
    log:
        "logs/find_barcodes_rna/{sample}_{sublibrary}.log",
    benchmark:
        "benchmarks/find_barcodes_rna/{sample}_{sublibrary}.tsv"
    message:
        "Finding barcodes for RNA {wildcards.sample} / {wildcards.sublibrary}"
    resources:
        mem=lambda wc, attempt: define_memory_requested(initial_value=1, attempts=attempt),
        runtime=lambda wc, attempt: define_time_requested(initial_value=4, attempts=attempt),
    threads: 8
    params:
        output_prefix="results/barcoded_rna/{sample}_{sublibrary}",
        n_mismatches=CONFIG.barcode_mismatches,
        enable_n_to_match=True,
        trim_r2=False,
    script:
        "../scripts/find_barcodes.py"


rule trim_rna:
    input:
        r1="results/barcoded_rna/{sample}_{sublibrary}_R1.fastq.gz",
    output:
        r1=temp("results/trimmed_rna/{sample}_{sublibrary}_trimmed.fq.gz"),
    log:
        "logs/trim_rna/{sample}_{sublibrary}.log",
    benchmark:
        "benchmarks/trim_rna/{sample}_{sublibrary}.tsv"
    message:
        "Trimming RNA reads for {wildcards.sample} / {wildcards.sublibrary}"
    container:
        "docker://ghcr.io/felixkrueger/trimgalore:latest"
    threads: 4
    resources:
        mem=lambda wc, attempt: define_memory_requested(initial_value=4, attempts=attempt),
        runtime=lambda wc, attempt: define_time_requested(initial_value=2, attempts=attempt),
    params:
        outdir="results/trimmed_rna/",
    shell:
        """
        trim_galore --poly_g -j {threads} --2colour 20 --nextera \
          --basename {wildcards.sample}_{wildcards.sublibrary} \
          --output_dir {params.outdir} \
          {input.r1} \
          > {log} 2>&1
        """


rule align_rna:
    input:
        r1="results/trimmed_rna/{sample}_{sublibrary}_trimmed.fq.gz",
        r2="results/barcoded_rna/{sample}_{sublibrary}_R2.fastq.gz",
        barcodes=lambda wc: get_barcode_whitelists(),
    output:
        bam=temp("results/star_rna/{sample}_{sublibrary}/Aligned.sortedByCoord.out.bam"),
        counts=temp(directory("results/star_rna/{sample}_{sublibrary}/Solo.out")),
    log:
        "logs/align_rna/{sample}_{sublibrary}.log",
    benchmark:
        "benchmarks/align_rna/{sample}_{sublibrary}.tsv"
    message:
        "Aligning RNA reads for {wildcards.sample} / {wildcards.sublibrary}"
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    threads: 12
    resources:
        mem=lambda wc, attempt: define_memory_requested(initial_value=40, attempts=attempt),
        runtime=lambda wc, attempt: define_time_requested(initial_value=8, attempts=attempt),
    params:
        star_index=CONFIG.star_index,
        outdir="results/star_rna/{sample}_{sublibrary}/",
    shell:
        """
        STAR --readFilesIn {input.r1} {input.r2} \
          --outFileNamePrefix {params.outdir} \
          --outTmpDir /tmp/star_{wildcards.sample}_{wildcards.sublibrary} \
          --genomeDir {params.star_index} \
          --readFilesCommand zcat \
          --runThreadN {threads} \
          --soloType CB_UMI_Complex \
          --soloCBwhitelist {input.barcodes} \
          --soloBarcodeReadLength 0 \
          --soloCBposition 1_0_1_7 1_38_1_45 1_76_1_83 1_114_1_118 \
          --soloUMIposition 1_119_1_126 \
          --soloCBmatchWLtype EditDist_2 \
          --soloFeatures GeneFull \
          --soloMultiMappers EM \
          --outSAMtype BAM SortedByCoordinate \
          --clipAdapterType CellRanger4 \
          --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
          > {log} 2>&1
        """


rule process_rna:
    input:
        counts=get_star_counts_for_sample,
    output:
        h5ad="results/processed_rna/{sample}.h5ad",
    log:
        "logs/process_rna/{sample}.log",
    benchmark:
        "benchmarks/process_rna/{sample}.tsv"
    message:
        "Processing RNA data for {wildcards.sample}"
    container:
        CONFIG.scnado_container
    resources:
        mem=lambda wc, attempt: define_memory_requested(initial_value=32, attempts=attempt),
        runtime=lambda wc, attempt: define_time_requested(initial_value=4, attempts=attempt),
    script:
        "../scripts/process_rna.py"
