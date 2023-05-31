rule star_index:
    input:
        fasta = "refs/genome.fasta",
        ref_annot = "refs/annotation.gtf"
    output:
        directory("refs/STAR/index")
    threads: config["star_mapping"]["star_threads"]
    params:
        extra = "--sjdbGTFfile refs/annotation.gtf --sjdbOverhang {}".format(config["star_mapping"]["star_index_Overhang"])
    log:
        "logs/ref/star_index.log"
    cache: False  # save space and time with between workflow caching (see Snakemake docs)
    resources: mem_mb=40000
    wrapper:
        "0.74.0/bio/star/index"



rule star_mapping:
    input:
        fq1 = get_fastq_fq1_star,
        fq2 = get_fastq_fq2_star,
        ref_fa = rules.star_index.output,
        ref_gtf = "refs/annotation.gtf",
    output:
        # see STAR manual for additional output files
        protected("results/star/{sample}/Aligned.sortedByCoord.out.bam"),
        tab="results/star/{sample}/ReadsPerGene.out.tab"
    log:
        "logs/star/pe/{sample}.log"
    params:
        # path to STAR reference genome index
        index = rules.star_index.output,
        # optional parameters
        extra = "--sjdbGTFfile refs/annotation.gtf {}".format(config["star_mapping"]["star_params"])
    threads: config["star_mapping"]["star_threads"]
    resources: mem_mb=40000
    wrapper:
        "0.74.0/bio/star/align"



rule bam_index:
    input:
        "results/star/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        protected("results/star/{sample}/Aligned.sortedByCoord.out.bam.bai")
    params:
        "" # optional params string
    wrapper:
        "v1.17.3/bio/samtools/index"


rule gene_count_matrix:
    input:
        expand("results/star/{sample}/ReadsPerGene.out.tab", sample=samples)
    output:
        gene_count = "results/gene_count.tsv"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/merge_gene_count.py"