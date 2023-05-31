rule get_sra:
    output:
        temp("sra/{accession}_1.fastq"),
        temp("sra/{accession}_2.fastq")
    log:
        "logs/get-sra/{accession}.log"
    threads: 8
    wrapper:
        "0.73.0/bio/sra-tools/fasterq-dump"



rule trimming:
    input:
        get_input
    output:
        fastq1 = temp("results/trimmed/{sample}/{unit}_R1.fastq.gz"),
        fastq2 = temp("results/trimmed/{sample}/{unit}_R2.fastq.gz"),
        qc = "results/trimmed/{sample}/{unit}.paired.qc.txt"
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    params:
        extra = config["trimming"]["cutadapt_params"],
        adapters = get_adapters,
    threads: 8
    wrapper:
        "0.73.0/bio/cutadapt/pe"



rule merge_fastqs:
    input:
        get_fastqs
    output:
        temp("results/merged/{sample}_{read}.fastq.gz")
    log:
        "logs/merge-fastqs/{sample}_{read}.log"
    wildcard_constraints:
        read = "R1|R2"
    shell:
        "cat {input} > {output} 2> {log}"


