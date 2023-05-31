rule get_genome:
    output:
        "refs/genome.fasta"
    params:
        species = "homo_sapiens",
        datatype = "dna",
        build = "GRCh38",
        release = config["Reference"]["release"]
    log:
        "logs/ref/get_genome.log"
    cache: False  # save space and time with between workflow caching (see Snakemake docs)
    wrapper:
        "0.73.0/bio/reference/ensembl-sequence"



rule get_annotation_gtf:
    output:
        "refs/annotation.gtf"
    params:
        species = "homo_sapiens",
        build = "GRCh38",
        release = config["Reference"]["release"],
        fmt = "gtf",
        flavor = "" # optional, e.g. chr_patch_hapl_scaff, see Ensembl FTP.
    log:
        "logs/ref/get_annotation.log"
    cache: False  # save space and time with between workflow caching (see Snakemake docs)
    wrapper:
        "0.73.0/bio/reference/ensembl-annotation"



rule get_annotation_gff3:
    output:
        "refs/annotation.gff3"
    params:
        species = "homo_sapiens",
        build = "GRCh38",
        release = config["Reference"]["release"],
        fmt = "gff3",
        flavor="" # optional, e.g. chr_patch_hapl_scaff, see Ensembl FTP.
    log:
        "logs/ref/get_annotation_gff3.log"
    cache: False  # save space and time with between workflow caching (see Snakemake docs)
    wrapper:
        "0.73.0/bio/reference/ensembl-annotation"



rule get_blacklist:
    input:
    output:
        blacklist = "refs/arriba_default/blacklist_hg38_GRCh38_v2.1.0.tsv.gz",
    conda:
        "../envs/arriba.yaml"
    threads: 1
    shell:
        "(find .snakemake/ -iname blacklist_hg38_GRCh38_v2.1.0.tsv.gz -type f | xargs cp -t refs/arriba_default)"



rule get_known_fusions:
    input:
    output:
        known_fusions = "refs/arriba_default/known_fusions_hg38_GRCh38_v2.1.0.tsv.gz",
    conda:
        "../envs/arriba.yaml"
    threads: 1
    shell:
        "(find .snakemake/ -iname known_fusions_hg38_GRCh38_v2.1.0.tsv.gz -type f | xargs cp -t refs/arriba_default)"



rule get_cytobands:
    input:
    output:
        cytobands = "refs/arriba_default/cytobands_hg38_GRCh38_v2.1.0.tsv",
    conda:
        "../envs/arriba.yaml"
    threads: 1
    shell:
        "(find .snakemake/ -iname cytobands_hg38_GRCh38_v2.1.0.tsv -type f | xargs cp -t refs/arriba_default)"



rule get_protein_domains:
    input:
    output:
        protein_domains = "refs/arriba_default/protein_domains_hg38_GRCh38_v2.1.0.gff3",
    conda:
        "../envs/arriba.yaml"
    threads: 1
    shell:
        "(find .snakemake/ -iname protein_domains_hg38_GRCh38_v2.1.0.gff3 -type f | xargs cp -t refs/arriba_default)"



rule get_arriba_draw_fusion_script:
    input:
    output:
        R_script = "workflow/scripts/draw_fusions.R",
    conda:
        "../envs/arriba.yaml"
    threads: 1
    shell:
        "(find .snakemake/ -iname draw_fusions.R -type f | xargs cp -t workflow/scripts)"


