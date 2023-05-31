rule arriba_fusion_calling:
    input:
        bam = "results/star/{sample}/Aligned.sortedByCoord.out.bam",
        bam_index = "results/star/{sample}/Aligned.sortedByCoord.out.bam.bai",
        # path to reference genome
        genome = "refs/genome.fasta",
        # path to annotation gtf
        annotation = "refs/annotation.gtf",
        blacklist = (rules.get_blacklist.output if config["arriba"]["default_lists"] else config["arriba"]["blacklist"]),
        known_fusions = (rules.get_known_fusions.output if config["arriba"]["default_lists"] else config["arriba"]["known_fusions"]),
        tag = (rules.get_known_fusions.output if config["arriba"]["default_lists"] else config["arriba"]["tag"])
    output:
        # gene fusions
        fusions = "results/arriba/tables/fusions/{sample}.tsv",
        # discarded gene fusions through arriba filters
        discarded = "results/arriba/tables/discarded_fusions/{sample}.discarded.tsv"
    log:
        "logs/arriba/{sample}.log"
    conda:
        "../envs/arriba.yaml"
    params:
        arriba_params = config["arriba"]["arriba_params"]
    threads: 2
    resources: mem_mb=20000
    shell:
        "(arriba -x {input.bam} -g {input.annotation} -a {input.genome} -b {input.blacklist} -k {input.known_fusions} "
        "{params.arriba_params} -o {output.fusions} -O {output.discarded} -t {input.tag}) &>{log}"



rule arriba_fusion_drawing:
    input:
        fusions = "results/arriba/tables/fusions/{sample}.tsv",
        bam = "results/star/{sample}/Aligned.sortedByCoord.out.bam",
        bam_index = "results/star/{sample}/Aligned.sortedByCoord.out.bam.bai",
        annotation = "refs/annotation.gtf",
        cytobands = (rules.get_cytobands.output if config["arriba"]["default_lists"] else config["arriba"]["cytobands"]),
        protein_domains = (rules.get_protein_domains.output if config["arriba"]["default_lists"] else config["arriba"]["protein_domains"]),
        r_script = rules.get_arriba_draw_fusion_script.output
    output:
        pdf = "results/arriba/drawed_fusions/{sample}.pdf"
    log:
        "logs/arriba/{sample}.visualization.log"
    conda:
        "../envs/arriba.yaml"
    params:
        arriba_draw_params = config["arriba"]["arriba_draw_params"]
    threads: 2
    resources: mem_mb=20000
    shell:
        "(workflow/scripts/./draw_fusions.R --fusions={input.fusions} --alignments={input.bam} --output={output.pdf} --annotation={input.annotation} "
        "{params.arriba_draw_params} --cytobands={input.cytobands} --proteinDomains={input.protein_domains}) &>{log}"



rule merge_arriba_results:
    input:
        expand("results/arriba/tables/fusions/{sample}.tsv", sample=samples)
    output:
        merged_tsv = report("results/arriba/tables/merged_arriba_results.tsv", category="arriba fusions", caption="../report/arriba_fusion_calling.rst")
    conda:
        "../envs/pandas.yaml"
    params:
        "results/fusions"
    script:
        "../scripts/merge_arriba_results.py"

