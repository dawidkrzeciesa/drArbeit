rule build_matrix_fusion:
    input:
        expand("results/arriba/tables/fusions/{sample}.tsv", sample=samples)
    output:
        matrix_fusions = "results/oncoprint/matrix/fusion_matrix_for_oncoprint.tsv"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/build_fusion_oncoprint_matrix.py"



rule build_temp_table_for_fusion_mates:
    input:
        tab = "results/arriba/tables/fusions/{sample}.tsv"
    output:
        temp_tab = temp("results/oncoprint/matrix/temp/{sample}.tsv")
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/build_fusion_mates_temp_tab.py"



rule build_matrix_fusion_mates:
    input:
        tab = expand("results/oncoprint/matrix/temp/{sample}.tsv", sample=samples)
    output:
        matrix_fus_mates = "results/oncoprint/matrix/fusion_mates_matrix_for_oncoprint.tsv"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/build_fusion_mates_oncoprint_matrix.py"



rule oncoprint_fusions:
    input:
        matrix = "results/oncoprint/matrix/fusion_matrix_for_oncoprint.tsv",
    output:
        report("results/oncoprint/fusions_oncoprint.pdf", category="oncoprints", caption="../report/oncoprint_fusions.rst")
    conda:
        "../envs/oncoprint.yaml"
    params:
        title = config["oncoprint"]["title"]
    script:
        "../scripts/oncoprint_fusion.R"



rule oncoprint_fusion_mates:
    input:
        matrix = "results/oncoprint/matrix/fusion_mates_matrix_for_oncoprint.tsv",
    output:
        report("results/oncoprint/fusion_mates_oncoprint.pdf", category="oncoprints", caption="../report/oncoprint_fusion_mates.rst")
    conda:
        "../envs/oncoprint.yaml"
    params:
        title = config["oncoprint"]["title"]
    script:
        "../scripts/oncoprint_fusion_mates.R"


