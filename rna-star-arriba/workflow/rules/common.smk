# input files for rule cutadapt_pe
def get_input(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit]

    if pd.isna(unit["fq1"]):
        # SRA sample -> always paired-end!!!
        accession = unit["sra"]
        return expand("sra/{accession}_{read}.fastq", accession=accession, read=[1, 2])
    else:
        return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()



# adapters for rule cutadapt_pe 
def get_adapters(wildcards):
    adapters = units.loc[wildcards.sample].loc[wildcards.unit, "adapters"]
    if isinstance(adapters, str):
        return adapters
    return ""



# input files for rule merge_fastqs
def get_fastqs(wildcards):
    return expand("results/trimmed/{sample}/{unit}_{read}.fastq.gz", unit=units.loc[wildcards.sample, "unit_name"], sample=wildcards.sample, read=wildcards.read)



# input.fq1 -> rule star_mapping
def get_fastq_fq1_star(wildcards):
    
    # triming skip: true = use raw reads
    if config["trimming"]["skip"]:
        unit=units.loc[wildcards.sample]
        
        # check for empty sra -> ONLY local fqs OR sra Data! 
        # Don't mix same sample_name with local fastq's and sra!!!
        if pd.isna(unit["sra"]).any():
            return units.loc[wildcards.sample]["fq1"]
        else:
            accession = unit["sra"]
            return expand("sra/{accession}_1.fastq", accession=accession)
    
    # triming skip: false = trimming
    else:
        return "results/merged/{sample}_R1.fastq.gz"



# input.fq2 -> rule star_mapping
def get_fastq_fq2_star(wildcards):
    
    # triming skip: true = use raw reads
    if config["trimming"]["skip"]:
        unit=units.loc[wildcards.sample]
        
        # check for empty sra -> ONLY local OR sra Data! 
        # Don't mix same sample_name with local fastq's and sra!!!
        if pd.isna(unit["sra"]).any():
            return units.loc[wildcards.sample]["fq2"]
        else:
            accession = unit["sra"]
            return expand("sra/{accession}_2.fastq", accession=accession)
    
    # triming skip: false = trimming
    else:
        return "results/merged/{sample}_R2.fastq.gz"

