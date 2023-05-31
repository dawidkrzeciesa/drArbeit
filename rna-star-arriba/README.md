# Snakemake workflow: detection of gene fusions


[![Snakemake](https://img.shields.io/badge/snakemake-≥6.1.0-brightgreen.svg)](https://snakemake.bitbucket.io)

This workflow detects gene fusions from human paired-end RNA sequencing data (__Research Use Only__). 

Fusion calling  is done by [Arriba](https://github.com/suhrig/arriba) after alignment of RNA-seq-data to reference genome [GRCh38](https://www.ensembl.org/Homo_sapiens/Info/Index) with [STAR](https://github.com/alexdobin/STAR). Subsequently, an [OncoPrint](https://github.com/jokergoo/ComplexHeatmap) of gene fusions will be generated. 


This workflow integrates download of [SRA](https://www.ncbi.nlm.nih.gov/sra) samples (if required) and genome reference data from [Ensembl](https://www.ensembl.org/index.html).

Optional adapter trimming (with [cutadapt](https://github.com/marcelm/cutadapt)) and [drawing](https://arriba.readthedocs.io/en/latest/visualization/) of detected fusions can be performed.


## License

This workflow is distributed under [MIT License](https://github.com/dawidkrzeciesa/rna-star-arriba/blob/main/LICENSE).

## Usage

If you use this workflow in a paper/publication, don't forget to give credits to the authors by citing the URL of this (original) repository. Thank you in advance.

### Hardware recommendations

[STAR](https://github.com/alexdobin/STAR/tree/master/doc) neads about 30 GB free RAM memory per sample. 64 GB RAM is recommended for this workflow.

### Step 1: Obtain a copy of this workflow

[Clone](https://help.github.com/en/articles/cloning-a-repository) this repository to your local system, at the place where data analysis will be performed.

### Step 2: Configure workflow

* Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config.yaml` to configure the workflow execution, and `sample_sheet.tsv` to specify your sample setup. 
* You can run it on local and/or SRA paired-end RNA seq data. Specify at least one or more sequencing units for each sample in the `unit_name` column (`sample_sheet.tsv`) such as replicates, lanes or runs.
* Adapter trimming is not necessary, but if desired, adapter sequences for each sample have to be provided in the `sample_sheet.tsv`  according to the [cutadapt manual](https://cutadapt.readthedocs.io/en/stable/guide.html#trimming-paired-end-reads).

### Step 3: Install Snakemake

Install [miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

Open new terminal and create conda environment: 

    conda create -n snakemake

Activate conda enviroment, [set up channels](https://bioconda.github.io/user/install.html#set-up-channels) and install mamba:

    conda activate snakemake

    conda config --add channels defaults
    
    conda config --add channels bioconda
    
    conda config --add channels conda-forge

    conda install mamba

Install [Snakemake](https://snakemake.readthedocs.io/en/stable/):

    mamba install snakemake

### Step 4: Execute workflow

Within your terminal: Navigate into the rna-star-arriba directory and activate the snakemake conda environment:

    conda activate snakemake

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N --resources mem_mb=$M

using `$N` cores and `$M` mb of RAM memory (at least 40000) 

Snakemake workflows can be executed in a cluster as well. For more details, see the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).

### Step 5: View and examine results

After successful execution, you can also create a self-contained interactive HTML [report](https://snakemake.readthedocs.io/en/stable/snakefiles/reporting.html) with all results via:

    snakemake --report report_name.zip

This report can, e.g., be forwarded to your collaborators.

### Known issues:
* _Rule star_mapping_ rises an error if the number of provided cores/threads is too high (That was the case for >20 threads per sample on our workstation).
* Generally, a single sample can be analyzed, but OncoPrint needs more than one sample and one gene fusion to generate a plot. Otherwise, _rule oncoprint_fusions_ rises an error. However, the Arriba results are saved into a table `(results/arriba/tables/merged_arriba_results.tsv)`, so please check the Arriba results.