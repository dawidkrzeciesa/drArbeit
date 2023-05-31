# Important

For this project, the [dna-seq-mtb](https://github.com/snakemake-workflows/dna-seq-mtb) pipeline (v1.3.3) was used. 

This pipeline utilizes the [dna-seq-varlociraptor](https://github.com/snakemake-workflows/dna-seq-varlociraptor) pipeline (v3.19.6) as a module, which is located in the workflow/resources/dna-seq-varlociraptor.v3.19.6 directory.

### The following are the significant modifications made to the original pipelines:

**dna-seq-varlociraptor**
- Updated varlociraptor version to 5.5.0

**dna-seq-mtb**
- Added a 3% contamination of tumor in normal tissue to the default variant calling scenario (workflow/resources/config/scenario.yaml).
- Added a filter to include events with a variant allele frequency of 25% or higher (workflow/resources/config/default.yaml). 



# Snakemake workflow: dna-seq-mtb

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/snakemake-workflows/dna-seq-mtb/workflows/Tests/badge.svg?branch=main)](https://github.com/snakemake-workflows/dna-seq-mtb/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for calling somatic variants in tumor/normal samples.
This workflow is under heavy development. 
A tumor-only mode will be added soon, as well as options for handling FFPE samples.

## Background

This workflow uses the [dna-seq-varlociraptor workflow](https://github.com/snakemake-workflows/dna-seq-varlociraptor) as a [module](https://github.com/snakemake-workflows/dna-seq-mtb/blob/main/workflow/Snakefile#L4), and configures it with reasonable defaults for running analysis in molecular tumor boards.


## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=snakemake-workflows%2Fdna-seq-mtb).

If you use this workflow in a paper, don't forget to give credits to the authors by citing https://github.com/snakemake-workflows/dna-seq-mtb.
