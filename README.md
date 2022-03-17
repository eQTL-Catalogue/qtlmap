# eQTL-Catalogue/qtlmap
**Portable eQTL analysis and statistical fine mapping workflow used by the eQTL Catalogue**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/kerimoff/qtlmap.svg)](https://hub.docker.com/r/kerimoff/qtlmap)
[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/2842)

### Introduction

**eQTL-Catalogue/qtlmap** is a bioinformatics analysis pipeline used for QTL Analysis.

The workflow takes phenotype count matrix (normalized and quality controlled) and genotype data as input, and finds associations between them with the help of sample metadata and phenotype metadata files (See [Input formats and preparation](docs/inputs_expl.md) for required input file details). To map QTLs, pipeline uses [QTLTools's](https://qtltools.github.io/qtltools/) PCA and RUN methods. For manipulation of files [BcfTools](https://samtools.github.io/bcftools/bcftools.html), [Tabix](http://www.htslib.org/doc/tabix.html) and custom [Rscript](https://www.rdocumentation.org/packages/utils/versions/3.5.3/topics/Rscript) scripts are used.

The pipeline is built using [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible.


### Documentation
The eQTL-Catalogue/qtlmap pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. Pipeline configuration
    * [Local installation](docs/configuration/local.md)
    * [Adding your own system](docs/configuration/adding_your_own.md)
3. [Input formats and preparation](docs/inputs_expl.md)
4. [Running the pipeline](docs/usage.md)
5. [Troubleshooting](docs/troubleshooting.md)

<!-- TODO nf-core: Add a brief overview of what the pipeline does and how it works -->

### Pipeline Description
Mapping QTLs is a process of finding statistically significant associations between phenotypes and genetic variants located nearby (within a specific window around phenotype, a.k.a cis window)
This pipeline is designed to perform QTL mapping. It is intended to add this pipeline to the nf-core framework in the future.
High level representation of the pipeline is shown below:

![High_level_schema](docs/images/metromap.png)

### Results
The output directory of the workflow contains the following subdirectories:

1. PCA - genotype and gene expression PCA values used as covariates for QTL analysis.
2. sumstats - QTL summary statistics from nominal and permutation passes.
3. susie - SuSiE fine mapping credible sets.
4. susie_full - full set of susie results for all tested variants (very large files).
5. susie_merged - susie credible sets merged with summary statistics from univariate QTL analysis.

Column names of the output files are explained [here](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tabix/Columns.md).


# Contributors
* Nurlan Kerimov
* Kaur Alasoo
* Masahiro Kanai
* Ralf Tambets
