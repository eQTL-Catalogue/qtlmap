# nf-core/qtlmap
**This pipeline mapping "Quantitive Trait Loci"-s**

[![Build Status](https://travis-ci.org/nf-core/qtlmap.svg?branch=master)](https://travis-ci.org/nf-core/qtlmap)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/qtlmap.svg)](https://hub.docker.com/r/nfcore/qtlmap)
![Singularity Container available](
https://img.shields.io/badge/singularity-available-7E4C74.svg)

### Introduction
The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible.


### Documentation
The nf-core/qtlmap pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. Pipeline configuration
    * [Local installation](docs/configuration/local.md)
    * [Adding your own system](docs/configuration/adding_your_own.md)
    * [Reference genomes](docs/configuration/reference_genomes.md)  
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

<!-- TODO nf-core: Add a brief overview of what the pipeline does and how it works -->

### Pipeline Description
Mapping QTLs is a process of finding statistically significant associations between phenotypes and genetic variants located nearby (within a specific window around phenotype, a.k.a cis window)
This pipeline is designed to perform QTL mapping. It is intended to add this pipeline to the nf-core framework in the future.
High level representation of the pipeline is shown below:

![High_level_schema](docs/images/QTLMap_pipeline_high_level_repr.png)


### Credits
kerimoff/qtlmap was originally written by Nurlan Kerimov under supervision of [Kaur Alasoo](https://github.com/kauralasoo)
