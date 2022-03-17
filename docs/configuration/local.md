# eQTL-Catalogue/qtlmap: Local Configuration

If running the pipeline in a local environment, we highly recommend using either Docker or Singularity.

## Docker
Docker is a great way to run `eQTL-Catalogue/qtlmap`, as it manages all software installations and allows the pipeline to be run in an identical software environment across a range of systems.

Nextflow has [excellent integration](https://www.nextflow.io/docs/latest/docker.html) with Docker, and beyond installing the two tools, not much else is required. The `eQTL-Catalogue/qtlmap` profile comes with a configuration profile for docker, making it very easy to use. 

First, install docker on your system: [Docker Installation Instructions](https://docs.docker.com/engine/installation/)

Then, simply run the analysis pipeline:
```bash
nextflow run eQTL-Catalogue/qtlmap -profile docker --studyFile ...
```

Nextflow will recognise `eQTL-Catalogue/qtlmap` and download the pipeline from GitHub. The `-profile docker` configuration lists the [eQTL-Catalogue/qtlmap](https://quay.io/repository/eqtlcatalogue/qtlmap?tag=v20.05.1) image that we have created and is hosted at quay.io, and this is downloaded.

## Singularity image
Many HPC environments are not able to run Docker due to security issues. [Singularity](http://singularity.lbl.gov/) is a tool designed to run on such HPC systems which is very similar to Docker. Even better, it can use create images directly from dockerhub or any other image repository.

```bash
git clone https://github.com/eQTL-Catalogue/qtlmap.git
cd qtlmap
nextflow run main.nf -profile singularity,test
```
