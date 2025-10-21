# eQTL-Catalogue/qtlmap
**Portable eQTL analysis and statistical fine mapping workflow used by the eQTL Catalogue**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)

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

# Test usage

1. download repository
```shell
git clone git@github.com:MiqG/qtlmap.git
```

2. download dbSNP datbase for finemapping
```shell
wget "https://zenodo.org/records/15170247/files/dbSNP_b151_GRCh38p7_parquet.tar.gz?download=1" -O dbSNP_b151_GRCh38p7_parquet.tar.gz
tar -zxvf dbSNP_b151_GRCh38p7_parquet.tar.gz
rm dbSNP_b151_GRCh38p7_parquet.tar.gz
```

3. run small test
```shell
nextflow run main.nf \
    -profile singularity \
    -resume \
    --studyFile testdata/multi_test.tsv \
    --rsid_map_file rsid_map/rsid_map_file.tsv \
    --max_memory 10.GB \
    --max_time 2.h \
    --max_cpus 1 \
    -c <(echo "process { withName: make_pca_covariates { cpus = 1 } }") \
    --sumstat_sort_cores 1 \
    --sumstat_sort_mem "4G" \
    --cis_window 1000000 \
    --mincisvariant 5 \
    --n_geno_pcs 3 \
    --n_pheno_pcs 3 \
    --run_permutation TRUE \
    --n_permutations 1000 \
    --vcf_has_R2_field FALSE \
    --n_batches 25 \
    --run_susie TRUE \
    --write_full_susie FALSE \
    --outdir testdata/test_results/
```

4. run with eQTL Catalogue settings
```shell
nextflow run main.nf \
    -profile singularity \
    -resume \
    --studyFile testdata/multi_test.tsv \
    --rsid_map_file rsid_map/rsid_map_file.tsv \
    --max_memory 10.GB \
    --max_time 2.h \
    --max_cpus 1 \
    -c <(echo "process { withName: make_pca_covariates { cpus = 1 } }") \
    --cis_window 1000000 \
    --n_geno_pcs 6 \
    --n_pheno_pcs 6 \
    --run_permutation TRUE \
    --n_batches 400 \
    --mincisvariant 5 \
    --n_permutations 1000 \
    --run_susie TRUE \
    --write_full_susie FALSE \
    --outdir testdata/test_results/
```

# Contributors
* Nurlan Kerimov
* Kaur Alasoo
* Masahiro Kanai
* Ralf Tambets
* Krista Freimann
