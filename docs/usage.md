# kerimoff/qtlmap: Usage

## Table of contents

* [Introduction](#general-nextflow-info)
* [Running the pipeline](#running-the-pipeline)
* [Main arguments](#main-arguments)
    * [`--expression_matrix`](#--expression_matrix)
    * [`--phenotype_metadata`](#--phenotype_metadata)
    * [`--sample_metadata`](#--sample_metadata)
    * [`--genotype_vcf`](#--genotype_vcf)
    * [`--cis_window`](#--cis_window)
    * [`--mincisvariant`](#--mincisvariant)
    * [`--n_batches`](#--n_batches)
    * [`--is_imputed`](#--is_imputed)
* [Using profiles](#-profile)
    * [`-profile`](#-profile-single-dash)
       * [`awsbatch`](#awsbatch)
       * [`conda`](#conda)
       * [`docker`](#docker)
       * [`singularity`](#singularity)
       * [`test`](#test)
* [Job resources](#job-resources)
* [Automatic resubmission](#automatic-resubmission)
* [Custom resource requests](#custom-resource-requests)
* [AWS batch specific parameters](#aws-batch-specific-parameters)
    * [`--awsqueue`](#--awsqueue)
    * [`--awsregion`](#--awsregion)
* [Other command line parameters](#other-command-line-parameters)
    * [`--outdir`](#--outdir)
    * [`--email`](#--email)
    * [`-name`](#-name)
    * [`-resume`](#-resume)
    * [`-c`](#-c)
    * [`--custom_config_version`](#--custom_config_version)
    * [`--max_memory`](#--max_memory)
    * [`--max_time`](#--max_time)
    * [`--max_cpus`](#--max_cpus)
    * [`--plaintext_email`](#--plaintext_email)


## General Nextflow info
Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline
The typical command for running the pipeline is as follows:
```bash
nextflow run kerimoff/qtlmap\
     --expression_matrix testdata/GEUVADIS_cqn.tsv\
     --phenotype_metadata testdata/GEUVADIS_phenotype_metadata.tsv\
     --sample_metadata testdata/GEUVADIS_sample_metadata.tsv\
     --genotype_vcf testdata/GEUVADIS_genotypes.vcf.gz
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

# Main arguments

## Mandatory Arguments
### `--expression_matrix`
Use this to specify the location of your phenotype matrix file (.tsv). Should contain at least the following columns: **_phenotype_id_**

For example:

```bash
nextflow run main.nf --expression_matrix testdata/GEUVADIS_cqn.tsv
```

```groovy
params {
    expression_matrix = "$baseDir/testdata/GEUVADIS_cqn.tsv"
}
```

### `--phenotype_metadata`
Use this to specify the location of your phenotype metadata file (.tsv). Should contain at least the following columns: **_phenotype_id, chromosome, phenotype_pos, strand_**

For example:

```bash
nextflow run main.nf --phenotype_metadata testdata/GEUVADIS_phenotype_metadata.tsv
```

```groovy
params {
    phenotype_metadata = "$baseDir/testdata/GEUVADIS_phenotype_metadata.tsv"
}
```

### `--sample_metadata`
Use this to specify the location of your sample metadata file (.tsv). Should contain at least the following columns: **_sample_id, genotype_id, qtl_group_**

For example:

```bash
nextflow run main.nf --sample_metadata testdata/GEUVADIS_sample_metadata.tsv
```

```groovy
params {
    sample_metadata = "$baseDir/testdata/GEUVADIS_sample_metadata.tsv"
}
```

### `--genotype_vcf`
Use this to specify the location of your genotype file (.vcf of .vcf.gz). For example:
Data column names of the VCF file should correspond to the _**genotype_id**_ column values of sample_metadata. See [Input files explanations](docs/inputs_expl.md) for more details

```bash
nextflow run main.nf --genotype_vcf testdata/GEUVADIS_genotypes.vcf.gz
```

```groovy
params {
    genotype_vcf = "$baseDir/testdata/GEUVADIS_genotypes.vcf.gz"
}
```
## Optional Arguments
### `--cis_window`
Use this to specify the cis-window length in bases. Deafult value is _**1,000,000 bases (1Mb)**_

```bash
nextflow run main.nf [mandatory arguments here] --cis_window 1500000
```

```groovy
params {
    cis_window = 1500000
}
```

### `--mincisvariant`
Use this to specify the threshold minimum variants in specified cis window. If for specific phenotype, there are less variants found in cis window than this threshold, the phenotype is filtered out and not processed further. Default value is _**5**_

```bash
nextflow run main.nf [mandatory arguments here] --mincisvariant 10
```

```groovy
params {
    mincisvariant = 10
}
```

### `--n_batches`
Use this to specify the number of batches used in QTL mapping run. QTLTools will split the genome into this number of chunks and perform the run in a parallel manner. The default value is _**400**_

```bash
nextflow run main.nf [mandatory arguments here] --n_batches 200
```

```groovy
params {
    n_batches = 200
}
```

### `--is_imputed`
Use this to specify if the provided genotype (VCF) file is imputed or not. The imputed genotype files should have the following columns _**CHROM, POS, ID, REF, ALT, TYPE, AC, AN, MAF, R2**_. The default value is _**TRUE**_

```bash
nextflow run main.nf [mandatory arguments here] --is_imputed FALSE
```

```groovy
params {
    is_imputed = false
}
```

## Using profiles
### `-profile`
Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile docker` - the order of arguments is important!

If `-profile` is not specified at all the pipeline will be run locally and expects all software to be installed and available on the `PATH`.

* `awsbatch`
    * A generic configuration profile to be used with AWS Batch.
* `conda`
    * A generic configuration profile to be used with [conda](https://conda.io/docs/)
    * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `docker`
    * A generic configuration profile to be used with [Docker](http://docker.com/)
    * Pulls software from dockerhub: [`nfcore/qtlmap`](http://hub.docker.com/r/nfcore/qtlmap/)
* `singularity`
    * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
    * Pulls software from singularity-hub
* `test`
    * A profile with a complete configuration for automated testing
    * Includes links to test data so needs no other parameters


## Job resources
### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests
Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [`Slack`](https://nf-core-invite.herokuapp.com/).

## AWS Batch specific parameters
Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use the `-awsbatch` profile and then specify all of the following parameters.
### `--awsqueue`
The JobQueue that you intend to use on AWS Batch.
### `--awsregion`
The AWS region to run your job in. Default is set to `eu-west-1` but can be adjusted to your needs.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

## Other command line parameters

<!-- TODO nf-core: Describe any other command line flags here -->

### `--outdir`
The output directory where the results will be saved.

### `--email`
Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to speicfy this on the command line for every run.

### `-name`
Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`
Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--custom_config_version`
Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default is set to `master`.

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

### `--max_memory`
Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`
Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`
Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`
Set to receive plain-text e-mails instead of HTML formatted.
