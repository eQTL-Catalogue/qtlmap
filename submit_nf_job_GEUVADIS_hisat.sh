#!/bin/bash

#SBATCH --time=36:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --job-name="HISAT_qtlmap"
#SBATCH --partition=amd

# Load needed system tools (Java 8 is required, one of singularity or anaconda - python 2.7 is needed,
# depending on the method for dependancy management). The exact names of tool modules might depend on HPC.
module load java-1.8.0_40
module load nextflow
module load singularity/3.5.3
module load squashfs/4.4

nextflow run main.nf -profile tartu_hpc,eqtl_catalogue -resume\
  --studyFile /gpfs/space/home/kerimov/qcnorm_fast/results_leafcutter_HISAT_GEUVADIS/GEUVADIS/GEUVADIS_qtlmap_inputs.tsv\
  --outdir /gpfs/space/home/kerimov/qtlmap_lc_hisat/GEUVADIS_HISAT_qtlmap_results/\
  --vcf_has_R2_field FALSE\
  --vcf_genotype_field DS\
  --covariates sex\
  --susie_skip_full TRUE


