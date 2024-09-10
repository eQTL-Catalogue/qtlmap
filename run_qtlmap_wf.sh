#!/bin/bash

#SBATCH --time=03:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5G
#SBATCH --job-name="qtl"
#SBATCH --partition=amd


module load any/jdk/1.8.0_265
module load nextflow
module load any/singularity/3.7.3
module load squashfs/4.4



 NXF_VER=22.04.3 nextflow run main.nf -profile tartu_hpc -resume\
    --studyFile testdata/multi_test.tsv\
     --vcf_has_R2_field FALSE\
     --run_permutation FALSE\
     --rsid_map_file /gpfs/helios/home/a82371/qtlmap_update/qtlmap/data/rsid_map_file.tsv\
     --n_batches 25