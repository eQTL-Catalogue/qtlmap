# kerimoff/qtlmap: Input files explanation

This pipeline requires 4 mandatory input files:

* **Phenotype Count Matrix (.tsv)**: Is the tab separated tabular data file which contains normalized and quality controlled phenotype counts. Should contain at least the following columns: **_phenotype_id_**

* **Sample Metadata (.tsv)**: Is the tab separated tabular data file which contains metadata of the samples represented in phenotype count matrix. Should contain at least the following columns: **_sample_id, genotype_id, qtl_group_**

* **Phenotype Metadata (.tsv)**: Is the tab separated tabular data file which contains metadata of the phenotypes represented in phenotype count matrix. Should contain at least the following columns: **_phenotype_id, chromosome, phenotype_pos, strand_**

* **Genotype data (VCF or VCF.gz)**: Is the VCF file which contains genotypic data of the all samples represented in phenotype count matrix. Data column names of the VCF file should correspond to the _**genotype_id**_ column values of sample_metadata.


## Relationships between input files are shown below.

![input_relationships](https://github.com/kerimoff/qtlmap/blob/master/docs/images/input_relations.svg)

