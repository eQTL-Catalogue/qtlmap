#!/usr/bin/env Rscript

message(" ## Loading libraries: devtools, dplyr, SummarizedExperiment, optparse")

# TODO: write dependecy description like installation of limma and lumi
suppressPackageStartupMessages(library("devtools"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("SummarizedExperiment"))
suppressPackageStartupMessages(library("optparse"))

#Parse command-line options
option_list <- list(
  #TODO look around if there is a package recognizing delimiter in dataset
  make_option(c("-g", "--genemeta"), type="character", default=NULL,
              help="Gene metadata file path of genes used in expression-matrix. Tab separated", metavar = "type"),
  make_option(c("-s", "--samplemeta"), type="character", default=NULL,
              help="Sample metadata file path of genes used in expression-matrix. Tab separated", metavar = "type"),
  make_option(c("-e", "--expression_matrix"), type="character", default=NULL,
              help="Expression matrix file path with gene phenotype-id in rownames and sample-is in columnnames", metavar = "type"),
  make_option(c("-v", "--variant-info"), type="character", default=NULL,
              help="Variant information file path.", metavar = "type"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="Path to the output directory.", metavar = "type"),
  #TODO change default of eQTL utils to github path
  make_option(c("--qtlutils"), type="character", default="../../eQTLUtils",
              help="eQTLUtils path to be loaded by devtools. [default \"%default\"]", metavar = "type"),
  make_option(c("-c","--cisdistance"), type="integer", default=1000000, 
              help="Cis distance in bases from center of gene. [default \"%default\"]", metavar = "number"),
  make_option(c("-m", "--mincisvariant"), type="integer", default=5,
              help="Minimum count of cis variants in cis-distance of gene to be taken into account. [default \"%default\"]", metavar = "number")
)

message(" ## Parsing options")
opt <- parse_args(OptionParser(option_list=option_list))

#Debugging
if (FALSE) {
  opt = list(c=1000000, m=6)
  opt$g="../metadata/gene_metadata/HumanHT-12_V4_gene_metadata.txt.gz"
  opt$s="../metadata/cleaned/Fairfax_2014.tsv"
  opt$e="../results/expression_matrices/HumanHT-12_V4/Fairfax_2014.tsv.gz"
  opt$v="../../temp/Fairfax_2014_GRCh38.variant_information.txt.gz"
  opt$qtlutils="../../eQTLUtils/"
  opt$o="../processed/Fairfax_2014/qtltools/input/array/" #-c 1000001 -m 6
}

gene_meta_path = opt$g
sample_meta_path = opt$s
expression_matrix_path = opt$e
variant_info_path = opt$v
output_dir = opt$o
cis_dist = opt$c
cis_min_var = opt$m
eqtl_utils_path = opt$qtlutils

load_all(eqtl_utils_path)

message("------ Options parsed ------")
message(paste0("eqtl_utils_path: ", eqtl_utils_path))
message(paste0("gene_meta_path: ", gene_meta_path))
message(paste0("sample_meta_path: ", sample_meta_path))
message(paste0("expression_matrix_path: ", expression_matrix_path))
message(paste0("variant_info_path: ", variant_info_path))
message(paste0("output_dir: ", output_dir))
message(paste0("cis_dist: ", cis_dist))
message(paste0("cis_min_var: ", cis_min_var))

#keep chromosomes
valid_chromosomes = c("1","10","11","12","13","14","15","16","17","18","19",
                      "2","20","21","22","3","4","5","6","7","8","9")

#Genrate Summarized Experiment object of the study
message(" ## Reading gene metadata file")
gene_meta = readr::read_tsv(gene_meta_path)

message(" ## Reading sample metadata file")
sample_metadata = readr::read_tsv(sample_meta_path)

message(" ## Reading expression matrix")
expression_matrix = read.table(expression_matrix_path, sep = "\t")

message(" ## Generating Summarized Experiment object of the study")
fairfax_2014_se = eQTLUtils::makeSummarizedExperiment(expression_matrix, gene_meta, sample_metadata, assay_name = "exprs")

#Filter SE to keep correct chromosomes and QC-passed samples
message(" ## Filterin out invalid chromosomes, RNA_QC failed samples and genotype_QC failed samples ")
#TODO add rna_qc and genotype_qc as optional parameters to optparse
fairfax_se_filtered = suppressWarnings(eQTLUtils::filterSummarizedExperiment(fairfax_2014_se, valid_chromosomes = valid_chromosomes, filter_rna_qc = TRUE, filter_genotype_qc = TRUE))

#Normalize and regress out batch effects
message(" ## Normalizing and regressing out the batch effects")
fairfax_norm = suppressWarnings(eQTLUtils::array_normaliseSE(fairfax_se_filtered, norm_method = "quantile", assay_name = "exprs", log_transform = TRUE, adjust_batch = TRUE, filter_quality = TRUE))

#Count the number of variants proximal to each gene and 
message(" ## Importing variant information")
var_info = eQTLUtils::importVariantInformation(variant_info_path)

#Remove genes without minimum variants in the region
message(" ## Filterin out genes which do not have enough (cis_min_var) variants in the region (cis_dist)")
fairfax_norm_filtered = eQTLUtils::checkCisVariants(fairfax_norm, var_info, cis_distance = cis_dist, min_cis_variant = cis_min_var)

#Export expression data to disk
message(" ## Exporting expression data to \'", output_dir, "\' to feed QTLTools with input")
eQTLUtils::studySEtoQTLTools(fairfax_norm_filtered, assay_name = "norm_exprs", output_dir)

message(" ## Input for QTLtools are generated in \'", output_dir, "\'")