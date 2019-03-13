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
              help="Minimum count of cis variants in cis-distance of gene to be taken into account. [default \"%default\"]", metavar = "number"),
  make_option(c("--quantification"), type="character", default="featureCounts",
              help="Quantification method used. Currently suppprted: featureCounts, array. [default \"%default\"]", metavar = "string"),
  make_option(c("--qtl_group_extra"), type="character", default="NULL",
              help="Name of the metadata column used in conjunction with qtl_group to split the dataset into independent QTLTools input files.", metavar = "string")
)

message(" ## Parsing options")
opt <- parse_args(OptionParser(option_list=option_list))

#Debugging
if (FALSE) {
  opt = list(c=1000000, m=7)
  opt$g="../metadata/gene_metadata/HumanHT-12_V4_gene_metadata.txt.gz"
  opt$s="../metadata/cleaned/Fairfax_2014.tsv"
  opt$e="../results/expression_matrices/HumanHT-12_V4/Fairfax_2014.tsv.gz"
  opt$v="../../temp/Fairfax_2014_GRCh38.variant_information.txt.gz"
  opt$qtlutils="../../eQTLUtils/"
  opt$o="../processed/Fairfax_2014/qtltools/input/array/" #-c 1000001 -m 6
  opt$quantification = "array"
}

if (FALSE) {
  opt = list(c=1000000, m=7)
  opt$g= "metadata/gene_metadata/featureCounts_Ensembl_92_gene_metadata.txt.gz"
  opt$s="metadata/cleaned/GENCORD.tsv"
  opt$e="results/expression_matrices/featureCounts/GENCORD.tsv.gz"
  opt$v="results/var_info/GENCORD_GRCh38.variant_information.txt.gz"
  opt$qtlutils="../eQTLUtils/"
  opt$o="processed/GENCORD/qtltools/input/featureCounts/"
  opt$quantification = "featureCounts"
}

if (FALSE) {
  opt = list(c=1000000, m=7)
  opt$g= "metadata/gene_metadata/HumanHT-12_V4_gene_metadata.txt.gz"
  opt$s="metadata/cleaned/Fairfax_2014.tsv"
  opt$e="results/expression_matrices/HumanHT-12_V4/Fairfax_2014.tsv.gz"
  opt$v="results/var_info/Fairfax_2014_GRCh38.variant_information.txt.gz"
  opt$qtlutils="../eQTLUtils/"
  opt$o="processed/Fairfax_2014/qtltools/input/array/"
  opt$quantification = "array"
  opt$qtl_group_extra = "sex"
}

gene_meta_path = opt$g
sample_meta_path = opt$s
expression_matrix_path = opt$e
variant_info_path = opt$v
output_dir = opt$o
cis_dist = opt$c
cis_min_var = opt$m
eqtl_utils_path = opt$qtlutils
quant_method = opt$quantification
if(opt$qtl_group_extra == "NULL"){
  extra_qtl_group = NULL
} else{
  extra_qtl_group = opt$qtl_group_extra
}

#Set input and output assay names based on quantification method
if(quant_method == "featureCounts"){
  input_assay_name = "counts"
  normalised_assay_name = "cqn"
} else if (quant_method == "array"){
  input_assay_name = "exprs"
  normalised_assay_name = "norm_exprs"
} else if (quant_method == "LeafCutter"){
  input_assay_name = "counts"
  normalised_assay_name = "qnorm"
}

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

#Genrate Summarized Experiment object of the study
message(" ## Reading gene metadata file")
gene_meta = read.table(gene_meta_path, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

message(" ## Reading sample metadata file")
sample_metadata = readr::read_tsv(sample_meta_path)

message(" ## Reading expression matrix")
expression_matrix = read.table(expression_matrix_path, sep = "\t", check.names = FALSE)

message(" ## Generating Summarized Experiment object of the study")
raw_se = eQTLUtils::makeSummarizedExperiment(expression_matrix, gene_meta, sample_metadata, assay_name = input_assay_name)

message(" ## Prepare SummarizedExperiment object for QTLtools")
norm_se = eQTLUtils::qtltoolsPrepareSE(raw_se, quant_method = quant_method)

#Count the number of variants proximal to each gene and 
message(" ## Importing variant information")
var_info = eQTLUtils::importVariantInformation(variant_info_path)

#Remove genes without minimum variants in the region
message(" ## Filterin out genes which do not have enough (cis_min_var) variants in the region (cis_dist)")
norm_filtered_se = eQTLUtils::checkCisVariants(norm_se, var_info, cis_distance = cis_dist, min_cis_variant = cis_min_var)

#Export expression data to disk
message(" ## Exporting expression data to \'", output_dir, "\' to feed QTLTools with input")
eQTLUtils::studySEtoQTLTools(norm_filtered_se, assay_name = normalised_assay_name, output_dir, extra_qtl_group = extra_qtl_group)

message(" ## Input for QTLtools are generated in \'", output_dir, "\'")