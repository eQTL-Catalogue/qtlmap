# divide count_matrix according to sample_metadata file of each qtlgroup
convertDFtoQTLtools <- function(sample_meta_qtlgroup, count_matrix, phenotype_data){
  #Make sure that all required columns are present
  assertthat::assert_that(assertthat::has_name(phenotype_data, "chromosome"))
  assertthat::assert_that(assertthat::has_name(phenotype_data, "phenotype_pos"))
  assertthat::assert_that(assertthat::has_name(phenotype_data, "phenotype_id"))
  assertthat::assert_that(assertthat::has_name(phenotype_data, "group_id"))
  assertthat::assert_that(assertthat::has_name(phenotype_data, "gene_id"))
  assertthat::assert_that(assertthat::has_name(phenotype_data, "strand"))
  
  assertthat::assert_that(assertthat::has_name(sample_meta_qtlgroup, "sample_id"))
  assertthat::assert_that(assertthat::has_name(sample_meta_qtlgroup, "genotype_id"))
  
  #Make sure that duplicate genotype ids are not present in sample metadata
  max_genotype_id_count = max(table(sample_meta_qtlgroup$genotype_id))
  assertthat::assert_that(max_genotype_id_count == 1)
  
  #Make genePos table for QTLTools
  pheno_data = dplyr::arrange(phenotype_data, chromosome, phenotype_pos) %>%
    dplyr::transmute(chromosome, left = phenotype_pos, right = phenotype_pos + 1, phenotype_id, group_id, strand) %>%
    dplyr::rename("#chr" = "chromosome") %>%
    dplyr::mutate(strand = ifelse(strand == 1, "+", "-"))
  
  #Extract phenotype and rename columns according to genotype id
  count_matrix_group <- count_matrix[, colnames(count_matrix) %in% 
     c(sample_meta_qtlgroup$sample_id %>% as.character(), "phenotype_id")]
  match_index <- match(names(count_matrix_group), sample_meta_qtlgroup$sample_id)
  names(count_matrix_group)[!is.na(match_index)] <- as.character(sample_meta_qtlgroup$genotype_id[na.omit(match_index)])
  count_matrix_group[,names(count_matrix_group) != "phenotype_id"] <- 
    round(count_matrix_group[,names(count_matrix_group) != "phenotype_id"], 3) #Round to three digits
  
  pheno_data_filtered <- pheno_data %>% dplyr::filter(phenotype_id %in% count_matrix_group$phenotype_id)

  #Make QTLtools phenotype table
  res = count_matrix_group %>%
    dplyr::select(phenotype_id, dplyr::everything()) %>%
    dplyr::left_join(pheno_data_filtered, ., by = "phenotype_id") %>%
    dplyr::arrange()
  
  message("Exclude phenotypes with zero variance")
  mat = as.matrix(res[,-(1:6)])
  var_vector = apply(mat, 1, var)
  print(length(var_vector))
  res = res[var_vector > 0,]

  return(res)
}

importVariantInformation <- function(path){
  info_col_names = c("chr","pos","snp_id","ref","alt","type","AC","AN", "MAF", "R2")
  into_col_types = c("character", "integer", "character", 
                     "character", "character", "character", "integer", 
                     "integer", "double", "double")
  
  # snp_info = readr::read_delim(path, delim = "\t", col_types = into_col_types, col_names = info_col_names)
  snp_info = utils::read.delim(path, colClasses = into_col_types, col.names = info_col_names, header = FALSE) %>% as.data.frame()
  assertthat::assert_that(nrow(snp_info)>1, msg="variant info is empty")

  snp_info = dplyr::mutate(snp_info, indel_length = pmax(nchar(alt), nchar(ref))) %>%
    dplyr::mutate(is_indel = ifelse(indel_length > 1, TRUE, FALSE)) %>%
    dplyr::mutate(MAF = pmin(AC/AN, 1-(AC/AN)))
  return(snp_info)
}


performQTLGroupPCA <- function(ph_count_matrix_group, output_dir){
  # perfrom PCA for phenotype matrix
  pheno_pca <- stats::prcomp(t(ph_count_matrix_group[-(1:6)]), center=TRUE, scale. = TRUE)
  pheno_pca_x <- t(pheno_pca$x) %>% as.data.frame() 
  
  # change PC column values as into pheno_PC
  pheno_pca_x <- cbind(SampleID = paste0("pheno_", rownames(pheno_pca_x)), pheno_pca_x)
  
  # write phenotype PCA matrix into file
  message(" ## write phenotype PCA matrix to ", file.path(output_dir, "pheno_PCA.tsv"))
  write.table(pheno_pca_x, file = file.path(output_dir, "pheno_PCA.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

# main script which gets inputs and produces outputs
# inputs:
# count_matrix: Tab separated file containing normalized phenotype counts
#             M Columns: 
#                 1st column_name: phenotype_id
#                 starting from 2nd column: sampleIds
#             N Rows: Each row respresents one phenotype
#
# phenotype_data: Tab separated file containing phenotype metadata 
#   *Should have at least the following columns: chromosome, phenotype_pos, strand, phenotype_id
#
# sample_metadata: Tab separated file containing sample metadata
#   *Should have at least the following columns: sample_id, genotype_id, qtl_group
#
# var_info: variant info of the 
#
# cis_distance: Distance to search for variants around the phenotype (look to both sides)
# cis_min_var: Minimum variants needed to be found in cis_distance for further analysis of the phenotype. 
# output_dir: Output directory to write QTLTools input files.

suppressPackageStartupMessages(library("dplyr"))

option_list <- list(
  #TODO look around if there is a package recognizing delimiter in dataset
  optparse::make_option(c("-p", "--phenometa"), type="character", default=NULL,
              help="Phenotype metadata file path of genes used in expression-matrix. Tab separated", metavar = "type"),
  optparse::make_option(c("-s", "--samplemeta"), type="character", default=NULL,
              help="Sample metadata file path of genes used in expression-matrix. Tab separated", metavar = "type"),
  optparse::make_option(c("-e", "--expression_matrix"), type="character", default=NULL,
              help="Expression matrix file path with gene phenotype-id in rownames and sample-is in columnnames", metavar = "type"),
  optparse::make_option(c("-v", "--variant-info"), type="character", default=NULL,
              help="Variant information file path.", metavar = "type"),
  optparse::make_option(c("-o", "--outdir"), type="character", default="./QTLTools_input_files",
              help="Path to the output directory.", metavar = "type"),
  optparse::make_option(c("-c","--cisdistance"), type="integer", default=1000000, 
              help="Cis distance in bases from center of gene. [default \"%default\"]", metavar = "number"),
  optparse::make_option(c("-m", "--mincisvariant"), type="integer", default=5,
              help="Minimum count of cis variants in cis-distance of gene to be taken into account. [default \"%default\"]", metavar = "number"),
  optparse::make_option(c("-a", "--add-covariates"), type="character", default="null",
                        help="Comma-separate list of additonal covariate columns names to be included in the analysis. [default \"%default\"]", metavar = "type")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

#Deubgging
if(FALSE){
  opt = list(
    p = "testdata/GEUVADIS_phenotype_metadata.tsv",
    s = "testdata/GEUVADIS_sample_metadata.tsv",
    e = "testdata/GEUVADIS_cqn.tsv",
    v = "testdata/GEUVADIS_test_ge_LCL.variant_information.txt.gz",
    c = 1000000,
    m = 5,
    o = "./QTLTools_input_files",
    a = "sex")
}

phenotype_meta_path = opt$p
sample_meta_path = opt$s
expression_matrix_path = opt$e
variant_info_path = opt$v
output_dir = opt$o
cis_distance = opt$c
cis_min_var = opt$m
covariate_names = opt$a

message("------ Options parsed ------")
message(paste0("gene_meta_path: ", phenotype_meta_path))
message(paste0("sample_meta_path: ", sample_meta_path))
message(paste0("expression_matrix_path: ", expression_matrix_path))
message(paste0("variant_info_path: ", variant_info_path))
message(paste0("output_dir: ", output_dir))
message(paste0("cis_distance: ", cis_distance))
message(paste0("cis_min_var: ", cis_min_var))
message(paste0("covariate_names: ", covariate_names))

message(" ## Reading gene metadata file")
phenotype_data <- utils::read.delim(phenotype_meta_path, quote = "", header = TRUE, stringsAsFactors = FALSE) %>% base::as.data.frame()

message(" ## Reading sample metadata file")
sample_metadata <- utils::read.delim(sample_meta_path, quote = "", header = TRUE, stringsAsFactors = FALSE) %>% base::as.data.frame()

message(" ## Reading expression matrix")
count_matrix <- utils::read.delim(expression_matrix_path, quote = "", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE) %>% base::as.data.frame() 
print(count_matrix[1:6,1:6])

message(" ## Importing variant info")
var_info = importVariantInformation(variant_info_path)

#Check that all required columns are there
#Check phenotype metadata
assertthat::assert_that(assertthat::has_name(phenotype_data, "chromosome"))
assertthat::assert_that(assertthat::has_name(phenotype_data, "phenotype_pos"))
assertthat::assert_that(assertthat::has_name(phenotype_data, "strand"))
assertthat::assert_that(assertthat::has_name(phenotype_data, "phenotype_id"))
assertthat::assert_that(assertthat::has_name(phenotype_data, "gene_id"))

#Check variant information
assertthat::assert_that(assertthat::has_name(var_info, "chr"))
assertthat::assert_that(assertthat::has_name(var_info, "pos"))

#Check count matrix
assertthat::assert_that(assertthat::has_name(count_matrix, "phenotype_id"))

#Check sample metadata
assertthat::assert_that(assertthat::has_name(sample_metadata, "genotype_id"))
assertthat::assert_that(assertthat::has_name(sample_metadata, "sample_id"))
assertthat::assert_that(assertthat::has_name(sample_metadata, "qtl_group"))

message("Keep samples that are present both in the expression data and sample metdata")
shared_samples = intersect(sample_metadata$sample_id, colnames(count_matrix)[-1])
message(paste0("Number of samples included in the analysis: ", length(shared_samples)))
sample_metadata = dplyr::filter(sample_metadata, sample_id %in% shared_samples)
count_matrix = count_matrix[,c("phenotype_id", shared_samples)]

message(" ## Making gene ranges")
#Make GRanges objects
gene_ranges = GenomicRanges::GRanges(
  seqnames = phenotype_data$chromosome,
  ranges = IRanges::IRanges(
    start = phenotype_data$phenotype_pos - cis_distance,
    end = phenotype_data$phenotype_pos + cis_distance
  ),
  strand = "*"
)
var_ranges = GenomicRanges::GRanges(
  seqnames = var_info$chr,
  ranges = IRanges::IRanges(start = var_info$pos, end = var_info$pos),
  strand = "+"
)
olap_count = GenomicRanges::countOverlaps(gene_ranges, var_ranges, ignore.strand = TRUE)
count_df = dplyr::select(phenotype_data, phenotype_id, chromosome, phenotype_pos) %>%
  dplyr::mutate(snp_count = olap_count) %>%
  dplyr::arrange(snp_count) %>%
  dplyr::filter(snp_count >= cis_min_var)

#Keep phenotypes that contain enough variants nearby
message(" ## Filtering variants for cis_min_var")
phenotype_data <- phenotype_data[phenotype_data$phenotype_id %in% count_df$phenotype_id,]
count_matrix_cis_filter <- count_matrix[count_matrix$phenotype_id %in% count_df$phenotype_id,]

#Split the SE into list based on qtl_group
message(" ## Make sure that the data contains only a single QTL group")
qtl_groups = unique(sample_metadata$qtl_group)
assertthat::assert_that(length(qtl_groups) == 1, msg = "Number of distinct qtl_groups exceeds 1. Please check your sample metadata.")
qtl_group = qtl_groups[1]

message(" ## Prepare data for QTLtools/FastQTL")
qtltools_df = convertDFtoQTLtools(sample_metadata, count_matrix = count_matrix_cis_filter, phenotype_data = phenotype_data)
file_path = file.path(output_dir, paste(qtl_group, "bed", sep = "."))
write.table(qtltools_df, file_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#Extract sample names
sample_names = colnames(qtltools_df)[-(1:6)]
file_path = file.path(output_dir, paste(qtl_group, "sample_names.txt", sep = "."))
write.table(sample_names, file_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#Perform phenotype PCA
transposed_mat <- t(qtltools_df[-(1:6)])
transposed_mat <- transposed_mat[, apply(transposed_mat, 2, var) != 0]
pheno_pca <- stats::prcomp(transposed_mat, center=TRUE, scale. = TRUE)
pheno_pca_x <- t(pheno_pca$x) %>% as.data.frame()

# change PC column values as into pheno_PC
pheno_pca_x <- cbind(SampleID = paste0("pheno_", rownames(pheno_pca_x)), pheno_pca_x)
utils::write.table(pheno_pca_x, file = file.path(output_dir, "phenoPCA.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#Prepare covariates
if(!(covariate_names == "null")){
  covariate_values = strsplit(covariate_names, ",")[[1]]
  rownames(sample_metadata) = sample_metadata$genotype_id
  cov_meta = sample_metadata[sample_names,covariate_values, drop = FALSE]
  
  #Convert character columns to factors levels
  for (i in seq_along(colnames(cov_meta))){
    if(class(cov_meta[,i]) == "character"){
      cov_meta[,i] = as.integer(as.factor(cov_meta[,i]))
    }
  }
  cov_df = t(cov_meta)
  cov_df <- cbind(SampleID = paste0("cov_", rownames(cov_df)), cov_df)
  utils::write.table(cov_df, file = file.path(output_dir, "additional_covariates.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
} else {
  #Write an empty file with just the header
  utils::write.table(pheno_pca_x[0,], file = file.path(output_dir, "additional_covariates.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}


