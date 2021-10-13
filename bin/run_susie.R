suppressPackageStartupMessages(library("devtools"))
suppressPackageStartupMessages(library("susieR"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("Rsamtools"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("dplyr"))

option_list <- list(
  #TODO look around if there is a package recognizing delimiter in dataset
  optparse::make_option(c("--phenotype_meta"), type="character", default=NULL,
                        help="Phenotype metadata file path of genes used in expression-matrix. Tab separated", metavar = "type"),
  optparse::make_option(c("--sample_meta"), type="character", default=NULL,
                        help="Sample metadata file path of genes used in expression-matrix. Tab separated", metavar = "type"),
  optparse::make_option(c("--expression_matrix"), type="character", default=NULL,
                        help="Expression matrix file path with gene phenotype-id in rownames and sample-is in columnnames", metavar = "type"),
  optparse::make_option(c("--phenotype_list"), type="character", default=NULL,
                        help="Path to the phenotype list file.", metavar = "type"),
  optparse::make_option(c("--genotype_matrix"), type="character", default=NULL,
                        help="Genotype dosage matrix extracted from VCF.", metavar = "type"),
  optparse::make_option(c("--covariates"), type="character", default=NULL,
                        help="Path to covariates file in QTLtools format.", metavar = "type"),
  optparse::make_option(c("--out_prefix"), type="character", default="./finemapping_output",
                        help="Prefix of the output files.", metavar = "type"),
  optparse::make_option(c("--qtl_group"), type="character", default=NULL,
                        help="Value of the current qtl_group.", metavar = "type"),
  optparse::make_option(c("--cisdistance"), type="integer", default=1000000, 
                        help="Cis distance in bases from center of gene. [default \"%default\"]", metavar = "number"),
  optparse::make_option(c("--chunk"), type="character", default="1 1", 
                        help="Perform analysis in chunks. Eg value 5 10 would indicate that phenotypes are split into 10 chunks and the 5th one of those will be processed. [default \"%default\"]", metavar = "type"),
  optparse::make_option(c("--eqtlutils"), type="character", default=NULL,
              help="Optional path to the eQTLUtils R package location. If not specified then eQTLUtils is assumed to be installed in the container. [default \"%default\"]", metavar = "type"),
  optparse::make_option(c("--permuted"), type="character", default="true",
                        help="If 'false', lead variants were extracted from nominal p-value files. [default \"%default\"]", metavar = "type")
  
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

#Debugging
if(FALSE){
  opt = list(phenotype_list = "testdata/Kasela_2017.T-cell_CD4_microarray.permuted.txt.gz",
             cisdistance = 200000,
             genotype_matrix = "testdata/GEUVADIS_genotypes.dose.tsv.gz",
             covariates = "testdata/GEUVADIS_test_ge.covariates.txt",
             expression_matrix = "testdata/GEUVADIS_cqn.tsv",
             sample_meta = "testdata/GEUVADIS_sample_metadata.tsv",
             phenotype_meta = "testdata/GEUVADIS_phenotype_metadata.tsv",
             chunk = "1 1",
             out_prefix = "./finemapping_output",
             eqtlutils = "../eQTLUtils/",
             qtl_group = "LCL",
             permuted = "true"
  )
}

#Load eQTLUtils
if(opt$eqtlutils == "null"){
  opt$eqtlutils = NULL
}
if (!is.null(opt$eqtlutils)){
  devtools::load_all(opt$eqtlutils)
}

#Print all options
print(opt)

#Define helper functions
importQtlmapCovariates <- function(covariates_path){
  pc_matrix = read.table(covariates_path, check.names = F, header = T, stringsAsFactors = F)
  pc_transpose = t(pc_matrix[,-1])
  colnames(pc_transpose) = pc_matrix$SampleID
  pc_df = dplyr::mutate(as.data.frame(pc_transpose), genotype_id = rownames(pc_transpose)) %>%
    dplyr::as_tibble() %>% 
    dplyr::select(genotype_id, dplyr::everything())
  
  #Make PCA matrix
  pc_matrix = as.matrix(dplyr::select(pc_df,-genotype_id))
  rownames(pc_matrix) = pc_df$genotype_id
  return(pc_matrix)
}

importQtlmapPermutedPvalues <- function(perm_path){
  tbl = read.table(perm_path, check.names = F, header = T, stringsAsFactors = F) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(p_fdr = p.adjust(p_beta, method = "fdr")) %>%
    dplyr::mutate(group_id = molecular_trait_object_id)
  return(tbl)
}

splitIntoBatches <- function(n, batch_size){
  n_batches = ceiling(n/batch_size)
  batch_ids = rep(seq(1:n_batches), each = batch_size)[1:n]
  return(batch_ids)
}

splitIntoChunks <- function(chunk_number, n_chunks, n_total){
  chunk_size = max(1,floor(n_total/(n_chunks)))
  batches = splitIntoBatches(n_total,chunk_size)
  batches[batches > n_chunks] = n_chunks
  selected_batch = batches == chunk_number
  return(selected_batch)
}

finemapPhenotype <- function(phenotype_id, se, genotype_file, covariates, cis_distance){
  message("Processing phenotype: ", phenotype_id)
  
  #Extract phenotype from SE
  gene_vector = eQTLUtils::extractPhentypeFromSE(phenotype_id, se, "counts") %>%
    dplyr::mutate(phenotype_value_std = qnorm((rank(phenotype_value, na.last = "keep") - 0.5) / sum(!is.na(phenotype_value))))
  selected_phenotype = phenotype_id
  gene_meta = dplyr::filter(SummarizedExperiment::rowData(se) %>% as.data.frame(), phenotype_id == selected_phenotype)

  #Rearrange samples in the covariates matrix
  covariates_matrix = cbind(covariates[gene_vector$genotype_id,], 1)
  
  #Import genotype matrix
  genotype_matrix = eQTLUtils::extractGenotypeMatrixFromDosage(
    chr = gene_meta$chromosome, 
    start = gene_meta$phenotype_pos - cis_distance, 
    end = gene_meta$phenotype_pos + cis_distance, 
    dosage_file = genotype_file)

  #Residualise gene expression and genotype matrix
  hat = diag(nrow(covariates_matrix)) - covariates_matrix %*% solve(crossprod(covariates_matrix)) %*% t(covariates_matrix)
  expression_vector = hat %*% gene_vector$phenotype_value_std
  names(expression_vector) = gene_vector$genotype_id
  
  gt_matrix = genotype_matrix[,names(expression_vector)]
  
  #Exclude variants with no alternative alleles
  gt_matrix = gt_matrix[rowSums(round(gt_matrix,0), na.rm = TRUE) != 0,]
  
  #Replace missing values with row means
  gt_matrix = t(gt_matrix) %>% zoo::na.aggregate() %>% t()

  #Standardise genotypes
  gt_std = t(gt_matrix - apply(gt_matrix, 1, mean))
  gt_hat = hat %*% gt_std
  
  # Fit finemapping model
  fitted <- susieR::susie(gt_hat, expression_vector,
                          L = 10,
                          estimate_residual_variance = TRUE, 
                          estimate_prior_variance = TRUE,
                          scaled_prior_variance = 0.1,
                          verbose = TRUE,
                          compute_univariate_zscore = TRUE,
                          min_abs_corr = 0)
  fitted$variant_id = rownames(gt_matrix)
  return(fitted)
}

extractResults <- function(susie_object){
  credible_sets = susie_object$sets$cs
  cs_list = list()
  susie_object$sets$purity = dplyr::as_tibble(susie_object$sets$purity) %>%
    dplyr::mutate(
      cs_id = rownames(susie_object$sets$purity),
      cs_size = NA,
      cs_log10bf = NA,
      overlapped = NA
    )
  added_variants = c()
  for (index in seq_along(credible_sets)){
    cs_variants = credible_sets[[index]]
    cs_id = susie_object$sets$cs_index[[index]]

    is_overlapped = any(cs_variants %in% added_variants)
    susie_object$sets$purity$overlapped[index] = is_overlapped
    susie_object$sets$purity$cs_size[index] = length(cs_variants)
    susie_object$sets$purity$cs_log10bf[index] = log10(exp(susie_object$lbf[cs_id]))
    if (!is_overlapped) {
      cs_list[[index]] = dplyr::tibble(cs_id = paste0("L", cs_id),
                                       variant_id = susie_object$variant_id[cs_variants])
      added_variants = append(added_variants, cs_variants)
    }
  }
  df = purrr::map_df(cs_list, identity)

  #Extract purity values for all sets
  purity_res = susie_object$sets$purity

  #Sometimes all the PIP values are 0 and there are no purity values, then skip this step
  if(nrow(purity_res) > 0){
    purity_df = dplyr::as_tibble(purity_res) %>%
      dplyr::filter(!overlapped) %>%
      dplyr::mutate(
        cs_avg_r2 = mean.abs.corr^2,
        cs_min_r2 = min.abs.corr^2,
        low_purity = min.abs.corr < 0.5
      )  %>%
      dplyr::select(cs_id, cs_log10bf, cs_avg_r2, cs_min_r2, cs_size, low_purity) 
  } else{
    purity_df = dplyr::tibble()
  }
  
  #Extract betas and standard errors and lbf_variables
  mean_vec = susieR::susie_get_posterior_mean(susie_object)
  sd_vec = susieR::susie_get_posterior_sd(susie_object)
  alpha_mat = t(susie_object$alpha)
  colnames(alpha_mat) = paste0("alpha", seq(ncol(alpha_mat)))
  mean_mat = t(susie_object$alpha * susie_object$mu) / susie_object$X_column_scale_factors
  colnames(mean_mat) = paste0("mean", seq(ncol(mean_mat)))
  sd_mat = sqrt(t(susie_object$alpha * susie_object$mu2 - (susie_object$alpha * susie_object$mu)^2)) / (susie_object$X_column_scale_factors)
  colnames(sd_mat) = paste0("sd", seq(ncol(sd_mat)))
  lbf_variable_mat = t(susie_object$lbf_variable)
  colnames(lbf_variable_mat) = paste0("lbf_variable", seq(ncol(lbf_variable_mat)))
  posterior_df = dplyr::tibble(variant_id = names(mean_vec), 
                               pip = susie_object$pip,
                               z = susie_object$z,
                               posterior_mean = mean_vec, 
                               posterior_sd = sd_vec) %>%
                 dplyr::bind_cols(purrr::map(list(alpha_mat, mean_mat, sd_mat, lbf_variable_mat), dplyr::as_tibble))

  if(nrow(df) > 0 & nrow(purity_df) > 0){
    cs_df = purity_df
    variant_df = dplyr::left_join(posterior_df, df, by = "variant_id") %>%
      dplyr::left_join(cs_df, by = "cs_id")
  } else{
    cs_df = NULL
    variant_df = NULL
  }

  return(list(cs_df = cs_df, variant_df = variant_df))
}

#Import all files
expression_matrix = readr::read_tsv(opt$expression_matrix)
sample_metadata = utils::read.csv(opt$sample_meta, sep = '\t', stringsAsFactors = F)
phenotype_meta = utils::read.csv(opt$phenotype_meta, sep = "\t", stringsAsFactors = F)
covariates_matrix = importQtlmapCovariates(opt$covariates)

#Exclude covariates with zero variance
exclude_cov = apply(covariates_matrix, 2, sd) != 0
covariates_matrix = covariates_matrix[,exclude_cov]

#Import list of phenotypes for finemapping
if (opt$permuted == "true"){
  phenotype_table = importQtlmapPermutedPvalues(opt$phenotype_list)
  filtered_list = dplyr::filter(phenotype_table, p_fdr < 0.01)
  phenotype_list = dplyr::semi_join(phenotype_meta, filtered_list, by = "group_id")
  message("Number of phenotypes included for analysis: ", nrow(phenotype_list))
} else {
  nominal_leads = read.table(opt$phenotype_list, stringsAsFactors = FALSE) %>%
    dplyr::transmute(phenotype_id = V1, n_variants = V6, p_nominal= V12) %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(phenotype_id) %>%
    dplyr::mutate(p_bonferroni = p.adjust(p_nominal, n = n_variants)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(p_fdr = p.adjust(p_bonferroni, method = "fdr")) %>%
    dplyr::filter(p_fdr < 0.1) #More lenient threshold for bonferroni correction
  phenotype_list = dplyr::semi_join(phenotype_meta, nominal_leads, by = "phenotype_id")
  message("Number of phenotypes included for analysis: ", nrow(phenotype_list))
}

#Keep only those phenotypes that are present in the expression matrix
phenotype_list = dplyr::filter(phenotype_list, phenotype_id %in% expression_matrix$phenotype_id)

#Set parameters
cis_distance = opt$cisdistance
genotype_file = opt$genotype_matrix
study_id = sample_metadata$study[1]

#Make a SummarizedExperiment of the expression data
se = eQTLUtils::makeSummarizedExperimentFromCountMatrix(assay = expression_matrix, 
                                                         row_data = phenotype_meta, 
                                                         col_data = sample_metadata, 
                                                         quant_method = "gene_counts",
                                                         reformat = FALSE)

#If qtl_group is not specified, then use the first value in the qtl_group column of the sample metadata
if(is.null(opt$qtl_group)){
  opt$qtl_group = se$qtl_group[1]
}

#Split phenotype list into chunks
chunk_vector = strsplit(opt$chunk, split = " ") %>% unlist() %>% as.numeric()
chunk_id = chunk_vector[1]
n_chunks = chunk_vector[2]
selected_chunk = splitIntoChunks(chunk_id, n_chunks, length(phenotype_list$phenotype_id))
selected_phenotypes = phenotype_list$phenotype_id[selected_chunk] %>% setNames(as.list(.), .)

#Only proceed if the there are more than 0 phenotypes
if(length(selected_phenotypes) > 0){
  #Check that the qtl_group is valid and subset
  assertthat::assert_that(opt$qtl_group %in% unique(se$qtl_group))
  selected_qtl_group = eQTLUtils::subsetSEByColumnValue(se, "qtl_group", opt$qtl_group)
  
  #Apply finemapping to all genes
  results = purrr::map(selected_phenotypes, ~finemapPhenotype(., selected_qtl_group, 
                                                              genotype_file, covariates_matrix, cis_distance))
  
  #Define fine-mapped regions
  region_df = dplyr::transmute(phenotype_list, phenotype_id, finemapped_region = paste0("chr", chromosome, ":", 
                                                                                        phenotype_pos - cis_distance, "-",
                                                                                        phenotype_pos + cis_distance))
  #Extract credible sets from finemapping results
  res = purrr::map(results, extractResults) %>%
    purrr::transpose()
  
  #Extract information about all variants
  variant_df <- purrr::map_df(res$variant_df, identity, .id = "phenotype_id")
  if(nrow(variant_df) > 0){
    variant_df <- variant_df %>%
      dplyr::left_join(region_df, by = "phenotype_id") %>%
      tidyr::separate(variant_id, c("chr", "pos", "ref", "alt"),sep = "_", remove = FALSE) %>%
      dplyr::mutate(chr = stringr::str_remove_all(chr, "chr")) %>%
      dplyr::mutate(cs_index = cs_id) %>%
      dplyr::mutate(cs_id = paste(phenotype_id, cs_index, sep = "_"))
  }
  
  #Extraxt information about credible sets
  cs_df <- purrr::map_df(res$cs_df, identity, .id = "phenotype_id")
} else { #Define empty data frames
  cs_df = dplyr::tibble()
  variant_df = dplyr::tibble()
}

if(nrow(cs_df) > 0){
  cs_df = dplyr::left_join(cs_df, region_df, by = "phenotype_id") %>%
    dplyr::mutate(cs_index = cs_id) %>%
    dplyr::mutate(cs_id = paste(phenotype_id, cs_index, sep = "_")) %>%
    dplyr::transmute(molecular_trait_id = phenotype_id, cs_id, cs_index, finemapped_region, cs_log10bf, cs_avg_r2, cs_min_r2, cs_size, low_purity)
  
  #Extract information about variants that belong to a credible set
  in_cs_variant_df <- dplyr::filter(variant_df, !is.na(cs_index) & !low_purity) %>%
    dplyr::transmute(molecular_trait_id = phenotype_id, variant = variant_id, chromosome = chr, position = pos, 
                     ref, alt, cs_id, cs_index, finemapped_region, pip, z, cs_min_r2, cs_avg_r2, cs_size, posterior_mean, posterior_sd, cs_log10bf)
} else{
  #Initialize empty tibbles with correct column names
  in_cs_variant_df = dplyr::tibble(
    molecular_trait_id = character(),
    variant = character(),
    chromosome = numeric(),
    position = numeric(),
    ref = numeric(),
    alt = numeric(),
    cs_id = numeric(),
    cs_index = numeric(),
    finemapped_region = numeric(),
    pip = numeric(),
    z = numeric(),
    cs_min_r2 = numeric(),
    cs_avg_r2 = numeric(),
    cs_size = numeric(),
    posterior_mean = numeric(),
    posterior_sd = numeric(),
    cs_log10bf = numeric()
  )
  
  cs_df = dplyr::tibble(
    molecular_trait_id = numeric(),
    cs_id = numeric(),
    cs_index = numeric(),
    finemapped_region = numeric(),
    cs_log10bf = numeric(),
    cs_avg_r2 = numeric(),
    cs_min_r2 = numeric(),
    cs_size = numeric(),
    low_purity = numeric()
  )
}

#Extract information about all variants
if(nrow(variant_df) > 0){
  variant_df_transmute <- dplyr::transmute(variant_df, molecular_trait_id = phenotype_id, variant = variant_id, chromosome = chr, position = pos, ref, alt, cs_id, cs_index, low_purity, finemapped_region, pip, z, posterior_mean, posterior_sd)  
  variant_df <- dplyr::bind_cols(variant_df_transmute, dplyr::select(variant_df,alpha1:lbf_variable10))
} else{
  variant_df = dplyr::tibble(
    molecular_trait_id = character(),
    variant = character(),
    chromosome = numeric(),
    position = numeric(),
    ref = character(),
    alt = character(),
    cs_id = character(),
    cs_index = character(),
    low_purity = character(),
    finemapped_region = character(),
    pip = numeric(),
    z = numeric(),
    posterior_mean = numeric(),
    posterior_sd = numeric(),
    alpha1 = numeric(),
    alpha2 = numeric(),
    alpha3 = numeric(),
    alpha4 = numeric(),
    alpha5 = numeric(),
    alpha6 = numeric(),
    alpha7 = numeric(),
    alpha8 = numeric(),
    alpha9 = numeric(),
    alpha10 = numeric(),
    mean1 = numeric(),
    mean2 = numeric(),
    mean3 = numeric(),
    mean4 = numeric(),
    mean5 = numeric(),
    mean6 = numeric(),
    mean7 = numeric(),
    mean8 = numeric(),
    mean9 = numeric(),
    mean10 = numeric(),
    sd1 = numeric(),
    sd2 = numeric(),
    sd3 = numeric(),
    sd4 = numeric(),
    sd5 = numeric(),
    sd6 = numeric(),
    sd7 = numeric(),
    sd8 = numeric(),
    sd9 = numeric(),
    sd10 = numeric(),
    lbf_variable1 = numeric(),
    lbf_variable2 = numeric(),
    lbf_variable3 = numeric(),
    lbf_variable4 = numeric(),
    lbf_variable5 = numeric(),
    lbf_variable6 = numeric(),
    lbf_variable7 = numeric(),
    lbf_variable8 = numeric(),
    lbf_variable9 = numeric(),
    lbf_variable10 = numeric()
  )
}

#Export high purity credible set results only
write.table(in_cs_variant_df, paste0(opt$out_prefix, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#Export all other results
write.table(cs_df, paste0(opt$out_prefix, ".cred.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(variant_df, paste0(opt$out_prefix, ".snp.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

