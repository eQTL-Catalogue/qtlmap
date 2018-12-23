#' Save a list of matrices into a suitable format for QTLTools
#'
#' Works with expression and covariates matrices.
#'
#' @param data_list list of matrices
#' @param output_dir relative path to the output dir
#' @param file_suffix suffix added to each file after their name in the list.
#' @return None
#' @author Kaur Alasoo
#' @export
saveQTLToolsMatrices <- function(data_list, output_dir, file_suffix = "bed", col_names = TRUE){

  #Check if the output dir exists and if not then create one
  if(!file.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
  }

  #Save each matrix as a separate  txt file
  for (sn in names(data_list)){
    file_path = file.path(output_dir, paste(sn, file_suffix, sep = "."))
    print(file_path)
    write.table(data_list[[sn]], file_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = col_names)
  }
}

#' Import QTLtools output table from permutation run into R.
#'
#'
#' @param file_path Path to the QTLtools output file.
#' @return data_frame containing gene_ids, snp ids and p-values.
#' @author Kaur Alasoo
#' @export
importQTLtoolsTable <- function(file_path){
  col_names = c("group_id", "pheno_chr", "pheno_start", "pheno_end", "strand", "phenotype_id", "group_size", "n_cis_snps",
                "distance", "snp_id", "snp_chr", "snp_start", "snp_end", "df", "dummy", "beta1",
                "beta2", "p_nominal","slope","p_perm","p_beta")
  col_types = "cciicciiicciiiddddddd"
  table = readr::read_delim(file_path, col_names = col_names, delim = " ", col_types = col_types) %>%
    dplyr::filter(!is.na(p_beta)) %>%
    dplyr::mutate(p_bonferroni = p_nominal*group_size*n_cis_snps) %>%
    dplyr::mutate(p_bonferroni = pmin(p_bonferroni,1)) %>%
    dplyr::mutate(p_fdr = p.adjust(p_beta, method = "fdr")) %>%
    dplyr::mutate(qvalue = qvalue::qvalue(p_beta)$qvalues) %>%
    dplyr::arrange(p_fdr)
  return(table)
}

#' Fetch particular genes from tabix indexed QTLtools output file.
#'
#' @param gene_ranges GRanges object with coordinates of the cis regions around genes.
#' @param tabix_file Tabix-indexed fastqtl output file.
#'
#' @return List of data frames containing QTLtools results for each gene.
#' @export
qtltoolsTabixFetchPhenotypes <- function(phenotype_ranges, tabix_file){

  #Assertions about input
  assertthat::assert_that(class(phenotype_ranges) == "GRanges")
  assertthat::assert_that(assertthat::has_name(GenomicRanges::elementMetadata(phenotype_ranges), "phenotype_id"))

  #Set column names for rasqual
  fastqtl_columns = c("phenotype_id","pheno_chr","pheno_start", "pheno_end",
                      "strand","n_snps", "distance", "snp_id", "snp_chr",
                      "snp_start", "snp_end", "p_nominal","beta", "is_lead")
  fastqtl_coltypes = "cciiciicciiddi"

  result = list()
  for (i in seq_along(phenotype_ranges)){
    selected_phenotype_id = phenotype_ranges[i]$phenotype_id
    print(i)
    tabix_table = scanTabixDataFrame(tabix_file, phenotype_ranges[i],
                                     col_names = fastqtl_columns, col_types = fastqtl_coltypes)[[1]] %>%
      dplyr::filter(phenotype_id == selected_phenotype_id)

    #Add additional columns
    result[[selected_phenotype_id]] = tabix_table
  }
  return(result)
}


#' Post-process QTLTools mbv results to find the best matching individual for each sample
#'
#' @param mbv_df Data frame with MBV results for one sequencing sample.
#'
#' @return Data frame with one row identifying the best matching invidivual for this sample
#' @export
mbvFindBestMatch <- function(mbv_df){
  res = dplyr::transmute(mbv_df, mbv_genotype_id = SampleID,
                         het_consistent_frac = n_het_consistent/n_het_covered,
                         hom_consistent_frac = n_hom_consistent/n_hom_covered)

  #Identify best het
  best_het = dplyr::arrange(res, -het_consistent_frac) %>% dplyr::filter(dplyr::row_number() == 1)
  other_het = dplyr::arrange(res, -het_consistent_frac) %>% dplyr::filter(dplyr::row_number() > 1)
  best_row = dplyr::mutate(best_het, het_min_dist = min(best_het$het_consistent_frac - other_het$het_consistent_frac),
                           hom_min_dist = min(best_het$hom_consistent_frac - other_het$hom_consistent_frac))

  #Compare against best hom
  best_hom = dplyr::arrange(res, -hom_consistent_frac) %>% dplyr::filter(dplyr::row_number() == 1)
  if(best_row$mbv_genotype_id != best_hom$mbv_genotype_id){
    best_row = dplyr::mutate(best_row, het_consistent_frac = as.numeric(NA), hom_consistent_frac = as.numeric(NA),
                             het_min_dist = as.numeric(NA), hom_min_dist = as.numeric(NA))
  }
  return(best_row)
}


#' Convert SummarizedExperiment object into a bed file suitable for QTLTools
#'
#' @param se SummarizedExperiment object
#' @param assay_name Assay of the SummarizedExperiment object used for QTL mapping
#'
#' @return Assay converted into bed format suitable for QTLTools
#' @export
convertSEtoQTLtools <- function(se, assay_name = "cqn"){

  #Extract rowData from the SE
  phenotype_data = rowData(se) %>%
    as.data.frame() %>%
    dplyr::as_tibble()

  #Make sure that all required columns are present
  assertthat::assert_that(assertthat::has_name(phenotype_data, "chromosome"))
  assertthat::assert_that(assertthat::has_name(phenotype_data, "phenotype_pos"))
  assertthat::assert_that(assertthat::has_name(phenotype_data, "phenotype_id"))
  assertthat::assert_that(assertthat::has_name(phenotype_data, "group_id"))
  assertthat::assert_that(!is.null(se$genotype_id))

  #Make genePos table for QTLTools
  pheno_data = dplyr::arrange(phenotype_data, chromosome, phenotype_pos) %>%
    dplyr::transmute(chromosome, left = phenotype_pos, right = phenotype_pos, phenotype_id, group_id, strand) %>%
    dplyr::rename_("#chr" = "chromosome") %>%
    dplyr::mutate(strand = ifelse(strand == 1, "+", "-"))

  #Exptract phenotype and rename columns according to genotype id
  assay = assays(se)[[assay_name]]
  colnames(assay) = se$genotype_id
  assay = round(assay, 3) #Round to three digits

  #Make QTLtools phenotype table
  res = dplyr::mutate(as.data.frame(assay), phenotype_id = rownames(assay)) %>%
    dplyr::select(phenotype_id, dplyr::everything()) %>%
    dplyr::left_join(pheno_data, ., by = "phenotype_id") %>%
    dplyr::arrange()

  return(res)
}


importQTLtoolsPCA <- function(pca_path){
  naive_pca = readr::read_delim(pca_path, delim = " ") %>%
    dplyr::select(-SampleID)
  sample_ids = colnames(naive_pca)
  pca_df = t(naive_pca) %>%
    as.data.frame() %>%
    dplyr::as_tibble()
  colnames(pca_df) = paste0("PC", 1:ncol(pca_df))
  pca_df$sample_id = sample_ids
  pca_df = dplyr::select(pca_df, sample_id, everything())
  return(pca_df)
}

studySEtoQTLTools <- function(se, assay_name, out_dir){

  #Make assertions
  assertthat::assert_that(assertthat::has_name(SummarizedExperiment::colData(se), "qtl_group"))
  assertthat::assert_that(assertthat::has_name(SummarizedExperiment::assays(se), assay_name))

  #Split the SE into list based on qtl_group
  qtl_groups = unique(se$qtl_group)
  group_list = setNames(as.list(qtl_groups), qtl_groups)
  group_se_list = purrr::map(group_list, ~subsetSEByColumnValue(se, "qtl_group", .))

  #Convert SE onbjects to QTLtools
  qtltools_list = purrr::map(group_se_list, ~convertSEtoQTLtools(., assay_name = assay_name))
  saveQTLToolsMatrices(qtltools_list, output_dir = out_dir, file_suffix = "bed")

  #Extract sample names
  sample_names = purrr::map(qtltools_list, ~colnames(.)[-(1:6)])
  saveQTLToolsMatrices(sample_names, output_dir = out_dir, file_suffix = "sample_names.txt", col_names = FALSE)
}
