#' Generate MBV_results plots
#'
#' @param mbv_files_path Path where mbv_output files live
#' @param output_path Path where the plots will ve saved
#' @author Nurlan Kerimov
#' @export

plot_mbv_results <- function(mbv_files_path, output_path){
  mbv_files = list.files(mbv_files_path, full.names = T)
  
  #Make sample names
  sample_names = stringr::str_replace_all(basename(mbv_files), ".mbv_output.txt", "")
  sample_list = setNames(mbv_files, sample_names)
  
  #Import mbv files
  mbv_results_list = purrr::map(sample_list, ~readr::read_delim(., delim = " ", col_types = "ciiiiiiiiii"))
  
  for (i in c(1:length(mbv_results_list))) {
    df <- mbv_results_list[[i]]
    mbv_sample_name <- names(mbv_results_list)[i]
    
    df <- df %>% dplyr::mutate(het_consistent_frac = n_het_consistent/n_het_covered, 
                        hom_consistent_frac = n_hom_consistent/n_hom_covered,
                        match = FALSE) %>% 
                 dplyr::arrange(-het_consistent_frac)
    df$match[1] <- TRUE  
    
    plot_mbv <- ggplot2::ggplot(df, aes(x = het_consistent_frac, y = hom_consistent_frac, label = SampleID, color=match)) + 
      ggplot2::geom_point() + 
      ggplot2::scale_color_manual(values=c("red4","darkgreen"))
    
    if (!dir.exists(output_path)) { 
      dir.create(paste0(output_path, "/jpeg/"), recursive = TRUE)
      dir.create(paste0(output_path, "/plotly/dependencies/"), recursive = TRUE)
    }
    ggplot2::ggsave(filename = paste0(mbv_sample_name, "_mbv_plot.jpeg"), plot = plot_mbv, path = paste0(output_path, "/jpeg"), device = "jpeg")
    
    MDS_ggplotly_plot <- plotly::ggplotly(plot_mbv)
    htmlwidgets::saveWidget(widget = plotly::as_widget(MDS_ggplotly_plot), 
                            file = file.path(normalizePath(paste0(output_path, "/plotly")), paste0(mbv_sample_name, "_plotly.html")),
                            libdir = "dependencies")
  }
}


#' Generate Sex dependent QC Plot.
#'
#' @param study_data SummarizedExperiment file to be analysed
#' @param html_output Boolean value if html output should be created (Default:FALSE)
#' @param output_dir html file output dir, if html_output is TRUE (Default:current directory)
#' @return plotly object of SEX QC
#' @author Nurlan Kerimov
#' @export
plotSexQC <- function(study_data, html_output=FALSE, output_dir="./"){
  joined <- calculateSexQCDataFrame(study_data)
  
  if (all(is.na(joined$sex))) { joined$sex <- "NA"}
  
  study_name <- study_data$study %>% unique()
  # generate the plot
  Sex_QC_plot <- ggplot2::ggplot(joined, 
    ggplot2::aes(x=(ENSG00000229807+1) %>% log2(), y=(Y_chrom_mean+1) %>% log2(), label = sample_id)) + 
    ggplot2::geom_point(ggplot2::aes(col=sex)) +
    ggplot2::labs(x="Expression XIST", y="Expression genes on Y", title = paste0(study_name, " DS - TPM normalized, log2 | Sample Size: ", nrow(joined))) 
  MDS_ggplotly_plot <- plotly::ggplotly()

  if (html_output) {
    if (!dir.exists(output_dir)) { dir.create(output_dir) }
    htmlwidgets::saveWidget(plotly::as_widget(MDS_ggplotly_plot), file.path(normalizePath(output_dir), paste0(study_name, "_sex_QC.html")))
  }
  
  return(MDS_ggplotly_plot)
}


#' Generate Sex dependent QC DataFrame.
#'
#' @param study_data SummarizedExperiment file to be analysed
#' @return DataFrame object of SEX QC
#' @author Nurlan Kerimov
#' @export
calculateSexQCDataFrame <- function(study_data){
  # get the rowData of SummarizedExperiment
  study_rowdata <- SummarizedExperiment::rowData(study_data) %>% SummarizedExperiment::as.data.frame()
  XIST_Count_study <-study_data[study_rowdata$gene_id=="ENSG00000229807",] %>% SummarizedExperiment::assay()
  
  #get the Y chromosome genes only from rowdata
  study_Y_chrom_data <- dplyr::filter(study_rowdata, chromosome=="Y", gene_type=="protein_coding") %>% dplyr::arrange(gene_start)
  
  #get feature counts of Y chromosome related genes
  counts_df <- study_data %>% SummarizedExperiment::assay() %>% SummarizedExperiment::as.data.frame() 
  study_Y_gene_counts <- counts_df %>% subset(rownames(counts_df) %in% study_Y_chrom_data$gene_id)
  
  # get XIST gene (ENSG00000229807) and all Y chromosome genes
  selected_genes_for_TPM <- study_Y_gene_counts %>% rownames() %>% c("ENSG00000229807")
  
  #normalise the counts Transcripts Per Million (TPM)
  study_data <- eQTLUtils::normaliseSE_tpm(study_data)
  tpms_data <- (study_data %>% SummarizedExperiment::assays())$tpms 
  
  normalized_counts <- tpms_data[intersect(rownames(tpms_data), selected_genes_for_TPM),] %>% SummarizedExperiment::as.data.frame()
  normalized_counts_xist <- normalized_counts["ENSG00000229807",]
  normalized_counts_Y <- normalized_counts[-which(rownames(normalized_counts) %in% "ENSG00000229807"),]
  
  #get the mean of Y chromosome related gene feature counts (after tpm normlisation)
  study_Y_gene_mean_counts <- normalized_counts_Y %>% apply(2, mean) %>% SummarizedExperiment::as.data.frame() 
  study_Y_gene_mean_counts <- study_Y_gene_mean_counts %>% mutate(sample_id = rownames(study_Y_gene_mean_counts))
  colnames(study_Y_gene_mean_counts)[which(names(study_Y_gene_mean_counts) == ".")] <- c("Y_chrom_mean")
  
  normalized_counts_xist <- t(normalized_counts_xist) %>% SummarizedExperiment::as.data.frame()
  normalized_counts_xist <- normalized_counts_xist %>% mutate(sample_id = rownames(normalized_counts_xist))
  
  # join Y_chrom and Xist datasets for plot
  joined <- inner_join(study_Y_gene_mean_counts, normalized_counts_xist) 
  joined <- study_data %>% SummarizedExperiment::colData() %>% SummarizedExperiment::as.data.frame() %>% select(sex, sample_id) %>% inner_join(joined)
  
  return(joined)
}


#' Generate MDS QC Plot for different cell types
#'
#' @param study_data_se SummarizedExperiment file to be analysed
#' @param condition Boolean value if html output should be created (Default:FALSE)
#' @param html_output Boolean value if html output should be created (Default:FALSE)
#' @param output_dir html file output dir, if html_output is TRUE (Default:current directory)
#' @return MDS plot of study 
#' @author Nurlan Kerimov
#' @export
plotMDSAnalysis <- function(study_data_se, condition = "all", html_output=FALSE, output_dir="./"){
  mds_matrix = calculateMDSMatrix(study_data_se, condition)
  study_name <- study_data_se$study %>% unique()
  
  mds_plot = ggplot2::ggplot(mds_matrix, ggplot2::aes(x = V1, y = V2, color = cell_type, shape = study, label = sample_id)) + 
    ggplot2::geom_point() + 
    ggplot2::scale_shape_manual(values=seq(0,6)) +
    ggplot2::labs(x="Expression XIST", y="Expression genes on Y", title = paste0(study_name, " DS - TPM normalized, log2 | Sample Size: ", nrow(study_data_se %>% SummarizedExperiment::colData()))) 
  
  MDS_ggplotly_plot <- plotly::ggplotly()

  if (html_output) {
    if (!dir.exists(output_dir)) { dir.create(output_dir) }
    htmlwidgets::saveWidget(plotly::as_widget(MDS_ggplotly_plot), file.path(normalizePath(output_dir), paste0(study_name, "_MDS_plot.html")))
  }
  
  return(MDS_ggplotly_plot)
}


#' Generate MDS Matrix for SummarizedExperiment
#'
#' @param study_data_se SummarizedExperiment file to be analysed
#' @param condition Boolean value if html output should be created (Default:FALSE)
#' @return MDS Matrix of study 
#' @author Nurlan Kerimov
#' @export
calculateMDSMatrix <- function(study_data_se, condition = "all"){
  if (condition!="all") {
    study_data_se = study_data_se[,study_data_se$condition == condition]
  }
  
  # choose only protein coding genes and TPM normalise
  processed_se = eQTLUtils::filterSE_gene_types(study_data_se, valid_gene_types = "protein_coding") %>% eQTLUtils::normaliseSE_tpm()
  processed_se = processed_se[apply(SummarizedExperiment::assays(processed_se)$tpms, 1, median) > 1,]
  
  #Perform MDS
  matrix = log(SummarizedExperiment::assays(processed_se)$tpms+0.1,2)
  dist = cor(matrix, method = "pearson")
  fit <- MASS::isoMDS(1-dist, k=2)
  
  mds_matrix = SummarizedExperiment::as.data.frame(fit$points) %>%
    as_tibble() %>%
    dplyr::mutate(sample_id = rownames(fit$points)) %>%
    dplyr::left_join(SummarizedExperiment::as.data.frame(SummarizedExperiment::colData(processed_se)), by = "sample_id")
  
  return(mds_matrix)
}

#' Generate PCA QC Plot for different cell types
#'
#' @param study_data_se SummarizedExperiment file to be analysed
#' @param condition Boolean value if html output should be created (Default:FALSE)
#' @param html_output Boolean value if html output should be created (Default:FALSE)
#' @param output_dir html file output dir, if html_output is TRUE (Default:current directory)
#' @return PCA plot of study 
#' @author Nurlan Kerimov
#' @export
plotPCAAnalysis <- function(study_data_se, condition = "all", html_output=FALSE, output_dir="./"){
  pca_matrix = calculatePCAMatrix(study_data_se, condition)
  study_name <- study_data_se$study %>% unique()
  
  PCA.plot <- ggplot2::ggplot(pca_matrix, ggplot2::aes(x = PC1, y = PC2, color = cell_type, shape = study, label = sample_id)) + 
    ggplot2::geom_point() + 
    ggplot2::scale_shape_manual(values=seq(0,6)) + 
    ggplot2::labs(x="PC 1", y="PC 2", title = paste0(study_name, " PCA analysis | Sample Size: ", nrow(study_data_se %>% SummarizedExperiment::colData())))
  
  PCA_ggplotly_plot <- plotly::ggplotly()
  
  if (html_output) {
    if (!dir.exists(output_dir)) { dir.create(output_dir) }
    htmlwidgets::saveWidget(plotly::as_widget(PCA_ggplotly_plot), file.path(normalizePath(output_dir), paste0(study_name, "_PCA_plot.html")))
  }
  
  return(PCA_ggplotly_plot)
}


#' Generate PCA Matrix for SummarizedExperiment
#'
#' @param study_data_se SummarizedExperiment file to be analysed
#' @param condition if needed PCA matrix of specific condition
#' @return PCA Matrix of study 
#' @author Nurlan Kerimov
#' @export
calculatePCAMatrix <- function(study_data_se, condition = "all"){
  if (condition!="all") {
    study_data_se = study_data_se[,study_data_se$condition == condition]
  }
  
  # choose only protein coding genes and TPM normalise
  processed_se = eQTLUtils::filterSE_gene_types(study_data_se, valid_gene_types = "protein_coding") %>% eQTLUtils::normaliseSE_tpm()
  processed_se = processed_se[apply(SummarizedExperiment::assays(processed_se)$tpms, 1, median) > 1,]
  
  #Perform PCA
  pca_res = eQTLUtils::transformSE_PCA(processed_se, assay_name = "tpms", n_pcs = 10, log_transform = TRUE, center = TRUE, scale. = TRUE)
  return(pca_res$pca_matrix)
}