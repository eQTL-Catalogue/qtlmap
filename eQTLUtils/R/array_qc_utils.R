# Microarray data QC helpers
# by Liis Kolberg

# Filter SE for gene types and/or chromosome

array_filterSE_gene_types <- function(se, valid_chromosomes = c("1","10","11","12","13","14","15","16","17","18","19",
                                                          "2","20","21","22","3","4","5","6","7","8","9","MT","X","Y"),
                                valid_gene_types = c("lincRNA","protein_coding","IG_C_gene","IG_D_gene","IG_J_gene",
                                                     "IG_V_gene", "TR_C_gene","TR_D_gene","TR_J_gene", "TR_V_gene",
                                                     "3prime_overlapping_ncrna","known_ncrna", "processed_transcript",
                                                     "antisense","sense_intronic","sense_overlapping")){
  selected_genes = dplyr::filter(as.data.frame(rowData(se)),
                                 gene_type %in% valid_gene_types, chromosome %in% valid_chromosomes)

  return(se[as.character(selected_genes$phenotype_id),])
}

#' Generate PCA QC Plot
#'
#' @param se SummarizedExperiment file to be analysed
#' @param assay_name Assay name to be analysed (exprs, norm_exprs or batch_exprs)
#' @param filter_quality Boolean value if samples with RNA_QC_passed = FALSE should be excluded (Default:TRUE)
#' @param show_quality Boolean value if samples with RNA_QC_passed = FALSE should be highlighted (Default:FALSE)
#' @param html_output Boolean value if html output should be created (Default:FALSE)
#' @param output_dir html file output dir, if html_output is TRUE (Default:current directory)
#' @return PCA plot of study
#' @author Liis Kolberg
#' @export
#'
array_plotPCA <- function(se, assay_name = "exprs", filter_quality = TRUE, show_quality = FALSE, html_output = FALSE, output_dir = "./"){
  processed_se = se[, se$cell_type %in% c("monocytes", "neutrophils", "T-cell", "B-cell", "platelet")]
  processed_se = processed_se %>% array_filterSE_gene_types(valid_gene_types="protein_coding")

  if(filter_quality){
    processed_se = processed_se[,processed_se$RNA_QC_passed == TRUE]
  }

  study_name <- processed_se$study %>% unique()
  study_name <- paste0(study_name, collapse="_")

  # Extract fields
  row_data = SummarizedExperiment::rowData(processed_se)
  col_data = SummarizedExperiment::colData(processed_se)
  assays = SummarizedExperiment::assays(processed_se)

  expr_matrix = assays[[assay_name]]
  if(assay_name=="exprs"){
    expr_matrix = log2_transform(expr_matrix)
  }
  expr_matrix = as.matrix(expr_matrix)
  pca_res = prcomp(t(expr_matrix), center = T, scale. = T, rank. = 10)

  pca_matrix = as.data.frame(pca_res$x) %>%
    as_tibble() %>%
    dplyr::mutate(sample_id = rownames(pca_res$x)) %>%
    dplyr::left_join(as.data.frame(col_data), by = "sample_id")

  pca_matrix$grp = paste(pca_matrix$cell_type, pca_matrix$marker, pca_matrix$study)
  pca_matrix$grp_condition = paste(pca_matrix$condition, pca_matrix$timepoint)

  conditions = table(pca_matrix$grp_condition)
  types = table(pca_matrix$grp)

  if(length(types) > length(conditions)){
    pcaplt <- ggplot2::ggplot(pca_matrix, aes(x = PC1, y = PC2, color = grp, shape = grp_condition, label = sample_id)) + geom_point() + scale_shape_manual(values=seq(0,6)) + theme_bw()
  }
  else{
    pcaplt <- ggplot2::ggplot(pca_matrix, aes(x = PC1, y = PC2, color = grp_condition, shape = grp, label = sample_id)) + geom_point() + scale_shape_manual(values=seq(0,6)) + theme_bw()
  }
  if(show_quality){
    pcaplt <- pcaplt + geom_point(data = pca_matrix[pca_matrix$RNA_QC_passed==FALSE, ], aes(x = PC1, y = PC2), color="darkred", show.legend = FALSE)
  }
  PCA_ggplotly_plot = plotly::ggplotly()

  if (html_output) {
    if (!dir.exists(output_dir)) { dir.create(output_dir) }
    htmlwidgets::saveWidget(plotly::as_widget(PCA_ggplotly_plot), file.path(output_dir, paste0(study_name,"_", assay_name, "_PCA", ".html")))
  }
  return(PCA_ggplotly_plot)
}

#' Generate Boxplot QC Plot for n random samples
#'
#' @param se SummarizedExperiment file to be analysed
#' @param n Number of random samples (Default:10)
#' @param assay_name Assay name to be analysed (exprs, norm_exprs or batch_exprs)
#' @param filter_quality Boolean value if samples with RNA_QC_passed = FALSE should be excluded (Default:TRUE)
#' @param pdf_output Boolean value if pdf output should be created (Default:FALSE)
#' @param output_dir html file output dir, if html_output is TRUE (Default:current directory)
#' @return Boxplot plot of sample expression values
#' @author Liis Kolberg
#' @export
#'
array_Boxplot <- function(se, n = 10, assay_name = "exprs", filter_quality = TRUE, pdf_output = FALSE, output_dir = "./"){
  processed_se = se[, se$cell_type %in% c("monocytes", "neutrophils", "T-cell", "B-cell", "platelet")]
  processed_se = processed_se %>% array_filterSE_gene_types(valid_gene_types="protein_coding")
  if(filter_quality){
    processed_se = processed_se[,processed_se$RNA_QC_passed == TRUE]
  }

  study_name <- processed_se$study %>% unique()
  study_name <- paste0(study_name, collapse="_")
  mat <- SummarizedExperiment::assays(processed_se)[[assay_name]]
  if(assay_name=="exprs"){
    mat = log2_transform(mat)
    ylab_name = "log2(exprs)"
  }else{
    ylab_name = "normalized(log2(exprs))"
  }

  mat <- mat[,sample(ncol(mat), n)]

  meltData <- reshape::melt(mat)
  names(meltData) <- c("gene", "sample_id", "value")
  meltData$sample_id <- as.character(meltData$sample_id)
  box_dat <- meltData %>%
    dplyr::left_join(as.data.frame(colData(processed_se)), by = "sample_id")

  bxplot <- ggplot2::ggplot(box_dat, aes(x=sample_id, y=value, fill=marker)) + geom_boxplot() + xlab("Sample ID") +
    ylab(ylab_name) + ggtitle(study_name) + theme_bw()

  if (pdf_output) {
    if (!dir.exists(output_dir)) { dir.create(output_dir) }
    ggsave(filename=file.path(output_dir, paste0(study_name,"_", assay_name, "_boxplot", ".pdf")))
  }
  return(bxplot)
}

#' Log2 transform intensity values (same approach as in function lumiB)
#'
#' @param mat Matrix of intensity values
#' @return Log2 transformed matrix
#' @author Liis Kolberg
#' @export
log2_transform <- function(mat){
  # forcePositives
  offset <- apply(mat, 2, min, na.rm = TRUE)
  offset[offset <= 0] <- offset[offset <= 0] - 1.01
  offset[offset > 0] <- 0
  offset <- rep(1, nrow(mat)) %*% t(offset)
  mat <- mat - offset
  return(log2(mat))
}

#' Normalise microarray data using lumiN function
#'
#' @param se SummarizedExperiment file to be normalised
#' @param norm_method Normalisation method from quantile, rsn, ssn, loess, vsn, rankinvariant (Default:quantile)
#' @param assay_name Assay name to be normalised (exprs or batch_exprs)
#' @param log_transform Boolean value if intensities should be log2 transformed (Default:TRUE)
#' @param adjust_batch Boolean value if data should be adjusted for batch effects before normalisation (Default:FALSE)
#' @param filter_quality Boolean value if should exclude samples with RNA_QC_passed=FALSE (Default:TRUE)
#' @return SummarizedExperiment with assay 'norm_exprs'
#' @author Liis Kolberg
#' @importFrom dplyr "%>%"
#' @export
#'
array_normaliseSE <- function(se, norm_method = "quantile", assay_name = "exprs", log_transform = TRUE, adjust_batch = FALSE, filter_quality = TRUE){
  valid_methods <- c("quantile", "rsn", "ssn", "loess", "vsn", "rankinvariant") # Quantile normalization, Robust spline normalization

  if (!(norm_method %in% valid_methods))
    stop(paste0("Not a valid method. Please use value from the following list: ", paste0(valid_methods, collapse=", ")))

  processed_se = se
  if(filter_quality){
    processed_se = processed_se[,processed_se$rna_qc_passed==TRUE]
  }
  # Extract fields
  row_data = SummarizedExperiment::rowData(processed_se)
  col_data = SummarizedExperiment::colData(processed_se)
  assays = SummarizedExperiment::assays(processed_se)

  # Batch adjustment
  if(adjust_batch & length(table(col_data$batch)) > 1){
    adj_se = array_RemoveBatch(processed_se, assay_name=assay_name)
    assays = SummarizedExperiment::assays(adj_se)
    expr_matrix = assays[["batch_exprs"]]
  }
  else{
    expr_matrix = assays[[assay_name]]
    if(log_transform){
      expr_matrix = expr_matrix %>% log2_transform()
    }
  }
  # Perform normalisation
  expr_matrix = as.matrix(expr_matrix)
  expr_norm = lumi::lumiN(expr_matrix, method = norm_method)

  # Update assays
  assays[["norm_exprs"]] = expr_norm

  # Make an update se object
  processed_se = SummarizedExperiment::SummarizedExperiment(
    assays = assays,
    colData = col_data,
    rowData = row_data)
  return(processed_se)
}

#' Generate Sex dependent QC Plot.
#'
#' @param study_data SummarizedExperiment file to be analysed
#' @param assay_name Assay name to be analysed (exprs, norm_exprs or batch_exprs)
#' @param html_output Boolean value if html output should be created (Default:FALSE)
#' @param output_dir html file output dir, if html_output is TRUE (Default:current directory)
#' @return plotly object of SEX QC
#' @author Nurlan Kerimov, Liis Kolberg
#' @export
#'
array_PlotSexQC <- function(study_data, assay_name = "exprs", html_output=FALSE, output_dir="./"){
  joined <- array_CalculateSexQCDataFrame(study_data, assay_name=assay_name)

  if (all(is.na(joined$sex))) { joined$sex <- "NA"}

  study_name <- study_data$study %>% unique()

  # generate the plot
  Sex_QC_plot <- ggplot2::ggplot(joined, ggplot2::aes(x=XIST, y=Y_chrom_mean, label = sample_id)) + ggplot2::geom_point(ggplot2::aes(col=sex)) + ggplot2::labs(x="Expression XIST", y="Mean Expression genes on Y", title = paste0(study_name, " DS | Sample Size: ", nrow(joined))) + expand_limits(y=0, x=0) + ggplot2::theme_bw()
  Sex_QC_ggplotly_plot <- plotly::ggplotly()

  if (html_output) {
    if (!dir.exists(output_dir)) { dir.create(output_dir) }
    htmlwidgets::saveWidget(plotly::as_widget(Sex_QC_ggplotly_plot), file.path(output_dir, paste0(study_name, assay_name,  "_XIST", ".html")))
  }

  return(Sex_QC_ggplotly_plot)
}

#' Generate Sex dependent QC DataFrame.
#'
#' @param study_data SummarizedExperiment file to be analysed
#' @param assay_name Assay name to be analysed (exprs, norm_exprs or batch_exprs)
#' @return DataFrame object of SEX QC
#' @author Nurlan Kerimov, Liis Kolberg
#' @export
#'
array_CalculateSexQCDataFrame <- function(study_data, assay_name){
  # get the rowData of SummarizedExperiment
  study_rowdata <- SummarizedExperiment::rowData(study_data) %>% SummarizedExperiment::as.data.frame()

  # keep only naive samples
  study_data_filt = study_data[,study_data$condition=="naive" & study_data$RNA_QC_passed==TRUE & study_data$cell_type %in% c("monocytes", "neutrophils", "T-cell", "B-cell", "platelet")]

  # get the Y chromosome genes only from rowdata
  study_Y_chrom_data <- dplyr::filter(study_rowdata, chromosome=="Y", gene_type=="protein_coding") %>% dplyr::arrange(gene_start)

  # get expressions of Y chromosome related genes

  counts_df <- study_data_filt %>% SummarizedExperiment::assay() %>% SummarizedExperiment::as.data.frame()

  study_Y_gene_counts <- counts_df %>% subset(rownames(counts_df) %in% study_Y_chrom_data$phenotype_id)

  # Probe ID of XIST gene
  XIST_gene = as.character(study_rowdata[study_rowdata$gene_id=="ENSG00000229807",][["phenotype_id"]])

  # get XIST gene (ENSG00000229807) and all Y chromosome genes
  selected_genes <- study_Y_gene_counts %>% rownames() %>% c(XIST_gene)

  exprs_data <- (study_data_filt %>% SummarizedExperiment::assays())[[assay_name]]

  exprs <- exprs_data[intersect(rownames(exprs_data), selected_genes),] %>% SummarizedExperiment::as.data.frame()
  exprs_xist <- exprs[XIST_gene,]
  exprs_Y <- exprs[-which(rownames(exprs) %in% XIST_gene),]

  #get the mean of Y chromosome related gene expressions
  study_Y_gene_mean_exprs <- exprs_Y %>% apply(2, mean) %>% SummarizedExperiment::as.data.frame()
  study_Y_gene_mean_exprs <- study_Y_gene_mean_exprs %>% mutate(sample_id = rownames(study_Y_gene_mean_exprs))
  colnames(study_Y_gene_mean_exprs)[which(names(study_Y_gene_mean_exprs) == ".")] <- c("Y_chrom_mean")

  exprs_xist <- t(exprs_xist) %>% SummarizedExperiment::as.data.frame()
  colnames(exprs_xist) = c("XIST")
  exprs_xist <- exprs_xist %>% mutate(sample_id = rownames(exprs_xist))


  # join Y_chrom and Xist datasets for plot
  joined <- inner_join(study_Y_gene_mean_exprs, exprs_xist)
  joined <- study_data_filt %>% SummarizedExperiment::colData() %>% SummarizedExperiment::as.data.frame() %>% select(sex, sample_id) %>% inner_join(joined)

  return(joined)
}

#' Generate MDS QC Plot
#'
#' @param se SummarizedExperiment file to be analysed
#' @param assay_name Assay name to be analysed (exprs, norm_exprs or batch_exprs)
#' @param condition List of conditions to plot, e.g naive (Default:NULL)
#' @param filter_quality Boolean value if samples with RNA_QC_passed = FALSE should be excluded (Default:TRUE)
#' @param html_output Boolean value if html output should be created (Default:FALSE)
#' @param output_dir html file output dir, if html_output is TRUE (Default:current directory)
#' @return MDS plot of study
#' @author Liis Kolberg
#' @export
#'

array_PlotMDS <- function(se, assay_name="exprs", condition = NULL, filter_quality = TRUE, html_output = FALSE, output_dir = "./"){
  processed_se = se[,se$cell_type %in% c("monocytes", "neutrophils", "T-cell", "B-cell", "platelet")]
  if(!is.null(condition)){
    processed_se = processed_se[, processed_se$condition %in% condition]
  }
  if(filter_quality){
    processed_se = processed_se[, processed_se$RNA_QC_passed == TRUE]
  }
  processed_se = processed_se %>% array_filterSE_gene_types(valid_gene_types="protein_coding")

  study_name <- processed_se$study %>% unique()
  study_name <- paste0(study_name, collapse="_")

  matrix = SummarizedExperiment::assays(processed_se)[[assay_name]]
  if(assay_name == "exprs"){
    matrix = log2_transform(matrix)
  }
  # Perform MDS
  dist = cor(matrix, method = "pearson")

  fit <- MASS::isoMDS(1 - dist, k=2)

  mds_matrix = as.data.frame(fit$points) %>%
    as_tibble() %>%
    dplyr::mutate(sample_id = rownames(fit$points)) %>%
    dplyr::left_join(as.data.frame(colData(processed_se)), by = "sample_id")

  mds_matrix$grp = paste(mds_matrix$cell_type, mds_matrix$marker)
  mds_matrix$grp_condition = paste(mds_matrix$condition, mds_matrix$timepoint)

  conditions = table(mds_matrix$grp_condition)
  types = table(mds_matrix$grp)

  if(length(types) > length(conditions)){
    mds_plot = ggplot2::ggplot(mds_matrix, aes(x = V1, y = V2, color = grp, shape = grp_condition, label = sample_id)) +
      geom_point(size = 3) + scale_shape_manual(values=seq(0,6)) +
      xlab("MDS Coordinate 1") +
      ylab("MDS Coordinate 2") + theme_bw()
  }else{
    mds_plot = ggplot2::ggplot(mds_matrix, aes(x = V1, y = V2, color = grp_condition, shape = grp, label = sample_id)) +
      geom_point(size = 3) + scale_shape_manual(values=seq(0,6)) +
      xlab("MDS Coordinate 1") +
      ylab("MDS Coordinate 2") + theme_bw()
  }


  MDS_ggplotly_plot = plotly::ggplotly()
  if (html_output) {
    if (!dir.exists(output_dir)) { dir.create(output_dir) }
    htmlwidgets::saveWidget(plotly::as_widget(MDS_ggplotly_plot), file.path(output_dir, paste0(study_name,"_", assay_name, "_MDS", ".html")))
  }
  return(MDS_ggplotly_plot)
}

#' Adjust for batch effects using limma removeBatchEffect
#'
#' @param se SummarizedExperiment file to be analysed
#' @param assay_name Assay name to be analysed (exprs, norm_exprs or batch_exprs)
#' @return Batch effect adjusted SE with new assay called batch_exprs
#' @author Liis Kolberg
#' @export
#'

array_RemoveBatch <- function(se, assay_name="exprs"){
  #processed_se = se[,se$RNA_QC_passed==TRUE & se$cell_type %in% c("monocytes", "neutrophils", "T-cell", "B-cell", "platelet") ]
  processed_se = se
  batches = factor(colData(processed_se)$batch)

  mat = as.matrix(assays(processed_se)[[assay_name]])
  # removeBatchEffect assumes log transformed data
  mat = log2_transform(mat)

  sample_metadata = SummarizedExperiment::colData(processed_se)
  nr_of_cell_types = table(sample_metadata$cell_type)
  # merge condition and timepoint
  condition_timepoint = paste(sample_metadata$condition, sample_metadata$timepoint)
  nr_of_conditions = table(condition_timepoint)

  if(dim(nr_of_cell_types) > 1 & dim(nr_of_conditions) > 1){
    design = model.matrix(~sample_metadata$cell_type + condition_timepoint)
  }
  else if(dim(nr_of_cell_types) > 1){
      design = model.matrix(~sample_metadata$cell_type)
    }
  else if(dim(nr_of_conditions) > 1){
    design = model.matrix(~condition_timepoint)
  }
  else{
    design = matrix(1, ncol(mat), 1)
  }
  mat_wobatch = limma::removeBatchEffect(mat, batch = batches, design = design)
  SummarizedExperiment::assays(processed_se)[["batch_exprs"]] = data.frame(mat_wobatch)
  processed_se = SummarizedExperiment::SummarizedExperiment(
    assays = assays(processed_se),
    colData = colData(processed_se),
    rowData = rowData(processed_se))
  return(processed_se)
}
