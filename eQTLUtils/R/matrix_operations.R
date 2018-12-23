calculateMean <- function(matrix, design, factor, sample_id_col = "sample_id", na.rm = FALSE){
  #Calculate the mean value in matrix over all possible factor values.

  #If the factor is not a factor then make it a factor.
  if(!is.factor(design[,factor])){
    design[,factor] = factor(design[,factor])
  }

  #Set sample_id column as rownames
  rownames(design) = design[,sample_id_col]
  factor = design[,factor]
  levs = levels(factor)
  result = c()
  for (lev in levs){
    filter = factor == lev
    samples = rownames(design[filter,])
    mat = matrix[,samples]
    mat = rowMeans(mat, na.rm)
    result = cbind(result, mat)
  }
  colnames(result) = levs
  return(data.frame(result))
}

zScoreNormalize <- function(matrix){
  #Normalize expression matrix by z-score
  matrix = matrix - rowMeans(matrix)
  matrix = matrix / apply(matrix, 1, sd)
  return(matrix)
}

replaceNAsWithRowMeans <- function(matrix){
  #replace with row means
  na_pos = which(is.na(matrix), arr.ind = TRUE)
  matrix[na_pos] = rowMeans(matrix, na.rm=TRUE)[na_pos[,1]]

  #If there are addional NAs left (whole row NAs) then replace with 0
  matrix[is.na(matrix)] = 0
  return(matrix)
}

#' Force a vector of values into standard normal distribution
#'
#' @param x numeric vector with arbitrary distribution
#'
#' @return Vector with a standard normal distribution
#' @export
quantileNormaliseVector = function(x){
  qnorm(rank(x,ties.method = "random")/(length(x)+1))
}

quantileNormaliseMatrix <- function(matrix){
  quantile_matrix = matrix(0, nrow(matrix), ncol(matrix))
  for (i in seq_along(matrix[1,])){
    quantile_matrix[,i] = quantileNormaliseVector(matrix[,i])
  }
  #Add names
  rownames(quantile_matrix) = rownames(matrix)
  colnames(quantile_matrix) = colnames(matrix)
  return(quantile_matrix)
}

quantileNormaliseCols <- function(matrix,...){
  quantileNormaliseMatrix(matrix, ...)
}

quantileNormaliseRows <- function(matrix,...){
  t(quantileNormaliseMatrix(t(matrix), ...))
}

calculateTranscriptRatios <- function(expression_matrix, gene_transcript_map){

  #Check that gene_transcript_map has the correct columns
  assertthat::assert_that(assertthat::has_name(gene_transcript_map, "gene_id"))
  assertthat::assert_that(assertthat::has_name(gene_transcript_map, "transcript_id"))
  assertthat::assert_that(length(colnames(gene_transcript_map)) == 2)

  #Add gene ids to transcript expression matrix
  tpm_df = dplyr::mutate(as.data.frame(expression_matrix), transcript_id = rownames(expression_matrix)) %>%
    tbl_df() %>%
    dplyr::left_join(gene_transcript_map, by = "transcript_id") %>%
    dplyr::select(-transcript_id) %>%
    dplyr::select(gene_id, everything())

  #Calculate total expression per gene
  gene_total_tpms = purrrlyr::slice_rows(tpm_df, "gene_id") %>%
    purrrlyr::by_slice(~colSums(.) %>%
                      t() %>%
                      as.data.frame(), .collate = "rows")

  #Make matrix of gene expression values for each transcript
  tx_gene_expression = dplyr::left_join(gene_transcript_map, gene_total_tpms, by = "gene_id")
  tx_gene_matrix = dplyr::select(tx_gene_expression, -gene_id, -transcript_id) %>% as.matrix()
  rownames(tx_gene_matrix) = tx_gene_expression$transcript_id

  #calculate TPM ratios and add them to the SummarizedExperiment
  tpm_ratios = expression_matrix/tx_gene_matrix[rownames(expression_matrix),]
  return(tpm_ratios)
}
