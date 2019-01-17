qtltoolsPrepareSE <- function(se, quant_method){
  require("cqn")
  
  #Specify valid chromsomes and valid gene types
  valid_gene_types = c("lincRNA","protein_coding","IG_C_gene","IG_D_gene","IG_J_gene",
                       "IG_V_gene", "TR_C_gene","TR_D_gene","TR_J_gene", "TR_V_gene",
                       "3prime_overlapping_ncrna","known_ncrna", "processed_transcript",
                       "antisense","sense_intronic","sense_overlapping")
  valid_chromosomes = c("1","10","11","12","13","14","15","16","17","18","19",
                        "2","20","21","22","3","4","5","6","7","8","9")
  
  #Normalise featureCounts data
  if(quant_method == "featureCounts"){
    #Filter SE to keep correct chromosomes and QC-passed samples
    se_filtered = filterSummarizedExperiment(se, valid_chromosomes = valid_chromosomes,
                                             valid_gene_types = valid_gene_types,
                                             filter_rna_qc = TRUE, filter_genotype_qc = TRUE)
    
    #Identify expressed genes (At least 1 count in 10% of the samples)
    se_expressed = filterSE_expressedGenes(se_filtered, min_count = 1, min_fraction = 0.1, assay_name = "counts")
    
    #Normalise and make QTLtools matrix
    se_norm = normaliseSE_cqn(se_expressed, assay_name = "counts")
 
    } else if(quant_method == "array"){
    #Filter SE to keep correct chromosomes and QC-passed samples
    message(" ## Filterin out invalid chromosomes, RNA_QC failed samples and genotype_QC failed samples ")
    se_filtered = eQTLUtils::filterSummarizedExperiment(se, valid_chromosomes = valid_chromosomes, 
                                                            filter_rna_qc = TRUE, filter_genotype_qc = TRUE)
    
    #Normalize and regress out batch effects
    message(" ## Normalizing and regressing out the batch effects")
    se_norm = eQTLUtils::array_normaliseSE(se_filtered, norm_method = "quantile", assay_name = "exprs", 
                                           log_transform = TRUE, adjust_batch = TRUE, filter_quality = TRUE)
  }

  return(se_norm)
}