importVariantInformationFromGDS <- function(gdsfile){
  
  #Import individual columns
  snp_pos = GDSArray::GDSArray(gdsfile, "snp.position")
  snp_chromosome = GDSArray::GDSArray(gdsfile, "snp.chromosome")
  snp_id = GDSArray::GDSArray(gdsfile, "snp.rs.id")
  
  #Make a data frame
  snp_df = dplyr::data_frame(gds_snp_id = as.integer(names(snp_id)), 
                             chromosome = as.vector(snp_chromosome), 
                             pos = as.vector(snp_pos), 
                             snp_id = as.vector(snp_id))
  return(snp_df)
}

extractVariantGenotypeFromGDS <- function(variant_id, variant_information, gdsfile){
  assertthat::assert_that(length(variant_id) == 1)
  
  #Extract variant gds_id from variant infromation
  var_filtered = dplyr::filter(variant_information, snp_id == variant_id)
  selected_gds_id = var_filtered$gds_snp_id
  assertthat::assert_that(length(selected_gds_id) == 1)
  
  #Extract genotype from the gds file
  geno = GDSArray::GDSArray(gdsfile, "genotype")
  genotype = geno[selected_gds_id,]
  genotype_df = dplyr::data_frame(genotype_id = colnames(geno), genotype_value = genotype, snp_id = variant_id)
  return(genotype_df)
}