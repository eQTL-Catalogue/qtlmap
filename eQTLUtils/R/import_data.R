#' A general function to quickly import tabix indexed tab-separated files into data_frame
#'
#' @param tabix_file Path to tabix-indexed text file
#' @param param A instance of GRanges, RangedData, or RangesList
#' provide the sequence names and regions to be parsed. Passed onto Rsamtools::scanTabix()
#' @param ... Additional parameters to be passed on to readr::read_delim()
#'
#' @return List of data_frames, one for each entry in the param GRanges object.
#' @export
scanTabixDataFrame <- function(tabix_file, param, ...){
  tabix_list = Rsamtools::scanTabix(tabix_file, param = param)
  df_list = lapply(tabix_list, function(x,...){
    if(length(x) > 0){
      if(length(x) == 1){
        #Hack to make sure that it also works for data frames with only one row
        #Adds an empty row and then removes it
        result = paste(paste(x, collapse = "\n"),"\n",sep = "")
        result = readr::read_delim(result, delim = "\t", ...)[1,]
      }else{
        result = paste(x, collapse = "\n")
        result = readr::read_delim(result, delim = "\t", ...)
      }
    } else{
      #Return NULL if the nothing is returned from tabix file
      result = NULL
    }
    return(result)
  }, ...)
  return(df_list)
}


#Import transcript metadata from biomart web export
importBiomartMetadata <- function(biomart_path){
  transcript_meta = readr::read_tsv(biomart_path)
  col_df = dplyr::data_frame(column_name = c('Gene stable ID', 'Transcript stable ID', 'Chromosome/scaffold name', 'Gene start (bp)', 'Gene end (bp)', 'Strand', 'Transcript start (bp)', 'Transcript end (bp)', 'Transcription start site (TSS)', 'Transcript length (including UTRs and CDS)', 'Transcript support level (TSL)', 'APPRIS annotation', 'GENCODE basic annotation', 'Gene name', 'Transcript name', 'Transcript count', 'Transcript type', 'Gene type', 'Gene % GC content', 'Version (gene)', 'Version (transcript)'),
                             column_id = c('gene_id', 'transcript_id', 'chromosome', 'gene_start', 'gene_end', 'strand', 'transcript_start', 'transcript_end', 'tss', 'transcript_length', 'transcript_tsl', 'transcript_appris', 'is_gencode_basic', 'gene_name', 'transcript_name', 'transcript_count', 'transcript_type', 'gene_type', 'gene_gc_content', 'gene_version', 'transcript_version'))
  transcript_meta = transcript_meta[,col_df$column_name] %>% dplyr::distinct()
  colnames(transcript_meta) = col_df$column_id
  return(transcript_meta)
}

#' Import variant information extracted from VCF file into R
#'
#' The variant information text file can be generated from the VCF using the following
#' bcftools command:
#' bcftools query -f '\%CHROM\\t\%POS\\t\%ID\\t\%REF\\t\%ALT\\t\%TYPE\\t\%AC\\t\%AN\\n' path/to/vcf_file.vcf.gz | bgzip > path/to/variant_infromation_file.txt.gz
#'
#' @param path Path to the the variant information text file.
#' @export
importVariantInformation <- function(path){
  info_col_names = c("chr","pos","snp_id","ref","alt","type","AC","AN", "MAF", "R2")
  into_col_types = "cicccciidd"
  snp_info = readr::read_delim(path, delim = "\t", col_types = into_col_types, col_names = info_col_names)
  snp_info = dplyr::mutate(snp_info, indel_length = pmax(nchar(alt), nchar(ref))) %>%
    dplyr::mutate(is_indel = ifelse(indel_length > 1, TRUE, FALSE)) %>%
    dplyr::mutate(MAF = pmin(AC/AN, 1-(AC/AN)))
  return(snp_info)
}
