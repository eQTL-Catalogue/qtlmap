suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  optparse::make_option(c("--cs_results"), type="character", default=NULL,
                        help="Path to fine mapping credible sets file. Tab separated", metavar = "type"),
  optparse::make_option(c("--sumstats"), type="character", default=NULL,
                        help="Path the summary statistics file. Tab separated", metavar = "type"),
  optparse::make_option(c("--output"), type="character", default=NULL,
                        help="Path to the output file.", metavar = "type")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

#Import data
credible_sets = read.table(opt$cs_results, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = "\t") %>%
  dplyr::as_tibble()
sumstats = read.table(opt$sumstats, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = "\t") %>%
  dplyr::as_tibble()

#Merge the two tables
cs_table = dplyr::transmute(credible_sets, molecular_trait_id, variant, cs_id, pip, cs_size, z, cs_min_r2, finemapped_region) %>% 
  dplyr::left_join(sumstats, by = c("molecular_trait_id", "variant")) %>%
  dplyr::select(molecular_trait_id, gene_id, cs_id, variant, rsid, cs_size, pip, pvalue, beta, se, z, cs_min_r2, finemapped_region)

#Save file to disk
file_handle = gzfile(opt$putput,"w")
write.table(cs_table, file_handle, sep = "\t", row.names = F, col.names = T, quote = FALSE)
close(file_handle)



