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

message(" ## Reading credible sets.")
credible_sets = read.table(opt$cs_results, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = "\t") %>%
  dplyr::as_tibble()
message(" ## Read ", nrow(credible_sets), " rows of credible sets.")

if (nrow(credible_sets) == 0) {
  credible_sets = dplyr::tibble(
    molecular_trait_id = numeric(),
    gene_id = numeric(),
    cs_id = numeric(),
    variant = numeric(),
    rsid = numeric(),
    cs_size = numeric(),
    pip = numeric(),
    pvalue = numeric(),
    beta = numeric(),
    se = numeric(),
    z = numeric(),
    cs_min_r2 = numeric(),
    finemapped_region = numeric()
  )
  
  #Save file to disk
  file_handle = gzfile(opt$output,"w")
  write.table(credible_sets, file_handle, sep = "\t", row.names = F, col.names = T, quote = FALSE)
  close(file_handle)
  
  message("Credible sets matrix is empty. Writing it as an output and stopping execution!")
  quit(save = "no", status = 0)
}

message(" ## Reading summary statistics.")
sumstats = read.table(opt$sumstats, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = "\t") %>%
  dplyr::as_tibble() %>%
  dplyr::filter(molecular_trait_id %in% credible_sets$molecular_trait_id)
message(" ## Read ", nrow(sumstats), " rows of summary statistics")

message(" ## Merge the two tables")
cs_table = dplyr::transmute(credible_sets, molecular_trait_id, variant, cs_id, pip, cs_size, z, cs_min_r2, finemapped_region) %>% 
  dplyr::left_join(sumstats, by = c("molecular_trait_id", "variant")) %>%
  dplyr::select(molecular_trait_id, gene_id, cs_id, variant, rsid, cs_size, pip, pvalue, beta, se, z, cs_min_r2, finemapped_region)
message(" ## After merge there are ", nrow(sumstats), " rows of merged matrix")

#Save file to disk
file_handle = gzfile(opt$output,"w")
write.table(cs_table, file_handle, sep = "\t", row.names = F, col.names = T, quote = FALSE)
close(file_handle)



