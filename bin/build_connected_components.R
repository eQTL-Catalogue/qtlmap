message(" ## Loading libraries: optparse")
suppressPackageStartupMessages(library("optparse"))

#Parse command-line options
option_list <- list(
  #TODO look around if there is a package recognizing delimiter in dataset
  make_option(c("-s", "--susie_file_path"), type="character", default=NULL,
              help="Purity filtered susie output. Tab separated file", metavar = "type"),
  make_option(c("-p", "--phenotype_meta"), type="character", default=NULL,
              help="Phenotype metadata file. Tab separated file", metavar = "type"),
  make_option(c("-q", "--qtl_group"), type="character", default=NULL,
              help="qtl group of the given dataset", metavar = "type"),
  make_option(c("-o", "--output_file"), type="character", default="./cc_lead_phenotypes.txt",
              help="Output file path and name [default \"%default\"]", metavar = "type"),
  make_option(c("-c", "--cs_size_threshold"), type="integer", default=200,
              help="Threshold to filter out credible sets bigger than given size", metavar = "type")
)

message(" ## Parsing options")
opt <- optparse::parse_args(OptionParser(option_list=option_list))

message(" ## Loading libraries: devtools, dplyr, SummarizedExperiment, cqn, data.table")
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))

make_connected_components_from_cs <- function(susie_all_df, z_threshold = 3, cs_size_threshold = 10) {
  # Filter the credible sets by Z-score and size
  susie_filt <- susie_all_df %>%
    dplyr::group_by(cs_uid) %>%
    dplyr::mutate(max_abs_z = max(abs(z))) %>%
    dplyr::filter(max_abs_z > z_threshold, cs_size < cs_size_threshold) %>%
    dplyr::ungroup()
  
  # make the ranges object in order to find overlaps
  cs_ranges = GenomicRanges::GRanges(
    seqnames = susie_filt$chromosome,
    ranges = IRanges::IRanges(start = susie_filt$position, end = susie_filt$position),
    strand = "*",
    mcols = data.frame(cs_uid = susie_filt$cs_uid, variant_id = susie_filt$variant, gene_id = susie_filt$molecular_trait_id)
  )
  
  # find overlaps and remove the duplicated
  olaps <-  GenomicRanges::findOverlaps(cs_ranges, cs_ranges, ignore.strand = TRUE) %>%
    GenomicRanges::as.data.frame() %>%
    dplyr::filter(queryHits <= subjectHits)
  
  # change variant sharing into credible set sharing 
  # not to have multiple connected components of variants but credible sets
  olaps <- olaps %>% dplyr::mutate(cs_uid_1 = cs_ranges$mcols.cs_uid[queryHits], cs_uid_2 = cs_ranges$mcols.cs_uid[subjectHits])
  edge_list <- olaps %>% dplyr::select(cs_uid_1, cs_uid_2) %>% BiocGenerics::unique() %>% base::as.matrix()
  
  # make the graph of connected components
  g <- igraph::graph_from_edgelist(edge_list, directed = F)
  g_cc <- igraph::components(g)
  
  # turn connected components graph into data frame
  cc_df <- data.frame(cc_membership_no = g_cc$membership, 
                      cs_uid = g_cc$membership %>% names()) %>% 
    dplyr::mutate(molecular_trait_id = base::sub(pattern = "_[^_]+$", replacement = "",x = base::gsub(".*\\%","",cs_uid)))
  
  return(cc_df)
}

#Debugging
if (FALSE) {
  opt = list()
  opt$s = "/Users/kerimov/Work/temp_files/leafcutter_data/Alasoo_2018_17Feb/susie/Alasoo_2018_leafcutter_macrophage_naive.purity_filtered.txt.gz"
  opt$p = "/Users/kerimov/Work/temp_files/leafcutter_data/Alasoo_2018_17Feb/leafcutter_metadata.txt.gz"
  opt$q = "Alasoo_2018_leafcutter_macrophage_IFNg"
  opt$o = "/Users/kerimov/Work/temp_files/leafcutter_data/Alasoo_2018_17Feb/cc_connected.tsv"
  opt$c = 200
}

susie_file_path = opt$s
phenotype_meta_path = opt$p
output_dir = opt$o
qtl_group_in = opt$q
cs_size_threshold_in = opt$c

message("######### Options: ######### ")
message("######### Working Directory  : ", getwd())
message("######### qtl_group          : ", qtl_group_in)
message("######### susie_file_path    : ", susie_file_path)
message("######### cs_size_threshold_in   : ", cs_size_threshold_in)
message("######### phenotype_meta_path: ", phenotype_meta_path)
message("######### output_dir         : ", output_dir)

message(" ## Reading susie file")
susie_naive <- readr::read_tsv(file = susie_file_path, col_types = "cccicccccdddddddd")

message(" ## Reading leafcutter metadata file")
# leafcutter_metadata <- readr::read_tsv(phenotype_meta_path) %>% 
#   dplyr::filter(!is.na(gene_id)) %>% 
#   dplyr::filter(gene_count == 1)

message(" ## Building Connected-Components")
susie_naive <- susie_naive %>% dplyr::mutate(cs_uid = paste0(qtl_group_in, "%", cs_id))
susie_naive_cc <- make_connected_components_from_cs(susie_all_df = susie_naive, cs_size_threshold = cs_size_threshold_in)

susie_naive_filt_cc <- susie_naive %>% 
  dplyr::filter(molecular_trait_id %in% susie_naive_cc$molecular_trait_id) %>% 
  dplyr::left_join(susie_naive_cc %>% dplyr::select(cc_membership_no, molecular_trait_id), by = "molecular_trait_id")

susie_naive_high_pip_var <- susie_naive_filt_cc %>% 
  dplyr::group_by(cc_membership_no) %>% 
  dplyr::arrange(-pip) %>% 
  dplyr::slice(1) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(-pip)

# susie_high_pip_with_gene <- susie_naive_high_pip_var %>% 
#   dplyr::left_join(leafcutter_metadata %>% dplyr::select(-chromosome), by = c("molecular_trait_id" = "phenotype_id")) %>% 
#   dplyr::filter(!is.na(gene_id))

message(" ## Writing: ", output_dir)
readr::write_tsv(x = susie_naive_high_pip_var %>% dplyr::select(molecular_trait_id) %>% dplyr::distinct(), file = output_dir, col_names = F)
