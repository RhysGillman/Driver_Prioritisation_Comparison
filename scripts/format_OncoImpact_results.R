#format_OncoImpact_results.R
# This code reformats the output from OncoImpact to be consistent with the other algorithms

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))


# Handling input arguments
option_list = list(
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network to use", metavar ="Network"),
  make_option(c("-c", "--celltype"), type="character", default="Liver", 
              help="cell type to analyse", metavar ="Cell Type")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

network_choice <- opt$network
cell_type <- opt$celltype



suppressWarnings(rm(OncoImpact_results))

for(result_file in list.files(paste0("results/CCLE_",network_choice,"/OncoImpact/",cell_type,"/sample_driver_list/"))){
  cell_ID <- gsub("\\.txt", "", result_file)
  result <- read_tsv(paste0("results/CCLE_",network_choice,"/OncoImpact/",cell_type,"/sample_driver_list/", result_file)) %>%
    dplyr::select(gene_ID = `#GENE`, impact = SAMPLE_IMPACT) %>%
    mutate(cell_ID) %>%
    relocate(cell_ID)
  if(!exists("OncoImpact_results")){
    OncoImpact_results <- result
  }else{
    OncoImpact_results <- rbind(OncoImpact_results, result)
  }
}

OncoImpact_results <- OncoImpact_results %>%
  arrange(desc(impact)) %>%
  group_by(cell_ID) %>%
  mutate(rank = row_number()) %>%
  dplyr::select(-impact) %>%
  arrange(cell_ID) %>%
  dplyr::mutate(lineage = cell_type) %>%
  dplyr::select(lineage, cell_ID, driver = gene_ID, rank)



write_csv(OncoImpact_results, paste0("results/CCLE_",network_choice,"/OncoImpact/",cell_type,".csv"))
