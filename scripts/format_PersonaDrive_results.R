#format_PersonaDrive_results.R
# This code reformats the output from PersonaDrive to be consistent with the other algorithms

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

samples <- read_tsv(paste0("results/CCLE_",network_choice,"/PersonaDrive/",cell_type, "/PersonaDrive.txt"), col_names = FALSE) %>% pull(X1)


rm(PersonaDrive_results)

for(i in 1:length(samples)){
  cell_drivers <- read_tsv(paste0("results/CCLE_",network_choice,"/PersonaDrive/",cell_type, "/PersonaDrive.txt"), col_names = FALSE, skip = (i-1), n_max = 1)
  colnames(cell_drivers) <- c("cell_ID", 1:(ncol(cell_drivers)-1))
  cell_drivers <- cell_drivers %>%
    gather("rank", "driver",-cell_ID) %>%
    na.omit()
  if(exists("PersonaDrive_results")){
    PersonaDrive_results <- rbind(PersonaDrive_results, cell_drivers)
  }else{
    
    PersonaDrive_results <- cell_drivers
  }
}

PersonaDrive_results <- PersonaDrive_results %>%
  mutate(lineage=cell_type) %>%
  dplyr::select(lineage, cell_ID, driver, rank)

write_csv(PersonaDrive_results, paste0("results/CCLE_",network_choice,"/PersonaDrive/",cell_type, ".csv"))


