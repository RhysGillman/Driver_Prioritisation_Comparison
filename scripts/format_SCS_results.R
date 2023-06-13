#format_SCS_results.R
# This code reformats the output from SCS to be consistent with the other algorithms

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



rm(SCS_results)
results_names <- read_csv(paste0("results/CCLE_",network_choice,"/SCS/",cell_type,"/sample_names.txt")) %>% deframe()
for(result in names(results_names)){
  path <- paste0("results/CCLE_",network_choice,"/SCS/",cell_type,"/result_sample_",result,".csv")
  if(file.exists(path)){
    tmp <- read_csv(path, skip = 1, col_names = c("gene_ID", "module", "impact")) %>%
      arrange(desc(impact)) %>%
      mutate(rank = row_number(), cell_ID=results_names[result])
    if(!exists("SCS_results")){
      SCS_results <- tmp
    }else{
      SCS_results <- rbind(SCS_results, tmp)
    }
  }
  
}

SCS_results <- SCS_results %>%
  dplyr::mutate(lineage = cell_type) %>%
  arrange(cell_ID) %>%
  dplyr::select(lineage, cell_ID, driver = gene_ID, rank)

write_csv(SCS_results, paste0("results/CCLE_",network_choice,"/SCS/",cell_type,".csv"))