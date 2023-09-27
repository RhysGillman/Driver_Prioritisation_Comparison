# format_sysSVM2_results.R

# This code reformats the output from PersonaDrive to be consistent with the other algorithms

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(biomaRt, quietly = T))

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


# Prepare entrez ID mapping file

if(!file.exists(paste0("validation_data/CCLE_", network_choice,"/entrez_ids.csv"))){
  genes <- data.table::fread(paste0("validation_data/CCLE_",network_choice,"/cnv.csv"), select = c("gene_ID")) %>% pull(gene_ID) %>% unique()
  mart <- biomaRt::useMart("ensembl")
  mart <- useDataset("hsapiens_gene_ensembl",mart=mart)
  attr <- listAttributes(mart)
  entrez_id_mapping <- getBM(attributes=c('hgnc_symbol','entrezgene_id'), mart = mart, filters = 'hgnc_symbol', values = genes)
  write_csv(entrez_id_mapping, paste0("validation_data/CCLE_", network_choice,"/entrez_ids.csv"))
}else{
  entrez_id_mapping <- read_csv(paste0("validation_data/CCLE_", network_choice,"/entrez_ids.csv"))
}

# Get sample list from results files
samples <- list.files(paste0("results/CCLE_",network_choice,"/sysSVM2/",cell_type,"/")) 
samples <- gsub(".csv","", samples)

suppressWarnings(rm(sysSVM2_results))

for(s in samples){
  cell_drivers <- read_csv(paste0("results/CCLE_",network_choice,"/sysSVM2/",cell_type, "/",s,".csv"), col_types = cols()) %>%
    dplyr::select(cell_ID=sample,entrez_ID=entrez,rank=sample_rank)
  if(exists("sysSVM2_results")){
    sysSVM2_results <- rbind(sysSVM2_results, cell_drivers)
  }else{
    sysSVM2_results <- cell_drivers
  }
}

sysSVM2_results <- sysSVM2_results %>%
  inner_join(entrez_id_mapping, by = c("entrez_ID"="entrezgene_id"), multiple="all") %>%
  dplyr::select(cell_ID,gene_ID=hgnc_symbol,rank) %>%
  group_by(cell_ID) %>%
  filter(!duplicated(gene_ID)) %>%
  arrange(rank) %>%
  mutate(rank=row_number()) %>%
  ungroup() %>%
  arrange(cell_ID, rank) %>%
  mutate(lineage=cell_type) %>%
  dplyr::select(lineage, cell_ID, driver=gene_ID, rank)

write_csv(sysSVM2_results, paste0("results/CCLE_",network_choice,"/sysSVM2/",cell_type, ".csv"))
