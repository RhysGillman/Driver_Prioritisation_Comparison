# run_DawnRank.r
# This master scripts uses the scripts within scripts/DawnRank to run DawnRank on the validation data

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(snowfall, quietly = T))
suppressPackageStartupMessages (library(maxstat, quietly = T))


# Handling input arguments
option_list = list(
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network to use", metavar ="Network"),
  make_option(c("-c", "--celltype"), type="character", default="Liver", 
              help="cell type to analyse", metavar ="Network")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

network_choice <- opt$network
cell_type <- opt$celltype

#############################
# Functions
#############################

fill_matrix <- function(matrix){
  filled_matrix <- matrix
  for(col_id in rownames(matrix)[!rownames(matrix) %in% colnames(matrix)]){
    column <- data.frame(rep(0, nrow(matrix)))
    colnames(column) <- col_id
    filled_matrix <- cbind(filled_matrix, column)
  }
  for(row_id in colnames(filled_matrix)[!colnames(filled_matrix) %in% rownames(filled_matrix)]){
    row <- data.frame(rep(0, ncol(filled_matrix)))
    colnames(row) <- row_id
    row <- row %>% t() %>% as.data.frame()
    colnames(row) <- colnames(filled_matrix)
    filled_matrix <- rbind(filled_matrix, row)
  }
  filled_matrix <- filled_matrix[,order(colnames(filled_matrix))]
  filled_matrix <- filled_matrix[order(rownames(filled_matrix)),]
  return(filled_matrix)
}


#############################
# Sample Info
#############################

sample_info <- read_csv(paste0("validation_data/CCLE_",network_choice,"/sample_info.csv")) %>% filter(lineage==cell_type)
samples <- sample_info$cell_ID

#############################
# Prepare Input Data
#############################

rna <- fread(paste0("validation_data/CCLE_",network_choice,"/tpm.csv"), select = c("gene_ID",samples)) %>% 
  arrange(gene_ID) %>%
  column_to_rownames("gene_ID") %>% 
  as.matrix()

mutation <-  fread(paste0("validation_data/CCLE_",network_choice,"/mutations.csv"), select = c("gene_ID",samples)) %>% 
  arrange(gene_ID) %>%
  column_to_rownames("gene_ID") %>% 
  as.matrix()

network <-  fread(paste0("validation_data/CCLE_",network_choice,"/network_directed.csv")) %>%
  dplyr::select(1,2) %>%
  table() %>%
  as.data.frame.matrix()

network <-  fread(paste0("validation_data/CCLE_",network_choice,"/network_directed.csv")) %>%
  pivot_wider(id_expand = T, names_from = colnames(network)[2], values_from = colnames(network)[3], values_fill = 0)


load("../DawnRank/data/brcaExamplePathway.rda")
isSymmetric(brcaExamplePathway)
