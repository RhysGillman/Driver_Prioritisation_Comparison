#format_PNC_results.R
# This code reformats the output from PNC to be consistent with the other algorithms

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

#############################
# Sample Info
#############################

sample_info <- read_csv(paste0("validation_data/CCLE_",network_choice,"/sample_info.csv")) %>% filter(lineage==cell_type)
samples <- sample_info$cell_ID %>% sort()


#############################
# Genetically Altered Genes
#############################
# Mutation

mutation <- fread(paste0("validation_data/CCLE_",network_choice,"/mutations.csv"), select = c("gene_ID",samples)) %>% 
  arrange(gene_ID) %>%
  column_to_rownames("gene_ID")

# CNV

cnv <- fread(paste0("validation_data/CCLE_",network_choice,"/cnv.csv"), select = c("gene_ID",samples)) %>%
  column_to_rownames("gene_ID")

# Get a list of altered genes for each sample

mutated <- which(mutation==1, arr.ind = TRUE) %>% 
  as.data.frame() %>%
  mutate(sample = colnames(mutation)[col]) %>%
  mutate(gene_ID = rownames(mutation)[row]) %>%
  dplyr::select(sample, gene_ID)

copy_altered <- which(cnv!=0, arr.ind = TRUE) %>% 
  as.data.frame() %>%
  mutate(sample = colnames(cnv)[col]) %>%
  mutate(gene_ID = rownames(cnv)[row]) %>%
  dplyr::select(sample, gene_ID)

altered_genes <- rbind(mutated,copy_altered)
rownames(altered_genes) <- NULL
altered_genes <- unique(altered_genes) %>%
  mutate(is_altered = TRUE)


#############################
# PNC Results
#############################

PNC <- read_csv(paste0("results/CCLE_",network_choice,"/PNC/",cell_type,"/result.csv"), show_col_types = FALSE) %>% 
  column_to_rownames("Row")

PNC_drivers <- which(PNC==1, arr.ind = TRUE) %>% 
  as.data.frame() %>%
  mutate(sample = colnames(PNC)[col]) %>%
  mutate(gene_ID = rownames(PNC)[row]) %>%
  dplyr::select(sample, gene_ID)

rownames(PNC_drivers) <- NULL

PNC_drivers <- PNC_drivers %>%
  unique() %>%
  left_join(altered_genes, by = c("sample", "gene_ID")) %>%
  mutate(is_altered = ifelse(is.na(is_altered), FALSE, is_altered))

out_deg <- read_csv(paste0("results/CCLE_",network_choice,"/PNC/",cell_type,"/out_deg.csv"),show_col_types = FALSE)
out_deg[is.na(out_deg)] <- 0
out_deg <- out_deg %>% pivot_longer(!gene_ID, names_to = "sample", values_to = "out_degree")

in_deg <- read_csv(paste0("results/CCLE_",network_choice,"/PNC/",cell_type,"/in_deg.csv"),show_col_types = FALSE)
in_deg[is.na(in_deg)] <- 0
in_deg <- in_deg %>% pivot_longer(!gene_ID, names_to = "sample", values_to = "in_degree")
degree <- full_join(out_deg,in_deg, by = c("gene_ID", "sample"))


PNC_results <- PNC_drivers %>% 
  left_join(degree, by = c("gene_ID", "sample")) %>%
  filter(is_altered) %>%
  arrange(desc(out_degree)) %>%
  group_by(sample) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  dplyr::mutate(lineage = cell_type) %>%
  arrange(sample, rank) %>%
  select(lineage, cell_ID = sample, driver = gene_ID, rank)

write_csv(PNC_results, paste0("results/CCLE_",network_choice,"/PNC/",cell_type,".csv"))