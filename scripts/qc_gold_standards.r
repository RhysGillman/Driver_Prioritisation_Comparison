suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(readxl, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(ggpubr, quietly = T))





# Handling input arguments
option_list = list(
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network being used", metavar ="Network")
  
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


network_choice <- opt$network



cell_list <- read_csv("data/cell_list.csv")

#############################
# Get cell line information
#############################

CCLE_sample_info <- read_csv("data/CCLE/Model.csv") %>%
  dplyr::select(DepMap_ID=ModelID, lineage=OncotreeLineage, cell_ID=StrippedCellLineName) %>%
  mutate(lineage=gsub(" |/","_",lineage)) %>%
  filter(cell_ID %in% cell_list$cell_ID) %>%
  unique()



#############################
# Data
#############################

gene_effect <- fread("data/CCLE/CRISPRGeneEffect.csv") %>%
  as.data.frame() %>%
  dplyr::rename(DepMap_ID=V1) %>%
  inner_join(CCLE_sample_info[c("DepMap_ID","cell_ID")], by = "DepMap_ID") %>%
  dplyr::select(-DepMap_ID) %>%
  column_to_rownames("cell_ID")

colnames(gene_effect) <- gsub(" \\(.*","", colnames(gene_effect))
gene_effect <- gene_effect[!duplicated(colnames(gene_effect))]


gold_standard <- read_csv(paste0("validation_data/CCLE_",network_choice,"/gold_standards.csv")) %>%
  dplyr::select(cell_ID,gene_ID) %>%
  unique()

#############################
# Check Normality
#############################


scaled_gene_effect <- gene_effect %>%
  rownames_to_column("cell_ID") %>%
  pivot_longer(cols = -cell_ID, names_to = "gene_ID", values_to = "gene_effect") %>%
  group_by(gene_ID) %>%
  mutate(scaled_gene_effect= (gene_effect - min(gene_effect))/(max(gene_effect) - min(gene_effect))  )

ggplot(scaled_gene_effect, aes(x=scaled_gene_effect)) +
  geom_histogram(bins = 20) +
  theme_minimal() +
  xlab("Min-Max Scaled Gene Effect") +
  ylab("Count")

ggqqplot(scaled_gene_effect$scaled_gene_effect)


#############################
# Check Gold Standards
#############################


gs_sample <- gold_standard %>%
  filter(gene_ID %in% sample(gold_standard$gene_ID, 12, replace = F)) %>% 
  mutate(is_gold_standard=T)

check_gene_effect_gs <- gene_effect[,unique(gs_sample$gene_ID)] %>% 
  as.data.frame() %>%
  rownames_to_column("cell_ID") %>%
  pivot_longer(-cell_ID, names_to = "gene_ID", values_to = "gene_effect") %>%
  left_join(gs_sample, by = c("cell_ID","gene_ID")) %>%
  mutate(is_gold_standard = ifelse(is.na(is_gold_standard), FALSE, is_gold_standard)) %>%
  left_join(CCLE_sample_info %>% dplyr::select(cell_ID, lineage), by = "cell_ID")

ggplot(check_gene_effect_gs, aes(x=gene_effect)) +
  geom_density() +
  geom_point(aes(y=is_gold_standard, colour = is_gold_standard)) +
  facet_wrap(~gene_ID, scales = "free")


