#!/usr/bin/env Rscript --vanilla

# choose_cells.r
# This script is run before running the driver prioritisation pipeline

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(DESeq2, quietly = T))
suppressPackageStartupMessages (library(ggrepel, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))

# Handling input arguments
option_list = list(
  make_option(c("-w", "--workdir"), type="character", default=NULL, 
              help="working directory", metavar ="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$workdir)){
  print_help(opt_parser)
  stop("Please provide a working directory", call.=FALSE)
}

# Set the working directory
setwd(opt$workdir)
#setwd("/Users/jc428796/OneDrive - James Cook University(1)/PhD/Bioinformatics/Driver_Prioritisation_Comparison_Pipeline")
#setwd("E:/JCU OneDrive/OneDrive - James Cook University/PhD/Bioinformatics/Driver_Prioritisation_Comparison_Pipeline")

dir.create("plots/QC",recursive = T)


#############################
# Get cell line information
###############################

CCLE_sample_info <- read_csv("data/CCLE/Model.csv") %>%
  filter(!is.na(OncotreeLineage)) %>%
  filter(!is.na(ModelID)) %>%
  filter(!is.na(StrippedCellLineName)) %>%
  dplyr::rename(DepMap_ID=ModelID, lineage=OncotreeLineage) %>%
  mutate(lineage=gsub(" |/","_",lineage))

CCLE_model_profile <- read_csv("data/CCLE/OmicsDefaultModelProfiles.csv") %>%
  dplyr::rename(DepMap_ID=ModelID)

# Only retain cell lines with all needed info
CCLE_sample_info <- CCLE_sample_info %>%
  filter(DepMap_ID %in% intersect(CCLE_sample_info$DepMap_ID, CCLE_model_profile$DepMap_ID))
CCLE_model_profile <- CCLE_model_profile %>%
  filter(DepMap_ID %in% intersect(CCLE_sample_info$DepMap_ID, CCLE_model_profile$DepMap_ID))

# Removing unknown lineages, keeping only lineages with n>10, removing normal cells
lineages <- CCLE_sample_info %>%
  filter(!is.na(lineage)) %>%
  group_by(lineage) %>%
  summarise(count = n()) %>%
  filter(count > 10) %>%
  filter(!lineage %in% c("Normal"))

# Filtering data to only retain selected lineages
CCLE_sample_info <- CCLE_sample_info %>%
  filter(lineage %in% lineages$lineage)

# Creating a df of DepMap ID / Cell Name mapping
CCLE_IDs <- CCLE_sample_info %>%
  dplyr::select(DepMap_ID, cell_ID = StrippedCellLineName)

CCLE_model_profile <- CCLE_model_profile %>%
  inner_join(CCLE_IDs, by = c("DepMap_ID"))

#all(CCLE_model_profile$cell_ID %in% CCLE_IDs$cell_ID)
#all(CCLE_IDs$cell_ID %in% CCLE_model_profile$cell_ID)

#########################
# Filtering Expression data to most variable genes
#########################
CCLE_counts <- data.table::fread("data/CCLE/OmicsExpressionGenesExpectedCountProfile.csv") %>%
  dplyr::mutate(across(-1, round)) %>%
  # Replacing DepMap IDs with cell names
  dplyr::rename(ProfileID = 1) %>%
  dplyr::right_join(CCLE_model_profile, by = "ProfileID") %>%
  na.omit() %>%
  dplyr::select(-c(ProfileID,ProfileType,DepMap_ID)) %>% dplyr::relocate(cell_ID) %>%
  # Fixing gene_IDs by removing everything after a space
  data.table::transpose(make.names = "cell_ID", keep.names = "gene_ID") %>%
  as.data.frame() %>%
  mutate(gene_ID = gsub(" .*","",gene_ID)) %>%
  dplyr::mutate(across(-gene_ID, as.numeric))

# Removing duplicate genes
dup <- rle(sort(CCLE_counts$gene_ID))
CCLE_counts <- CCLE_counts %>%
  filter(gene_ID %in% dup$values[dup$lengths==1])

rm(dup)
  
# Getting most variable genes
variance <- apply(CCLE_counts[,seq(2,ncol(CCLE_counts))], 1, var)
names(variance) <- CCLE_counts$gene_ID

CCLE_counts <- CCLE_counts %>% column_to_rownames("gene_ID")
CCLE_counts <- CCLE_counts[head(names(sort(variance,decreasing = T)),3000),]

rm(variance)

###############################
# DESeq2 Normalisation
###############################
# Here the count data is normalised using DESeq2 so that PCA plots can be generated

colData <- CCLE_sample_info %>%
  select(cell_ID = StrippedCellLineName, lineage) %>%
  filter(cell_ID %in% colnames(CCLE_counts)) %>%
  column_to_rownames("cell_ID")

CCLE_counts <- CCLE_counts[,rownames(colData)]

all(rownames(colData) == colnames(CCLE_counts))

dds <- DESeqDataSetFromMatrix(countData = CCLE_counts, 
                              colData = colData,
                              design = ~lineage)

dds <- DESeq(dds)

vst <- varianceStabilizingTransformation(dds, blind = FALSE)

###############################
# PCA Visualisation
###############################

DESeq2::plotPCA(vst, intgroup = "lineage")

ggsave("plots/QC/PCA_all_CCLE.png", width = 4000, height = 4000 ,units = "px")

pca_data <- DESeq2::plotPCA(vst, intgroup = "lineage", returnData = T)

ggplot(pca_data, aes(PC1, PC2)) +
  geom_point(size = 1) +
  facet_wrap(~lineage)

ggsave("plots/QC/PCA_all_CCLE_panel.png", width = 4000, height = 4000 ,units = "px")

###############################
# Manual Cell Selection Based on PCA
###############################
# Here we are selecting only cells that appear to be of a single cell type, rather than forming subclusters
# This is important so that a pseudo-normal cell line can be created based on the averaged expression in all cell lines for that lineage


cell_list <- pca_data %>%
  #Skin remove <0
  filter(case_when(lineage=="Skin"~PC2>0,T~is.numeric(PC2))) %>%
  # remove prostate
  # remove eye
  filter(!lineage %in% c("Prostate", "Eye")) %>%
  select(cell_ID = name, lineage) %>%
  arrange(lineage)

write_csv(cell_list, "data/cell_list.csv")

ggplot(pca_data %>% dplyr::mutate(discard = !name %in% cell_list$cell_ID), aes(x=PC1, y=PC2, colour=discard)) +
  geom_point(size = 1) +
  scale_colour_manual(values = c("black","red")) +
  facet_wrap(~group,)

ggsave("plots/QC/PCA_all_CCLE_panel_discarded.png", width = 4000, height = 4000 ,units = "px")

