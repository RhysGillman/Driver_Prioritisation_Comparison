# prepare_specific_data.r
# This script prepares all input data for running through the driver prioritisation pipeline

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(DESeq2, quietly = T))
suppressPackageStartupMessages (library(ggrepel, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(biomaRt, quietly = T))

# Handling input arguments
option_list = list(
  make_option(c("-w", "--workdir"), type="character", default=NULL, 
              help="working directory", metavar ="Working Directory"),
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network to use", metavar ="Network"),
  make_option(c("-c", "--confidence"), type="double", default=0.4, 
              help="lower bound confidence threshold to keep network edges, value between 0 and 1", metavar ="Confidence")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#if (is.null(opt$workdir)){
#  print_help(opt_parser)
#  stop("Please provide a working directory", call.=FALSE)
#}

# Set the working directory
#setwd(opt$workdir)
#setwd("/Users/jc428796/OneDrive - James Cook University(1)/PhD/Bioinformatics/Driver_Prioritisation_Comparison_Pipeline")
#setwd("E:/JCU OneDrive/OneDrive - James Cook University/PhD/Bioinformatics/Driver_Prioritisation_Comparison_Pipeline")

network_choice <- opt$network
network_confidence <- opt$confidence


#############################
# Get cell list
###############################
# This cell list is generated by the choose_cells.r script

cell_list <- read_csv("data/cell_list.csv")

#############################
# Get cell line information
###############################

CCLE_sample_info <- read_csv("data/CCLE/Model.csv") %>%
  dplyr::rename(DepMap_ID=ModelID, lineage=OncotreeLineage, cell_ID=StrippedCellLineName) %>%
  filter(cell_ID %in% cell_list$cell_ID) %>%
  mutate(lineage=gsub(" |/","_",lineage))

CCLE_model_profile <- read_csv("data/CCLE/OmicsDefaultModelProfiles.csv") %>%
  dplyr::rename(DepMap_ID=ModelID)

# Only retain cell lines with all needed info
CCLE_sample_info <- CCLE_sample_info %>%
  filter(DepMap_ID %in% intersect(CCLE_sample_info$DepMap_ID, CCLE_model_profile$DepMap_ID))
CCLE_model_profile <- CCLE_model_profile %>%
  filter(DepMap_ID %in% intersect(CCLE_sample_info$DepMap_ID, CCLE_model_profile$DepMap_ID))

# Creating a df of DepMap ID / Cell Name mapping
CCLE_IDs <- CCLE_sample_info %>%
  dplyr::select(DepMap_ID, cell_ID)

CCLE_model_profile <- CCLE_model_profile %>%
  inner_join(CCLE_IDs, by = c("DepMap_ID"))

#########################
# Filtering Expression data based on TPM to only include expressed genes
#########################

CCLE_tpm <- data.table::fread("data/CCLE/OmicsExpressionProteinCodingGenesTPMLogp1.csv") %>%
  # Replacing DepMap IDs with cell names
  dplyr::rename(DepMap_ID = 1) %>%
  dplyr::right_join(CCLE_IDs, by = "DepMap_ID") %>%
  na.omit() %>%
  dplyr::select(-DepMap_ID) %>% dplyr::relocate(cell_ID) %>%
  data.table::transpose(make.names = "cell_ID", keep.names = "gene_ID") %>%
  # Fixing gene_IDs by removing everything after a space
  mutate(gene_ID = gsub(" .*","",gene_ID)) %>%
  arrange(gene_ID) %>%
  column_to_rownames("gene_ID") %>%
  # Trimming genes with low reads
  
  # Filtering low expression based on TPM as this accounts for gene length
  # Cells contain approx 200 000 mRNA at any time, therefore tpm of 5 ~ 1 transcript per cell
  # Therefore, tpm of 5 is beging used a cutoff to select genes that are being expressed
  # Data is in log2(tpm+1), therefore filtering for expression abofe log2(5+1)
  dplyr::filter(rowMeans(dplyr::select(., where(is.numeric))) > log2(5+1)) %>%
  as.data.frame()

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
  dplyr::filter(gene_ID %in% rownames(CCLE_tpm))

# Removing duplicate genes
dup <- rle(sort(CCLE_counts$gene_ID))
CCLE_counts <- CCLE_counts %>%
  filter(gene_ID %in% dup$values[dup$lengths==1])

rm(dup)

CCLE_counts <- CCLE_counts %>%
  column_to_rownames("gene_ID")

CCLE_tpm <- CCLE_tpm[rownames(CCLE_counts),]

genes <- rownames(CCLE_counts)

#########################
# Network
#########################

if(network_choice=="STRINGv11"){
  
  #########################
  # STRING Aliases
  #########################
  # Setting up object to match STRING IDs with gene names
  
  STRING11_aliases <- read_tsv("data/STRINGv11/9606.protein.aliases.v11.0.txt", skip = 1, col_names = c("string_ID", "alias", "source")) %>%
    filter(alias %in% genes) %>%
    unique()
  
  # Duplicates exist with spurious matches, ordering sources by their frequency to prioritise the sources that most frequently match with the gene list
  
  sources <- STRING11_aliases %>%
    group_by(source) %>% 
    summarise(count = n()) %>%
    arrange(desc(count))
  
  dups <- STRING11_aliases %>%
    filter(string_ID %in% STRING11_aliases$string_ID[which(duplicated(STRING11_aliases$string_ID))])
  
  # Create empty df
  selected_ids <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(selected_ids) <- c("string_ID", "alias", "source")
  
  # loop through duplicated IDs and take the alias source which most highly reflects the gene list
  for(dup_id in unique(dups$string_ID)){
    temp <- dups %>%
      filter(string_ID == dup_id)
    s <- sources$source[which(sources$source %in% temp$source)[1]]
    temp <- temp %>%
      filter(source == s)
    selected_ids <- rbind(selected_ids, temp)
  }
  
  # Remove duplicated IDs in the original object and then rejoin with only the selected sources
  STRING11_aliases <- STRING11_aliases %>%
    filter(!string_ID %in% STRING11_aliases$string_ID[which(duplicated(STRING11_aliases$string_ID))]) %>%
    rbind(selected_ids)
  
  # The presence of duplicates at this point are due to cases (n = approx 200) where a single source has multiple aliases for a single string_ID, at this point these will be filtered out to prevent spurious connections in network which may bias centrality of nodes.
  STRING11_aliases <- STRING11_aliases %>%
    filter(!string_ID %in% STRING11_aliases$string_ID[which(duplicated(STRING11_aliases$string_ID))])
  
  rm(temp, s, dups, dup_id, sources, selected_ids)
  
  #########################
  # STRING Network
  #########################
  
  network <- read_delim("data/STRINGv11/9606.protein.actions.v11.0.txt") %>%
    # Taking only edges with known interaction direction
    dplyr::filter(is_directional) %>%
    # Taking only edges that are may affect gene expression, directly or indirectly (eg. ptmod, activation, inhibition)
    dplyr::filter(mode %in% c("inhibition", "activation", "expression", "ptmod")) %>%
    # Filtering for only high confidence interactions
    dplyr::filter(score >= network_confidence*1000) %>%
    # Changing protein names to match with rna data
    inner_join(STRING11_aliases, by = c("item_id_a"="string_ID")) %>%
    dplyr::select(-item_id_a) %>% dplyr::rename(protein_1 = alias) %>%
    inner_join(STRING11_aliases, by = c("item_id_b"="string_ID")) %>%
    dplyr::select(-item_id_b) %>% dplyr::rename(protein_2 = alias) %>%
    dplyr::select(protein_1, protein_2, score, a_is_acting) %>%
    group_by(protein_1, protein_2, a_is_acting) %>%
    summarise(score = max(score)) %>%
    na.omit()
  
  # Note, every single edge in the network has a reversed copy, ie. one where "a_is_acting" is TRUE and one where it is false
  # This can be checked below
  #forward_edges <- STRING11_network %>% mutate(edge = paste0(protein_1,"-",protein_2)) %>% pull(edge)
  #reverse_edges <- STRING11_network %>% mutate(edge = paste0(protein_2,"-",protein_1)) %>% pull(edge)
  #reverse_edges[which(reverse_edges %in% forward_edges)] %>% length() / nrow(STRING11_network)
  
  # To get directed edges, simply take "a_is_acting"==TRUE
  network_directed <- network %>%
    dplyr::filter(a_is_acting) %>%
    dplyr::select(protein_1, protein_2, confidence = score) %>%
    dplyr::filter(protein_1 != protein_2)
  
  # To get an undirected version of this network, take the highest confidence interaction from any reciprocal interactions
  network_undirected <- network_directed %>%
    group_by(edge = paste0(pmin(protein_1, protein_2), ":", pmax(protein_1, protein_2))) %>%
    summarise(confidence = max(confidence)) %>%
    ungroup() %>%
    separate(col = edge, into = c("protein_1", "protein_2"), sep = ":") %>%
    unique()
  
  rm(STRING11_aliases,network)
}

if(network_choice=="STRINGv11_undirected_only"){
  
  #########################
  # STRING Aliases
  #########################
  # Setting up object to match STRING IDs with gene names
  
  STRING11_aliases <- read_tsv("data/STRINGv11/9606.protein.aliases.v11.0.txt", skip = 1, col_names = c("string_ID", "alias", "source")) %>%
    filter(alias %in% genes) %>%
    unique()
  
  # Duplicates exist with spurious matches, ordering sources by their frequency to prioritise the sources that most frequently match with the gene list
  
  sources <- STRING11_aliases %>%
    group_by(source) %>% 
    summarise(count = n()) %>%
    arrange(desc(count))
  
  dups <- STRING11_aliases %>%
    filter(string_ID %in% STRING11_aliases$string_ID[which(duplicated(STRING11_aliases$string_ID))])
  
  # Create empty df
  selected_ids <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(selected_ids) <- c("string_ID", "alias", "source")
  
  # loop through duplicated IDs and take the alias source which most highly reflects the gene list
  for(dup_id in unique(dups$string_ID)){
    temp <- dups %>%
      filter(string_ID == dup_id)
    s <- sources$source[which(sources$source %in% temp$source)[1]]
    temp <- temp %>%
      filter(source == s)
    selected_ids <- rbind(selected_ids, temp)
  }
  
  # Remove duplicated IDs in the original object and then rejoin with only the selected sources
  STRING11_aliases <- STRING11_aliases %>%
    filter(!string_ID %in% STRING11_aliases$string_ID[which(duplicated(STRING11_aliases$string_ID))]) %>%
    rbind(selected_ids)
  
  # The presence of duplicates at this point are due to cases (n = approx 200) where a single source has multiple aliases for a single string_ID, at this point these will be filtered out to prevent spurious connections in network which may bias centrality of nodes.
  STRING11_aliases <- STRING11_aliases %>%
    filter(!string_ID %in% STRING11_aliases$string_ID[which(duplicated(STRING11_aliases$string_ID))])
  
  rm(temp, s, dups, dup_id, sources, selected_ids)
  
  #########################
  # STRING Network
  #########################
  
  network <- read_delim("data/STRINGv11/9606.protein.actions.v11.0.txt") %>%
    # Taking only edges that are may affect gene expression, directly or indirectly (eg. ptmod, activation, inhibition)
    dplyr::filter(mode %in% c("inhibition", "activation", "expression", "ptmod")) %>%
    # Filtering for only high confidence interactions
    dplyr::filter(score >= network_confidence*1000) %>%
    # Changing protein names to match with rna data
    inner_join(STRING11_aliases, by = c("item_id_a"="string_ID")) %>%
    dplyr::select(-item_id_a) %>% dplyr::rename(protein_1 = alias) %>%
    inner_join(STRING11_aliases, by = c("item_id_b"="string_ID")) %>%
    dplyr::select(-item_id_b) %>% dplyr::rename(protein_2 = alias) %>%
    dplyr::select(protein_1, protein_2, score) %>%
    group_by(protein_1, protein_2) %>%
    summarise(score = max(score)) %>%
    na.omit()
    
  
  # To get an undirected version of this network, take the highest confidence interaction from any reciprocal interactions
  network_undirected <- network %>%
    dplyr::filter(protein_1 != protein_2) %>%
    group_by(edge = paste0(pmin(protein_1, protein_2), ":", pmax(protein_1, protein_2))) %>%
    summarise(confidence = max(confidence)) %>%
    ungroup() %>%
    separate(col = edge, into = c("protein_1", "protein_2"), sep = ":") %>%
    unique()
  
  rm(STRING11_aliases,network)
}

#########################
# Mutations
#########################

# Reading in mutation data from CCLE
CCLE_mutations <- read_csv("data/CCLE/OmicsSomaticMutations.csv") %>% 
  # Adding CCLE names for cell lines of interest
  inner_join(CCLE_IDs, by = "DepMap_ID") %>%
  # Removing silent mutations
  filter(VariantInfo != "SILENT") %>%
  # Removing cell lines not included in analysis
  filter(!is.na(cell_ID))  %>%
  dplyr::relocate(cell_ID)


CCLE_mutation_matrix <- CCLE_mutations %>%
  dplyr::select(cell_ID, mutated_genes = HugoSymbol) %>%
  # Need to add all genes to include genes that are not mutated in any sample
  right_join((CCLE_counts %>% rownames_to_column("gene_ID") %>% dplyr::select(gene_ID)), by = c("mutated_genes" = "gene_ID")) %>%
  # Make unique so there is only one entry per mutated gene per patient, thus the matrix will be a 0,1 binary matrix
  unique()

print(paste0(
  CCLE_mutation_matrix %>%
    na.omit() %>%
    group_by(cell_ID) %>%
    summarise(count = n()) %>%
    pull(count) %>%
    mean() %>% round(),
  " non-silent mutations per patient on average"
))


# Convert 2 column DF into matrix
CCLE_mutation_matrix <- table(CCLE_mutation_matrix) %>% as.data.frame.matrix() %>% t() %>% as.data.frame()

#rm(CCLE_mutations)



#########################
# Copy Number
#########################
## DepMap CCLE data is in the form log2(CN ratio + 1)

CCLE_copy_number <- data.table::fread("data/CCLE/OmicsCNGene.csv")

CCLE_copy_number <- CCLE_copy_number %>%
  dplyr::rename(DepMap_ID = 1) %>%
  # Acquiring correct cell_IDs
  dplyr::inner_join(CCLE_IDs, by = "DepMap_ID") %>%
  dplyr::select(-DepMap_ID) %>% relocate(cell_ID) %>%
  filter(!is.na(cell_ID)) %>%
  data.table::transpose(make.names = "cell_ID", keep.names = "gene_ID") %>%
  mutate(gene_ID = gsub("\\s.*","", gene_ID)) %>%
  column_to_rownames("gene_ID") %>%
  as.data.frame() %>%
  na.omit()

# Need to correct copy number of X and Y chromosome genes


# Using biomaRt to get info on which chromosome each gene is on

ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
gene_chr_map <- getBM(attributes=c('hgnc_symbol','chromosome_name'), mart = ensembl) %>%
  filter(chromosome_name %in% c(1:23, "X", "Y")) %>%
  unique() %>%
  filter(hgnc_symbol != "", chromosome_name != "") %>%
  filter(hgnc_symbol %in% rownames(CCLE_copy_number)) %>%
  # Note, there are no Y-chromosome genes
  mutate(chromosome_name = ifelse(chromosome_name %in% c("X", "Y"), chromosome_name, "autosome"))

# Checking sex annotation based on X chromosome copy number

check_sex <- CCLE_copy_number %>%
  rownames_to_column("gene_ID") %>%
  pivot_longer(cols = -gene_ID, names_to = "cell_ID", values_to = "cn") %>%
  left_join(CCLE_sample_info %>% dplyr::select(cell_ID, Sex), by = "cell_ID") %>%
  left_join(gene_chr_map, by = c("gene_ID" = "hgnc_symbol")) %>%
  filter(chromosome_name == "X") %>%
  group_by(cell_ID, Sex) %>%
  summarise(mean_cn = mean(2^(cn)-1))

ggplot(check_sex, aes(x=mean_cn, colour = Sex)) +
  geom_density() +
  xlab("Mean X Chromosome Copy Number")

ggsave("plots/QC/CCLE_sex_annotation.png")

check_sex <- check_sex %>%
  mutate(corrected_sex = ifelse(mean_cn > 0.75, "Female", "Male"))


ggplot(check_sex, aes(x=mean_cn, colour = corrected_sex)) +
  geom_density() +
  xlab("Mean X Chromosome Copy Number")

ggsave("plots/QC/CCLE_sex_annotation_corrected.png")

# Changing sex annotations in sample_info

CCLE_sample_info <- CCLE_sample_info %>%
  left_join(check_sex[c("cell_ID","corrected_sex")] %>% unique() %>% na.omit(), by = "cell_ID")

rm(check_sex)

# Check Chromosome Expression

chr_expression <- CCLE_copy_number %>%
  rownames_to_column("gene_ID") %>%
  pivot_longer(cols = -gene_ID, names_to = "cell_ID", values_to = "cn") %>%
  left_join(CCLE_sample_info %>% dplyr::select(cell_ID, corrected_sex), by = "cell_ID") %>%
  left_join(gene_chr_map, by = c("gene_ID" = "hgnc_symbol")) %>%
  mutate(chromosome_name = ifelse(is.na(chromosome_name), "not_annotated", chromosome_name)) %>%
  group_by(cell_ID, corrected_sex, chromosome_name) %>%
  summarise(mean_cn = mean(2^(cn)-1))


ggplot(chr_expression, aes(x=mean_cn, colour = corrected_sex)) +
  geom_density() +
  xlab("Mean X Chromosome Copy Number") +
  facet_wrap(~chromosome_name)

ggsave("plots/QC/CCLE_CNV_not_corrected.png")

# Based on the generated plots, it seems that a correction has not been applied to X chromosome gene copy number based on sex
# Adding correction for this

x_genes <- gene_chr_map %>% filter(chromosome_name == "X") %>% pull(hgnc_symbol)
male_cells <- CCLE_sample_info %>% dplyr::select(cell_ID, corrected_sex) %>% filter(corrected_sex == "Male") %>% pull(cell_ID) %>% unique()

cnv_x <- CCLE_copy_number[x_genes,male_cells]
# Undo the log2(CN ratio + 1)
cnv_x <- 2^cnv_x -1
# Correct X chromosome genes by adding 0.5 ratio
cnv_x <- cnv_x + 0.5
# Redo normalisation step
cnv_x <- log2(cnv_x + 1)

CCLE_copy_number_corrected <- CCLE_copy_number
CCLE_copy_number_corrected[x_genes, male_cells] <- cnv_x

rm(x_genes, male_cells, cnv_x)

# Rechecking

chr_expression <- CCLE_copy_number_corrected %>%
  rownames_to_column("gene_ID") %>%
  pivot_longer(cols = -gene_ID, names_to = "cell_ID", values_to = "cn") %>%
  left_join(CCLE_sample_info %>% dplyr::select(cell_ID, corrected_sex), by = "cell_ID") %>%
  left_join(gene_chr_map, by = c("gene_ID" = "hgnc_symbol")) %>%
  mutate(chromosome_name = ifelse(is.na(chromosome_name), "not_annotated", chromosome_name)) %>%
  group_by(cell_ID, corrected_sex, chromosome_name) %>%
  summarise(mean_cn = mean(2^(cn)-1))


ggplot(chr_expression, aes(x=mean_cn, colour = corrected_sex)) +
  geom_density() +
  xlab("Mean X Chromosome Copy Number") +
  facet_wrap(~chromosome_name)

ggsave("plots/QC/CCLE_CNV_corrected.png")

rm(chr_expression, CCLE_copy_number, gene_chr_map, ensembl)

# Copy number is being converted to binary gain or loss values
## DepMap CCLE data is in the form log2(CN ratio + 1). 
## ie. a single deletion = log2(1/2 + 1) = 0.585, 
## homozygous deletion = 0. 
## Normal copy number = log2(2/2+1) = 1, 
## extra copy = log2(3/2+1) = 1.32, then 1.585, etc.
## Using following thresholds:
### log2(1.5) = 0.585 for deletions
### log2(3) = 1.585 for amplifications

CCLE_copy_number_corrected[CCLE_copy_number_corrected<log2(1.5)] <- -1
CCLE_copy_number_corrected[CCLE_copy_number_corrected>log2(1.5) & CCLE_copy_number_corrected<log2(3)] <- 0
CCLE_copy_number_corrected[CCLE_copy_number_corrected>log2(3)] <- 1
CCLE_copy_number_corrected <- round(CCLE_copy_number_corrected,0)


#########################
# Gold Standards
#########################

gold_standards <- read_csv("validation_data/all_gold_standards.csv")


#########################
# Filtering
#########################

# Because a diverse range of data types are being used for this analysis, it is necessary to make sure that 
# the same genes and samples are being used for each data type, and therefore only genes and cells for which 
# all of the data is available will be kept.

## Common Cells

cell_list <- Reduce(intersect, list(
  #RNA data
  colnames(CCLE_counts),
  colnames(CCLE_tpm),
  #Mutation data
  colnames(CCLE_mutation_matrix),
  #CNV data
  colnames(CCLE_copy_number_corrected)
))

# Additionally, need to have gold standard data for at least some cell lines in each lineage, though don't necessarily need gold standards for all cell lines
# This is because the "cohort" of cells is still useful for the algorithms when predicting drivers, and the cell lines lacking gold standard data can be
# Removed at a later time

CCLE_sample_info <- CCLE_sample_info %>%
  filter(lineage %in% gold_standards$lineage)

cell_list <- intersect(cell_list, CCLE_sample_info$cell_ID)




CCLE_counts <- CCLE_counts[,cell_list]
CCLE_tpm <- CCLE_tpm[,cell_list]
CCLE_mutation_matrix <- CCLE_mutation_matrix[,cell_list]
CCLE_mutations <- CCLE_mutations %>% filter(cell_ID %in% cell_list)
CCLE_copy_number_corrected <- CCLE_copy_number_corrected[,cell_list]
CCLE_sample_info <- CCLE_sample_info %>% filter(cell_ID %in% cell_list)
CCLE_IDs <- CCLE_IDs %>% filter(cell_ID %in% cell_list)

## Common Genes

gene_list <- Reduce(intersect, list(
  #RNA data
  rownames(CCLE_counts),
  rownames(CCLE_tpm),
  #Mutation data
  rownames(CCLE_mutation_matrix),
  #CNV data
  rownames(CCLE_copy_number_corrected),
  #STRING network
  #unique(append(network_directed$protein_1, network_directed$protein_2)),
  unique(append(network_undirected$protein_1, network_undirected$protein_2))
))
if(exists("network_directed")){
  # Need to ensure that BOTH ends of each interaction are in the gene list
  network_directed <- network_directed %>% 
    filter(protein_1 %in% gene_list & protein_2 %in% gene_list)
}

network_undirected <- network_undirected %>%
  filter(protein_1 %in% gene_list & protein_2 %in% gene_list)

gene_list <- unique(append(network_undirected$protein_1, network_undirected$protein_2))

# Now filter all other data to only include these genes

CCLE_counts <- CCLE_counts[gene_list,]
CCLE_tpm <- CCLE_tpm[gene_list,]
CCLE_mutation_matrix <- CCLE_mutation_matrix[gene_list,]
CCLE_mutations <- CCLE_mutations %>% filter(HugoSymbol %in% gene_list)
CCLE_copy_number_corrected <- CCLE_copy_number_corrected[gene_list,]
gold_standards <- gold_standards %>% filter(gene_ID %in% gene_list)


#########################
# Saving Data
#########################
DATA_DIR <- paste0("validation_data/CCLE_", network_choice)

dir.create(DATA_DIR)

write_csv(CCLE_counts %>% rownames_to_column("gene_ID"), paste0(DATA_DIR,"/counts.csv"))
write_csv(CCLE_tpm %>% rownames_to_column("gene_ID"), paste0(DATA_DIR,"/tpm.csv"))
write_csv(CCLE_mutation_matrix %>% rownames_to_column("gene_ID"), paste0(DATA_DIR,"/mutations.csv"))
write_csv(CCLE_mutations, paste0(DATA_DIR,"/mutations_MAF.csv"))
write_csv(CCLE_copy_number_corrected %>% rownames_to_column("gene_ID"), paste0(DATA_DIR,"/cnv.csv"))
write_csv(CCLE_sample_info, paste0(DATA_DIR,"/sample_info.csv"))
write_csv(network_directed, paste0(DATA_DIR,"/network_directed.csv"))
write_csv(network_undirected, paste0(DATA_DIR,"/network_undirected.csv"))
write_csv(gold_standards, paste0(DATA_DIR,"/gold_standards.csv"))
write.table(CCLE_sample_info %>% dplyr::select(lineage) %>% arrange(lineage) %>% unique, 
            sep="\t", 
            col.names=FALSE, 
            row.names=FALSE, 
            quote = FALSE,
            eol = "\n",
            file = paste0(DATA_DIR,"/lineages.txt"))

#########################
# Summaries
#########################

SUMMARY_DIR <- paste0(DATA_DIR, "/summary_info")

if(!dir.exists(SUMMARY_DIR)){
  dir.create(SUMMARY_DIR)
}


# Cell Counts
write_csv(
  CCLE_sample_info %>% group_by(lineage) %>% summarise(cell_count = n()),
  paste0(SUMMARY_DIR,"/cell_counts.csv")
)

# Mutations

mutation_summary <- CCLE_mutation_matrix %>%
  rownames_to_column("gene_ID") %>%
  pivot_longer(cols = -gene_ID, names_to = "cell_ID", values_to = "mutated") %>%
  group_by(cell_ID, mutated) %>%
  summarise(count = n()) %>%
  left_join(CCLE_sample_info[c("cell_ID","lineage")], by = "cell_ID") %>%
  group_by(lineage, mutated) %>%
  summarise(mean_count = mean(count), max_count = max(count), min_count = min(count))

write_csv(mutation_summary, paste0(SUMMARY_DIR,"/mutation_summary.csv"))
rm(mutation_summary)

#CNV

cnv_summary <- CCLE_copy_number_corrected %>% 
  rownames_to_column("gene_ID") %>%
  pivot_longer(cols = -gene_ID, names_to = "cell_ID", values_to = "copy_number") %>%
  group_by(cell_ID, copy_number) %>%
  summarise(count = n()) %>%
  left_join(CCLE_sample_info[c("cell_ID","lineage")], by = "cell_ID") %>%
  group_by(lineage, copy_number) %>%
  summarise(mean_count = mean(count), max_count = max(count), min_count = min(count))

write_csv(cnv_summary, paste0(SUMMARY_DIR,"/CNV_summary.csv"))

rm(cnv_summary)

# Network Summary

analyse_network <- function(network, name){
  # Layout
  cols <- ncol(network)
  
  if(cols <= 3){
    form <- "long"
  }else{
    form <- "matrix"
  }
  
  if(form == "long"){
    
    
    if(cols==2){
      weighted <- FALSE
      weight_min <- NA
      weight_max <- NA
    }else{
      if(is_empty(select_if(network, is.numeric))){
        weighted <- FALSE
      }else{
        weighted <- TRUE
        weight_min <- min(pull(dplyr::select_if(network,is.numeric)))
        weight_max <- max(pull(dplyr::select_if(network,is.numeric)))
      }
      
    }
    
    sub_network <- network %>% select_if(is.character)
    
    if(ncol(sub_network) != 2){
      print("Error: There are not 2 character columns. Check col_types")
      return(NULL)
    }
    
    colnames(sub_network) <- c("gene1", "gene2")
    
    if(nrow(network) == nrow(unique(sub_network))){
      unique_rows <- TRUE
    }else(
      unique_rows <- FALSE
    )
    
    n_edges <- nrow(unique(sub_network))
    
    if(identical(sort(unique(sub_network$gene1)), sort(unique(sub_network$gene2)))){
      equal_genes <- TRUE
    }else{
      equal_genes <- FALSE
    }
    
    forward_edges <- sub_network %>% mutate(edge = paste0(gene1,"-",gene2)) %>% pull(edge)
    reverse_edges <- sub_network %>% mutate(edge = paste0(gene2,"-",gene1)) %>% pull(edge)
    
    if(identical(sort(unique(forward_edges)), sort(unique(reverse_edges)))){
      equal_edges <- TRUE
    }else{
      equal_edges <- FALSE
    }
    n_uniq_gene1 <- length(unique(sub_network$gene1))
    n_uniq_gene2 <- length(unique(sub_network$gene2))
    n_uniq_genes_total <- length(unique(append(sub_network$gene1, sub_network$gene2)))
    reversed_edges <- paste0(round(length(intersect(forward_edges, reverse_edges))/length(forward_edges)*100,2),"%")
  }
  
  df_summary <- data.frame(name, form, weighted, weight_min, weight_max, unique_rows, equal_genes, equal_edges, n_uniq_genes_total, n_uniq_gene1, n_uniq_gene2, n_edges,reversed_edges)
  return(df_summary)
  
}

network_summary <- analyse_network(network_undirected, paste0(network_choice,"_undirected_network"))

if(exists("network_directed")){
  network_summary <- rbind(network_summary,analyse_network(network_directed, paste0(network_choice,"_directed_network")))
}

write_csv(network_summary, paste0(SUMMARY_DIR,"/network_summary.csv"))

rm(network_summary)

# Gold Standards Summary

gold_standards_summary <- gold_standards %>%
  group_by(lineage,cell_ID) %>%
  summarise(total = n(), drug_sensitive = sum(sensitive, na.rm = T), gene_dependent = sum(dependent, na.rm = T)) %>%
  group_by(lineage) %>%
  summarise(mean_total = mean(total),
            mean_drug_sensitive = mean(drug_sensitive),
            mean_gene_dependent = mean(gene_dependent)
            )

write_csv(gold_standards_summary, paste0(SUMMARY_DIR,"/gold_standards_summary.csv"))

rm(gold_standards_summary)



