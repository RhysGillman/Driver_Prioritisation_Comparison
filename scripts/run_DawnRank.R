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

for(s in list.files("scripts/DawnRank/")){
  source(paste0("scripts/DawnRank/",s))
}



#############################
# Sample Info
#############################

sample_info <- read_csv(paste0("validation_data/CCLE_",network_choice,"/sample_info.csv")) %>% filter(lineage==cell_type)
samples <- sample_info$cell_ID %>% sort()

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

# Need to convert 2-column dataframe of interactions into a complete matrix
network <-  fread(paste0("validation_data/CCLE_",network_choice,"/network_directed.csv")) %>%
  mutate(confidence = ifelse(confidence>0, 1, 0)) %>%
  pivot_wider(names_from = 2, values_from = 3) %>%
  complete(protein_1 = names(.)[-1]) %>%
  data.table::transpose(make.names = "protein_1", keep.names = "protein_2") %>%
  complete(protein_2 = names(.)[-1]) %>%
  data.table::transpose(make.names = "protein_2", keep.names = "protein_1") %>%
  column_to_rownames("protein_1") %>%
  mutate_all(~replace_na(., 0)) %>%
  as.matrix()

genes <- gene_list <- Reduce(intersect, list(
  rownames(rna),
  rownames(mutation),
  rownames(network),
  colnames(network)
))

# Make sure that all columns and rows display data in the same order
rna <- rna[genes, samples]
mutation <- mutation[genes, samples]
network <- network[genes,genes]

# Testing
#load("../DawnRank/data/brcaExamplePathway.rda")
#DawnRank_network <- brcaExamplePathway
#rm(brcaExamplePathway)
#genes_test <- intersect(genes, rownames(DawnRank_network))
#DawnRank_network <- which(DawnRank_network[genes_test,genes_test]==1)
#network_test <- which(network[genes_test,genes_test]==1)
#overlap <- intersect(network_test, DawnRank_network)
#length(overlap)/length(network_test)


#############################
# Running DawnRank
#############################

damping <- dawnDamping(network, 3)
dawnMatrix <- DawnMatrix(network) 

res_df <- data.frame(Gene = character(), Patient = character(), Rank = numeric(), PercentRank = numeric())

all(rownames(dawnMatrix) == rownames(mutation))
all(rownames(dawnMatrix) == rownames(rna))
all(rownames(rna) == rownames(mutation))
all(colnames(rna) == colnames(mutation))

for(sample in samples){
  print(sample)
  tumour_rna <- rna[,sample] %>% as.matrix()
  colnames(tumour_rna) <- sample
  normal_rna <- rna[,samples[samples != sample]]
  
  normalizedDawn <- DawnNormalize(tumorMat = tumour_rna, normalMat = normal_rna)
  colnames(normalizedDawn) <- sample
  
  dawn <- Dawn(dawnMatrix, normalizedDawn[,sample], mutation[,sample], damping, maxit = 100, patientTag = sample, epsilon = 1e-04)
  
  #mutated_genes <- dawn$summaryOutput
  
  indiv_res <- dawn$mutatedRanks %>% arrange(desc(PercentRank))
  
  res_df <- rbind(res_df, indiv_res)
}

res_df <- res_df %>%
  dplyr::select(cell_ID = Patient, gene_ID = Gene, PercentRank) %>%
  # Converting PercentRank into ordinal gene ranks
  arrange(desc(PercentRank)) %>%
  group_by(cell_ID) %>%
  mutate(rank = row_number()) %>%
  dplyr::select(-PercentRank) %>%
  arrange(cell_ID)

#############################
# Save Result
#############################


write_csv(res_df, paste0("results/CCLE_",network_choice,"/DawnRank_",cell_type,".csv"))






