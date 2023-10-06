#prepare_OncoImpact_data.R
# This code prepares the input data to run OncoImpact and stores it in tmp/

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))


# Handling input arguments
option_list = list(
  make_option(c("-w", "--workdir"), type="character", default=NULL, 
              help="working directory", metavar ="Working Directory"),
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network to use", metavar ="Network"),
  make_option(c("-c", "--celltype"), type="character", default="Liver", 
              help="cell type to analyse", metavar ="Cell Type"),
  make_option(c("-t", "--threads"), type="character", default=4, 
              help="Number of threads to use", metavar ="Threads")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

network_choice <- opt$network
cell_type <- opt$celltype
WORK_DIR <- opt$work
threads <- opt$threads


#############################
# Sample Info
#############################

sample_info <- read_csv(paste0("validation_data/CCLE_",network_choice,"/sample_info.csv")) %>% filter(lineage==cell_type)
samples <- sample_info$cell_ID %>% sort()

#############################
# Prepare Input Data
#############################

# RNA

rna <- fread(paste0("validation_data/CCLE_",network_choice,"/tpm.csv"), select = c("gene_ID",samples)) %>% 
  arrange(gene_ID) %>%
  column_to_rownames("gene_ID") 

#RNA data needs to be converted to log2 fold change data
#To do this, compare each sample vs the mean of all others
for(sample in samples){
  tumour_rna <- rna[,sample] %>% as.matrix()
  normal_rna <- rna[,samples[samples != sample]]
  mean_rna <- apply(normal_rna,1,mean)
  diff_rna <- tumour_rna - mean_rna %>% as.data.frame()
  colnames(diff_rna) <- sample
  if(sample == samples[1]){
    l2fc_rna <- diff_rna 
  }else{
    l2fc_rna <- cbind(l2fc_rna, diff_rna)
  }
}
rm(rna, tumour_rna, normal_rna, diff_rna, mean_rna)

# Mutations

mutation <- fread(paste0("validation_data/CCLE_",network_choice,"/mutations.csv"), select = c("gene_ID",samples)) %>% 
  arrange(gene_ID) %>%
  column_to_rownames("gene_ID")

# CNV

cnv <- fread(paste0("validation_data/CCLE_",network_choice,"/cnv.csv"), select = c("gene_ID",samples)) %>%
  column_to_rownames("gene_ID")

# Network

if(network_choice=="own"){
  network <- fread("data/own_networks/OncoImpact.txt", select = c("Gene1","Gene2"))
}else{
  network <- fread(paste0("validation_data/CCLE_",network_choice,"/network_undirected.csv")) %>%
    dplyr::select(1,2)
}



# Formatting

genes <- Reduce(intersect, list(
  rownames(l2fc_rna),
  rownames(mutation),
  append(network %>% pull(1), network %>% pull(2))
)) %>% sort

l2fc_rna <- l2fc_rna[genes, samples]
mutation <- mutation[genes, samples]


#############################
# Create Temporary Files
#############################

write_tsv(l2fc_rna %>% rownames_to_column("Genes"), "tmp/tmp_OncoImpact_EXP.txt")
write_tsv(mutation %>% rownames_to_column("Genes"), "tmp/tmp_OncoImpact_SNP.txt")
write_tsv(network, "tmp/tmp_OncoImpact_network.txt" ,col_names = FALSE)
write_tsv(cnv %>% rownames_to_column("Genes"), "tmp/tmp_OncoImpact_cnv.txt")



if(!dir.exists(paste0("results/CCLE_",network_choice,"/OncoImpact"))){
  dir.create(paste0("results/CCLE_",network_choice,"/OncoImpact"))
}



# Create Config File

#perl_dir <- WORK_DIR
  
  
perl_dir <- gsub("([A-Za-z]):", "/\\L\\1",WORK_DIR, perl = T)
  

cfg_lines <- c(
  paste0("outDir=",perl_dir,"/results/CCLE_",network_choice,"/OncoImpact/",cell_type),
  paste0("scriptDir=",perl_dir,"/scripts/OncoImpact"),
  paste0("numThreads=",threads),
  paste0("cnv=",perl_dir,"/tmp/tmp_OncoImpact_cnv.txt"),
  paste0("exp=",perl_dir,"/tmp/tmp_OncoImpact_EXP.txt"),
  paste0("snp=",perl_dir,"/tmp/tmp_OncoImpact_SNP.txt"),
  "dataType=RNA_SEQ",
  "testMode=0",
  paste0("network=",perl_dir,"/tmp/tmp_OncoImpact_network.txt")
  )

#cfg_lines_escaped <- gsub("=","_escaped=", cfg_lines, fixed = T)
#cfg_lines_escaped <- gsub(" ","\\ ", cfg_lines_escaped, fixed = T)
#cfg_lines_escaped <- gsub("(","\\(", cfg_lines_escaped, fixed = T)  
#cfg_lines_escaped <- gsub(")","\\)", cfg_lines_escaped, fixed = T)

write_lines(cfg_lines, "tmp/tmp_OncoImpact_config.cfg")
