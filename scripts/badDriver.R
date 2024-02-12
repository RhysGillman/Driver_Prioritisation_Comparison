# badDriver.R

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(foreach, quietly = T))
suppressPackageStartupMessages (library(doParallel, quietly = T))


# Handling input arguments
option_list = list(
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network to use", metavar ="Network"),
  make_option(c("-c", "--celltype"), type="character", default="Liver", 
              help="cell type to analyse", metavar ="Cell Type"),
  make_option(c("-t", "--threads"), type="integer", default=4, 
              help="Number of threads to use if running in parallel", metavar ="Rank Aggregation Method"),
  make_option(c("-r", "--refdrivers"), type="character", default="none", 
              help="reference driver set to prioritise. Leave as NULL (default) to not prioritise a reference set", metavar ="Reference Driver Set"),
  make_option(c("-s", "--simnumber"), type="integer", default=10, 
              help="the number of simulations to run", metavar ="Simulation Number"),
  make_option(c("-m", "--mode"), type="character", default="both", 
              help="Mode (options: driver, drug, both (Default))", metavar ="Mode")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

network_choice <- opt$network
cell_type <- opt$celltype
threads <- opt$threads
use_ref_drivers <- opt$refdrivers
n_sims <- opt$simnumber
mode <- opt$mode

if(threads>1){
  #registerDoParallel(cores=threads)
  cl <- makeCluster(threads, outfile = "log/evaluate_SL_partners.log")
  registerDoParallel(cl)
}

#mode="drug"
#cell_type="all"



if(cell_type=="all"){
  cell_types <- read_csv(paste0("validation_data/CCLE_",network_choice,"/sample_info.csv")) %>% pull(lineage) %>% unique()
}else{
  cell_types <- cell_type
}




for(cell_type in cell_types){
  



#############################
# Sample Info
#############################

sample_info <- read_csv(paste0("validation_data/CCLE_",network_choice,"/sample_info.csv")) %>% filter(lineage==cell_type)
samples <- sample_info$cell_ID %>% sort()


#############################
# Create Directories
#############################

if(!dir.exists(paste0("bad_driver_simulations/CCLE_",network_choice,"/",cell_type))){
  dir.create(paste0("bad_driver_simulations/CCLE_",network_choice,"/",cell_type), recursive = T)
}

if(!dir.exists(paste0("bad_driver_simulations/CCLE_",network_choice,"_drugs/",cell_type))){
  dir.create(paste0("bad_driver_simulations/CCLE_",network_choice,"_drugs/",cell_type), recursive = T)
}




if(mode %in% c("both","driver")){
  

#############################
# Mutation Data
#############################

# Mutations

mutation <- fread(paste0("validation_data/CCLE_",network_choice,"/mutations.csv"), select = c("gene_ID",samples)) %>% 
  arrange(gene_ID) %>%
  column_to_rownames("gene_ID")

cnv <- fread(paste0("validation_data/CCLE_",network_choice,"/cnv.csv"), select = c("gene_ID",samples)) %>%
  column_to_rownames("gene_ID")

mutation <- mutation[sort(rownames(mutation)),sort(colnames(mutation))]
cnv <- cnv[sort(rownames(cnv)),sort(colnames(cnv))]


if(!all(colnames(mutation)==colnames(cnv)) & all(rownames(mutation) == rownames(cnv))){
  stop("Error: colnames/rownames of cnv and mutation files do not match")
}


altered_genes <- cnv != 0 | mutation != 0
altered_genes[] <- as.integer(altered_genes)




#############################
# Reference Driver Set
#############################
#use_ref_drivers <- "CGC"

if(!use_ref_drivers=="none"){
  if(use_ref_drivers == "CGC"){
    ref_drivers <- read_csv("data/cancer_reference_genes/CGC/cancer_gene_census.csv") %>%
      dplyr::select(gene_ID = `Gene Symbol`) %>%
      unique() %>%
      pull(gene_ID)
  }
}




#############################
# badDriver prediction
#############################


for(i in 1:n_sims){
  
  bad_drivers <- foreach(sample=samples, .combine = "rbind", .packages = c("tidyverse","foreach")) %dopar% {
    
    if(use_ref_drivers=="none"){
    sample_altered_genes <- rownames(altered_genes)[which(altered_genes[,sample]==1)] %>% 
      # randomise order
      sample(replace = F)
    } else {
    sample_altered_genes <- rownames(altered_genes)[which(altered_genes[,sample]==1)]
    sample_ref_genes <- sample_altered_genes[which(sample_altered_genes %in% ref_drivers)] %>% 
      # randomise order
      sample(replace = F)
    sample_non_ref_genes <- sample_altered_genes[which(!sample_altered_genes %in% ref_drivers)] %>% 
      # randomise order
      sample(replace = F)
    sample_altered_genes <- append(sample_ref_genes,sample_non_ref_genes)
    }
    
    data.frame(lineage = cell_type,cell_ID = sample, driver = sample_altered_genes) %>% mutate(rank = row_number())
    
    
  }
  
  if(!use_ref_drivers=="none"){
  write_csv(bad_drivers, paste0("bad_driver_simulations/CCLE_",network_choice,"/",cell_type,"/bad_sim_no_ref_",i,".csv"))
  } else {
  write_csv(bad_drivers, paste0("bad_driver_simulations/CCLE_",network_choice,"/",cell_type,"/bad_sim_",use_ref_drivers,"_ref_",i,".csv"))
  }
  
}

}


if(mode %in% c("both","drug")){
  


#############################
# Drug Sensitivity Prediction
#############################

#############################
## Get Available Drugs
#############################


all_drugs <- read_csv("validation_data/drug_sensitivity.csv") %>%
  dplyr::select(cell_ID,drug_ID) %>%
  unique() %>%
  filter(cell_ID %in% samples) %>%
  group_by(cell_ID) %>%
  summarise(drugs = list(drug_ID)) %>%
  deframe()

#############################
## BadDrug Prediction
#############################



for(i in 1:n_sims){
  
  bad_drugs <- foreach(sample=samples[samples %in% names(all_drugs)], .combine = "rbind", .packages = c("tidyverse","foreach")) %dopar% {
    
    sample_drugs <- all_drugs[[sample]] %>%
      sample(replace=F)
    
    data.frame(lineage = cell_type,cell_ID = sample, drug_ID = sample_drugs) %>% mutate(rank = row_number())
    
    
  }
  
  
  write_csv(bad_drugs, paste0("bad_driver_simulations/CCLE_",network_choice,"_drugs/",cell_type,"/bad_drug_",i,".csv"))
  
}


}

}