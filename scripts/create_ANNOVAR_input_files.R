#create_ANNOVAR_input_files.R
# This code converts the CCLE variant calls into the input format required by ANNOVAR

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


if(!dir.exists(paste0("validation_data/CCLE_", network_choice, "/ANNOVAR_input/",cell_type))){
  dir.create(paste0("validation_data/CCLE_", network_choice, "/ANNOVAR_input/",cell_type), recursive = T)
}

#############################
# Sample Info
#############################

sample_info <- read_csv(paste0("validation_data/CCLE_",network_choice,"/sample_info.csv")) %>% filter(lineage==cell_type)
samples <- sample_info$cell_ID %>% sort()

#############################
# Create ANNOVAR Input Files
#############################

# Annovar format Chromosome ("chr" prefix is optional), Start, End, Reference Allele, Alternative Allele
# tab separated
# https://annovar.openbioinformatics.org/en/latest/user-guide/input/

mutations <- fread(paste0("validation_data/CCLE_",network_choice,"/mutations_MAF.csv"), sep = ",")

for(sample in samples){
  annovar_input <- mutations %>%
    filter(cell_ID == sample) %>%
    mutate(End = ifelse(nchar(Ref) > 2 & Alt=="-", Pos + nchar(Ref) - 1, Pos)) %>%
    dplyr::select(Chrom,Start=Pos,End,Ref,Alt)
  write_tsv(annovar_input, col_names = F, paste0("validation_data/CCLE_", network_choice, "/ANNOVAR_input/",cell_type,"/",sample,".avinput"))
  
}
