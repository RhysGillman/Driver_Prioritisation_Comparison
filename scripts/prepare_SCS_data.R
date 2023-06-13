#prepare_SCS_data.R
# This code prepares the input data to run SCS and stores it in tmp/

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
# Data needs to be in the form of log2fold-changes. This is the same as OncoImpact, and therefore the code has been adapted 
# from the OnocoImpact code to achieve the same purpose

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


# Formatting

genes <- Reduce(intersect, list(
  rownames(l2fc_rna),
  rownames(mutation)
))

l2fc_rna <- l2fc_rna[genes, samples]
mutation <- mutation[genes, samples]



#############################
# Create Temporary Files
#############################

write_tsv(l2fc_rna %>% rownames_to_column("Genes"), "tmp/tmp_SCS_EXP.txt")
write_tsv(mutation %>% rownames_to_column("Genes"), "tmp/tmp_SCS_SNP.txt")
write_tsv(cnv %>% rownames_to_column("Genes"), "tmp/tmp_SCS_cnv.txt")



if(!dir.exists(paste0("results/CCLE_",network_choice,"/SCS"))){
  dir.create(paste0("results/CCLE_",network_choice,"/SCS"))
}


