#prepare_combined_de_novo_methods_data.R
# This code prepares the input data to run all of the de novo network construction and control methods and stores it in tmp/

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
# Requires paired tumour and normal expression data
# Because this is not possible with the cell line data, we will make pseudonormal samples based on the expression of all 'other' cells iteratively



rna <- fread(paste0("validation_data/CCLE_",network_choice,"/tpm.csv"), select = c("gene_ID",samples)) %>% 
  arrange(gene_ID) %>%
  column_to_rownames("gene_ID") 

suppressWarnings(rm(pseudonormal_rna))
for(sample in samples){
  # Collect expression data for all samples other than the current sample
  other_rna <- rna[,samples[samples != sample]]
  # Collect mean expression for other samples as the pseudo-normal expression for each cell
  tmp_pseudonormal_rna <- enframe(rowMeans(other_rna), name = "gene_ID", value = sample) %>% column_to_rownames("gene_ID")
  
  
  if(!exists("pseudonormal_rna")){
    pseudonormal_rna <- tmp_pseudonormal_rna
  }else{
    pseudonormal_rna <- cbind(pseudonormal_rna, tmp_pseudonormal_rna)
  }
  
}

suppressWarnings(rm(tmp_pseudonormal_rna, other_rna))

# Formatting

genes <- rownames(rna) %>% unique() %>% sort()

rna <- rna[genes, samples]
pseudonormal_rna <- pseudonormal_rna[genes, samples]


#############################
# Create Temporary Files
#############################

write.table(rna, "tmp/tmp_combined_de_novo_methods_tumour_expression.txt", quote = F, col.names = NA, sep = "\t")
write.table(pseudonormal_rna, "tmp/tmp_combined_de_novo_methods_pseudonormal_expression.txt", quote = F, col.names = NA, sep = "\t")

