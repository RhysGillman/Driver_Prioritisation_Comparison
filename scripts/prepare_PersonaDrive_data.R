#prepare_PersonaDrive_data.R
# This data prepares the input data to run PersonaDrive and stores it in tmp/

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
# Functions
#############################

# Using code from https://github.com/shahcompbio/drivernet/blob/master/R/getPatientOutlierMatrix.R as the authors of PersonaDrive said they did.

# Personadrive uses a threshold of 0.5 instead of 2

getPatientOutlierMatrix <- function(patExpMatrix, th=2)
  {
  # Get sd of each column
  expSd   <- apply(patExpMatrix, 2, sd)
  # Only keep columns with sd > 0
  id <- expSd > 0
  expSd <- expSd[id]
  patExpMatrix <- patExpMatrix[, id]
  # Get the mean of each column
  expMean <- apply(patExpMatrix, 2, mean)
	
  # dnorma returns the probability density function for a normal distribution given the mean, sd
  # I believe this is the y axis value on a normal distribution
  # Overall, this part converts the values to their probability density value given the sd and mean of the column
  # Thus, higher value means closer to mean
  num <- dnorm(x=t(patExpMatrix), mean=expMean, sd=expSd, log=T)
  num <- t(num)
  # Note, adding a vector to a matrix adds each value to the entire row
  # Thus, if columns are genes, getting the mean + 0.5sd of each gene and getting the threshold probability density value corresponding to this
  numSd <- dnorm( x=t(expMean+th*expSd), mean=expMean, sd=expSd, log=T )
  
  # Creating a matrix that repeats the threshold values per gene for every sample
  y <- rep(numSd, each=dim(patExpMatrix)[1])
  y <- matrix(y, nrow=dim(patExpMatrix)[1], ncol=dim(patExpMatrix)[2])
	
  # TRUE values are outliers because they have probability densitites lower than the threshold
  patOutMatrix <- num <= y
  
  id <- colSums(patOutMatrix) > 0
  patOutMatrix <- patOutMatrix[, id]
  
  return(patOutMatrix)
}


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
  data.table::transpose(make.names = "gene_ID", keep.names = "cell_ID") %>%
  column_to_rownames("cell_ID")

rna_outliers <- getPatientOutlierMatrix(rna, th = 0.5)

# Mutations

mutation <- fread(paste0("validation_data/CCLE_",network_choice,"/mutations.csv"), select = c("gene_ID",samples)) %>% 
  arrange(gene_ID) %>%
  column_to_rownames("gene_ID")

# Network
#The default network used by PersonaDrive is weighted and seems to be non-directed, though the authors don't state this.
#It is in the form on a long table (3 columns protein 1, protein 2, weight.) It is a weighted network with weight values 
#ranging from 0.701 to 0.999. Every row in the table is unique, however every single edge is duplicated in reverse 
#(gene1-gene2 is also gene2-gene1), and thus every single gene is present in both columns of the table.

network <- fread(paste0("validation_data/CCLE_",network_choice,"/network_undirected.csv"))
colnames(network) <- c("protein_1", "protein_2", "confidence")
reverse_edges <- network %>% dplyr::select(protein_1 = protein_2, protein_2 = protein_1, confidence)
network <- rbind(network, reverse_edges) %>% unique()
rm(reverse_edges)
network <- network %>% mutate(confidence = confidence/1000)
colnames(network) <- c("genes1","genes2","Score")



# Formatting

genes <- Reduce(intersect, list(
  colnames(rna_outliers),
  rownames(mutation),
  append(network %>% pull(1), network %>% pull(2))
)) %>% sort()

rna_outliers <- rna_outliers[samples, genes]
mutation <- mutation[genes, samples]


#############################
# Create Temporary Files
#############################

write.table(rna_outliers, "tmp/tmp_PersonaDrive_DEGs.csv", sep = "," )
write_tsv(network, "tmp/tmp_PersonaDrive_network.tsv")
write.csv(mutation, "tmp/tmp_PersonaDrive_MUT.csv", quote = FALSE)


if(!dir.exists(paste0("results/CCLE_",network_choice,"/PersonaDrive/",cell_type))){
  dir.create(paste0("results/CCLE_",network_choice,"/PersonaDrive/",cell_type), recursive = T)
}



