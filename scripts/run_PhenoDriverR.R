# run_PhenoDriverR.R

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(AnnotationDbi, quietly = T))
suppressPackageStartupMessages (library(org.Hs.eg.db, quietly = T))
suppressPackageStartupMessages (library(DESeq2, quietly = T))
suppressPackageStartupMessages (library(gtools, quietly = T))
suppressPackageStartupMessages (library(stats, quietly = T))
suppressPackageStartupMessages (library(doParallel, quietly = T))
suppressPackageStartupMessages (library(foreach, quietly = T))
suppressPackageStartupMessages (library(clusterProfiler, quietly = T))


# Handling input arguments
option_list = list(
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network to use", metavar ="Network"),
  make_option(c("-c", "--celltype"), type="character", default="Liver", 
              help="cell type to analyse", metavar ="Network"),
  make_option(c("-t", "--threads"), type="numeric", default=4, 
              help="Number of threads to use", metavar ="Threads")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

network_choice <- opt$network
cell_type <- opt$celltype
threads <- opt$threads



#############################
# Functions
#############################

for(script in list.files("scripts/PhenoDriverR/")[which(str_detect(list.files("scripts/PhenoDriverR/"),"\\.R|\\.r"))]){
  source(paste0("scripts/PhenoDriverR/",script))
}

#############################
# Sample Info
#############################

sample_info <- read_csv(paste0("validation_data/CCLE_",network_choice,"/sample_info.csv")) %>% filter(lineage==cell_type)
samples <- sample_info$cell_ID %>% sort()

#############################
# Prepare Input Data
#############################
message("Loading Data")

# Network

if(network_choice=="own"){
  STN <- fread("data/own_networks/PhenoDriverR.txt")
}else{

  ## Only keeping edges with known excitatory vs inhibitory information
  original_network <- fread(paste0("validation_data/CCLE_",network_choice,"/network_directed.csv")) %>% select(Source_nodes = protein_1, Target_nodes = protein_2)
  original_STN <- read.delim("scripts/PhenoDriverR/SignalTransductionNet.txt")
  STN <- inner_join(original_STN, original_network, by = c("Source_nodes","Target_nodes"))

}

# RNA

rna <- fread(paste0("validation_data/CCLE_",network_choice,"/counts.csv"), select = c("gene_ID",samples)) %>% 
  arrange(gene_ID) %>%
  column_to_rownames("gene_ID") 

# Pseudonormal RNA

suppressWarnings(rm(pseudonormal_rna))
for(sample in samples){
  # Collect expression data for all samples other than the current sample
  other_rna <- rna[,samples[samples != sample]]
  # Collect mean expression for other samples as the pseudo-normal expression for each cell
  tmp_pseudonormal_rna <- enframe(round(rowMeans(other_rna),0), name = "gene_ID", value = sample) %>% column_to_rownames("gene_ID")
  
  
  if(!exists("pseudonormal_rna")){
    pseudonormal_rna <- tmp_pseudonormal_rna
  }else{
    pseudonormal_rna <- cbind(pseudonormal_rna, tmp_pseudonormal_rna)
  }
  
}

## Combine cancer and pseduonormal expression data
rna <- rna %>% rename_with(cols = everything(), ~ paste0(.x,"_t")) %>% rownames_to_column("gene_ID")
pseudonormal_rna <- pseudonormal_rna %>% rename_with(cols = everything(), ~ paste0(.x,"_n")) %>% rownames_to_column("gene_ID")
allexp <- merge(rna, pseudonormal_rna)


# Mutation

mutation <- fread(paste0("validation_data/CCLE_",network_choice,"/mutations_MAF.csv")) %>%
  filter(cell_ID %in% samples)

# Pathways

reactome <- read.delim('data/NCBI2Reactome.txt', header = F)

reactome <- reactome[reactome$V6 == 'Homo sapiens', 1:4]
reactome <- unique(reactome)
reactome$V1 <- AnnotationDbi::select(org.Hs.eg.db, keys = reactome$V1, columns = 'SYMBOL')$SYMBOL
reactome <- unique(reactome)
reactome <- reactome[!(is.na(reactome$V1)),]
reactome <- reactome[reactome$V4 %in% names(table(reactome$V4)[table(reactome$V4) > 10 &
                                                                 table(reactome$V4) < 200]),]


#############################
# PhenoDrive
#############################

message("Finding DEGs")
# DEGs

NTindex <- match(colnames(pseudonormal_rna)[-1:-1],
                 colnames(allexp))
TPindex <- match(colnames(rna)[-1:-1],
                 colnames(allexp))
allexp <- allexp[(rowSums(allexp[NTindex] > 20) >= (0.8 * length(NTindex))) |
                   (rowSums(allexp[TPindex] > 20) >= (0.8 * length(TPindex))), ]

expressionData <- list(Normal = allexp[, c(1, NTindex)], Tumor = allexp[, c(1, TPindex)])


diffexpgene <- DEA(exp = expressionData,
                   parallelworker = threads,
                   genenamecol = 1,
                   annotationcol = 1)

# Pathway Enrichment

message("Finding enriched pathways")

## Get personalized abnormal pathways
enrichreactomeRes <- pAbnormalPathway(diff = diffexpgene,
                                      reactome = reactome,
                                      parallelworker = threads)

## Calculating driving force matrix for each patient
genescores <- calDrivingForce(network = STN,
                              diffexpgene = diffexpgene,
                              enrichreactomeRes = enrichreactomeRes,
                              reactome = reactome,
                              parallelworker = threads)

# Getting Driver Genes

message("Prioritising Drivers")

drivingforcelist <- genescores$genescore
topgenes_ind = 100

maf <- mutation

genescore_individual <- drivingforcelist[!(sapply(drivingforcelist, is.null))]

sample_snp  <- mutation$cell_ID

sample_score <- gsub("_t","",names(genescore_individual))

genescore_individual_snp <- genescore_individual
for (i in 1:length(genescore_individual_snp)) {
  if (sample_score[i] %in% sample_snp) {
    temp <- maf[sample_snp %in% sample_score[i],]
    genescore_individual_snp[[i]] <- genescore_individual_snp[[i]][rownames(genescore_individual_snp[[i]]) %in%
                                                                     temp$HugoSymbol,]
  }
  else genescore_individual_snp[[i]] <- 1
}
genescore_individual_snp <- genescore_individual_snp[sample_score %in% sample_snp]
for (i in 1:length(genescore_individual_snp)) {
  if (is.null(dim(genescore_individual_snp[[i]]))) next
  else genescore_individual_snp[[i]] <- rowSums(genescore_individual_snp[[i]])
}
genescore_individual_snp <- lapply(genescore_individual_snp,
                                   function(x) return(x[order(abs(x),decreasing = T)]))

drivergene_ind <- lapply(genescore_individual_snp, function(x) return(x[1:topgenes_ind]))




topgenes <- drivergene_ind %>% 
  lapply(sort)

suppressWarnings(rm(result))
for(i in 1:length(topgenes)){
  indiv_result <- topgenes[i] %>%
    as.data.frame() %>%
    rownames_to_column("gene_ID") %>%
    mutate(rank=row_number(), lineage = cell_type, cell_ID = gsub("_t","",names(topgenes)[i])) %>%
    dplyr::select(lineage, cell_ID, driver = gene_ID, rank)
  if(!exists("result")){
    result <- indiv_result
  }else{
    result <- rbind(result, indiv_result)
  }
}

#############################
# Save Result
#############################
if(!dir.exists(paste0("results/CCLE_",network_choice,"/PhenoDriverR"))){
  dir.create(paste0("results/CCLE_",network_choice,"/PhenoDriverR"))
}

write_csv(result, paste0("results/CCLE_",network_choice,"/PhenoDriverR/",cell_type,".csv"))



