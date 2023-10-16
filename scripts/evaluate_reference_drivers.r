# evaluate_reference_drivers.r

# evaluate_SL_partners.r

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(foreach, quietly = T))
suppressPackageStartupMessages (library(doParallel, quietly = T))
suppressPackageStartupMessages (library(readxl, quietly = T))





# Handling input arguments
option_list = list(
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network being used", metavar ="Network"),
  make_option(c("-a", "--algorithms"), type="character", default="ALL", 
              help="algorithms to include in comparison separated by spaces, or 'ALL' (Default), or lead with '-' to exclude", metavar ="Algorithms"),
  make_option(c("-t", "--threads"), type="integer", default=4, 
              help="Number of threads to use if running in parallel", metavar ="Rank Aggregation Method"),
  make_option(c("-c", "--celltype"), type="character", default="all", 
              help="cell types to analyse, default = 'all'", metavar ="Cell Type")
  
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

max_SL <- opt$synleth_partners
network_choice <- opt$network
algorithms <- opt$algorithms
threads <- opt$threads
cell_type <- opt$celltype


#network_choice <- "own"
#cell_type <- "Liver"
threads <- 10

if(threads>1){
  #registerDoParallel(cores=threads)
  cl <- makeCluster(threads, outfile = "log/evaluate_reference_drivers.log")
  registerDoParallel(cl)
}

#############################
# Sample Info
#############################

sample_info <- read_csv(paste0("validation_data/CCLE_",network_choice,"/sample_info.csv"))
samples <- sample_info$cell_ID %>% sort()

#############################
# Gene List
#############################

genes <- fread(paste0("validation_data/CCLE_",network_choice,"/counts.csv"), select = "gene_ID") %>% pull(gene_ID)


#############################
# Altered Genes
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

altered_genes <- altered_genes %>%
  as.data.frame() %>%
  rownames_to_column("gene_ID") %>%
  pivot_longer(cols = -gene_ID, names_to = "cell_ID", values_to = "altered") %>%
  filter(altered) %>%
  dplyr::select(cell_ID, gene_ID) %>%
  arrange(cell_ID)


#############################
# Read In Results
#############################

#algorithms <- "-combined_de_novo_methods"
algorithms <- "ALL"

suppressWarnings(rm(aggregated_results))
for(alg in list.dirs(paste0("results/CCLE_",network_choice), recursive = F)){
  # First, get names of algorithms that have been run
  alg = str_extract(alg,"(?<=/)[^/]*$")
  # Skip consensus results. These will be added separately
  if(alg=="consensus"){next}
  # Check which algorithms are to be included
  if(algorithms[1]!="ALL"){
    if(all(substr(algorithms,1,1)=="-")){
      if(alg %in% gsub("-","",algorithms)){next}
    }else if(!alg %in% algorithms){
      next
    }else{
      warning("Invalid input for --algorithms. Cannot mix inclusions and exclusions in one statement.")
      stop()
    }
  }
  for(result_file in list.files(paste0("results/CCLE_", network_choice,"/",alg), pattern = "*.csv")){
    message(paste0("Reading result file for ", alg, "-", result_file))
    indiv_result <- data.table::fread(paste0("results/CCLE_", network_choice,"/",alg,"/",result_file))
    if(!"algorithm" %in% colnames(indiv_result)){
      indiv_result <- indiv_result %>% mutate(algorithm = alg)
    }
    if(!exists("aggregated_results")){
      aggregated_results <- indiv_result
    }else{
      aggregated_results <- rbind(aggregated_results,indiv_result)
    }
  }
}



#############################
# Consensus Results
#############################

consensus_drivers <- foreach(result=list.files(paste0("results/CCLE_",network_choice,"/consensus/")), .combine = "rbind") %do% {
  read_csv(paste0("results/CCLE_",network_choice,"/consensus/",result))
}

aggregated_results <- aggregated_results %>% rbind(consensus_drivers)


#############################
# Keywords
#############################

#CGC_keywords <- read_csv("data/cancer_reference_genes/CGC/cancer_gene_census.csv") %>%
#  mutate(keywords = paste(`Tumour Types(Somatic)`,`Tumour Types(Germline)`,`Cancer Syndrome`, sep = " "),
#         source = "CGC") %>%
#  dplyr::select(source, keywords) %>%
#  separate_longer_delim(cols = keywords, delim = " ") %>%
#  mutate(keywords = gsub(",|NA|\\.", "", keywords)) %>%
#  unique()


#NCG_keywords <- read_tsv("data/cancer_reference_genes/NCG/NCG_cancerdrivers_annotation_supporting_evidence.tsv") %>%
#  mutate(keywords = paste(cancer_type,primary_site, sep = " "),
#         source = "NCG") %>%
#  dplyr::select(source, keywords) %>%
#  separate_longer_delim(cols = keywords, delim = " ") %>%
#  mutate(keywords = gsub(",|NA|\\.", "", keywords)) %>%
#  unique()
  
  
#CancerMine_keywords <- read_tsv("data/cancer_reference_genes/CancerMine/cancermine_collated.tsv") %>%
  # Filtering for >1 citation count as done in PersonaDrive
#  filter(citation_count >1) %>%
#  mutate(keywords = paste(cancer_normalized, sep = " "),
#         source = "CancerMine") %>%
#  dplyr::select(source, keywords) %>%
#  separate_longer_delim(cols = keywords, delim = " ") %>%
#  mutate(keywords = gsub(",|NA|\\.", "", keywords)) %>%
#  unique()

#all_keywords <- rbind(CGC_keywords,NCG_keywords,CancerMine_keywords) %>%
#  dplyr::select(keywords) %>%
#  unique()


#write_csv(all_keywords,"data/cancer_reference_genes//raw_keywords.csv")

all_keywords <- read_xlsx("data/cancer_reference_genes/keywords.xlsx") %>%
  pivot_longer(cols = everything(), names_to = "lineage", values_to = "keyword", values_drop_na = T)



#############################
# Reference Drivers
#############################

CGC_all <- read_csv("data/cancer_reference_genes/CGC/cancer_gene_census.csv") %>%
  dplyr::select(gene_ID = `Gene Symbol`) %>%
  unique()

CGC_specific <- read_csv("data/cancer_reference_genes/CGC/cancer_gene_census.csv") %>%
  pivot_longer(cols = c(`Tumour Types(Somatic)`,`Tumour Types(Germline)`,`Cancer Syndrome`), names_to = "keyword_type", values_to = "keyword") %>%
  separate_longer_delim(cols = keyword, delim = ",") %>%
  mutate(keyword = gsub("^ ", "", keyword)) %>%
  dplyr::select(gene_ID = `Gene Symbol`, keyword) %>%
  rbind(
    read_csv("data/cancer_reference_genes/CGC/cancer_gene_census.csv") %>%
      mutate(keyword = paste(`Tumour Types(Somatic)`,`Tumour Types(Germline)`,`Cancer Syndrome`, sep = " ")) %>%
      separate_longer_delim(cols = keyword, delim = " ") %>%
      mutate(keyword = gsub(",|NA|\\.", "", keyword)) %>%
      dplyr::select(gene_ID = `Gene Symbol`, keyword)
  ) %>%
  left_join(all_keywords, by = "keyword", multiple = "all") %>%
  dplyr::select(gene_ID, lineage) %>%
  na.omit() %>%
  unique()


NCG_all <- read_tsv("data/cancer_reference_genes/NCG/NCG_cancerdrivers_annotation_supporting_evidence.tsv") %>%
  dplyr::select(gene_ID = symbol) %>%
  unique()

NCG_specific <- read_tsv("data/cancer_reference_genes/NCG/NCG_cancerdrivers_annotation_supporting_evidence.tsv") %>%
  pivot_longer(cols = c(cancer_type,primary_site), names_to = "keyword_type", values_to = "keyword") %>%
  separate_longer_delim(cols = keyword, delim = ",") %>%
  mutate(keyword = gsub("^ ", "", keyword)) %>%
  dplyr::select(gene_ID = symbol, keyword) %>%
  rbind(
    read_tsv("data/cancer_reference_genes/NCG/NCG_cancerdrivers_annotation_supporting_evidence.tsv") %>%
      mutate(keyword = paste(cancer_type,primary_site, sep = " ")) %>%
      separate_longer_delim(cols = keyword, delim = " ") %>%
      mutate(keyword = gsub(",|NA|\\.", "", keyword)) %>%
      dplyr::select(gene_ID = symbol, keyword)
    ) %>%
  left_join(all_keywords, by = "keyword", multiple = "all") %>%
  dplyr::select(gene_ID, lineage) %>%
  na.omit() %>%
  unique()


CancerMine_all <- read_tsv("data/cancer_reference_genes/CancerMine/cancermine_collated.tsv") %>%
  # Filtering for >1 citation count as done in PersonaDrive
  filter(citation_count >1) %>%
  dplyr::select(gene_ID = gene_normalized) %>%
  unique()

CancerMine_specific <- read_tsv("data/cancer_reference_genes/CancerMine/cancermine_collated.tsv") %>%
  # Filtering for >1 citation count as done in PersonaDrive
  filter(citation_count >1) %>%
  dplyr::select(gene_ID = gene_normalized, keyword = cancer_normalized) %>%
  rbind(
    read_tsv("data/cancer_reference_genes/CancerMine/cancermine_collated.tsv") %>%
      separate_longer_delim(cancer_normalized, delim = " ") %>%
      mutate(keyword = gsub(",|NA|\\.", "", cancer_normalized)) %>%
      dplyr::select(gene_ID = gene_normalized, keyword)
  ) %>%
  left_join(all_keywords, by = "keyword", multiple = "all") %>%
  dplyr::select(gene_ID, lineage) %>%
  na.omit() %>%
  unique()

reference_drivers_all <- rbind(
  CGC_all %>% mutate(source="CGC"),
  NCG_all %>% mutate(source = "NCG"),
  CancerMine_all %>% mutate(source = "CancerMine")
) %>%
  filter(gene_ID %in% genes)

reference_drivers_specific <- rbind(
  CGC_specific %>% mutate(source="CGC"),
  NCG_specific %>% mutate(source = "NCG"),
  CancerMine_specific %>% mutate(source = "CancerMine")
) %>%
  filter(gene_ID %in% genes)

rm(CGC_all,CGC_specific,NCG_all,NCG_specific,CancerMine_all,CancerMine_specific,all_keywords)


#############################
# Results
#############################
if(cell_type != "all"){
  aggregated_results <- aggregated_results %>% filter(lineage == cell_type)
}


run_evaluation <- function(ref_set,gold_standard_type){
  
  
  if(gold_standard_type == "all"){
    gs <- reference_drivers_all %>% filter(source==ref_set) %>% pull(gene_ID) %>% unique()
    samples <- aggregated_results$cell_ID %>% unique()
  } else if(gold_standard_type == "specific"){
    
    if(ref_set == "all"){
      gold_standard <- reference_drivers_specific
    }else{
      gold_standard <- reference_drivers_specific %>% filter(source==ref_set)
    }
    
    samples <- aggregated_results %>% filter(lineage %in% gold_standard$lineage) %>% pull(cell_ID) %>% unique()
  }
  
  
  
  
  if(!file.exists(paste0("results/CCLE_",network_choice,"/full_results_", ref_set,"_", gold_standard_type, "_reference_driver_level.csv"))){
    
    prediction_stats <- foreach(sample=samples, .combine = "rbind", .packages = c("tidyverse","foreach"), .export = c("aggregated_results","altered_genes")) %dopar% {
      
      lin <- aggregated_results %>% filter(cell_ID==sample) %>% pull("lineage") %>% unique()
      
      alt_genes <- altered_genes %>% filter(cell_ID==sample) %>% pull(gene_ID)
      
      if(gold_standard_type == "specific"){
        gs_specific <- gold_standard %>% filter(lineage==lin) %>% pull(gene_ID) %>% unique()
        gs_specific <- gs_specific[which(gs_specific %in% alt_genes)]
      }else{
        gs_specific <- gs[which(gs %in% alt_genes)]
      }
      
      
      
      
      
      tmp1 <- aggregated_results %>% filter(cell_ID == sample)
      
      foreach(alg=unique(tmp1$algorithm), .combine = "rbind", .packages = c("tidyverse","foreach"), .export = "run_evaluation") %do% {
        
        tmp2 <- tmp1 %>% filter(algorithm == alg)
        if(nrow(tmp2)==0){return(NULL)}else{
        max_drivers <- max(tmp2 %>% pull(rank))
        
        print(paste0("Calculating stats for ", sample, "(",which(samples == sample),"/",length(samples),")", " targets using ", alg))
        
        foreach(n=1:max_drivers, .combine = "rbind", .packages = c("tidyverse"), .export = "run_evaluation") %do% {
          
          predicted <- tmp2 %>% filter(rank <= n) %>% pull(driver)
          correct <- predicted[which(predicted %in% gs_specific)]
          wrong <- predicted[which(!predicted %in% gs_specific)]
          TP <- length(correct)
          FP <- length(wrong)
          precision <- TP/n
          recall <- TP/length(gs_specific)
          F1 <- 2*((precision*recall)/(precision + recall))
          if(is.nan(F1)){
            F1 <- 0
          }
          n_gs <- length(gs_specific)
          
          data.frame(lineage=lin,
                     cell_ID = sample, 
                     algorithm = alg, 
                     n = n,
                     correct = paste0(correct, collapse = ";"), 
                     TP = TP,
                     FP = FP,
                     precision = precision,
                     recall = recall,
                     F1 = F1,
                     n_gs = n_gs)
          
        }
        }
        
        
        
      }
      
    }
    write_csv(prediction_stats, paste0("results/CCLE_",network_choice,"/full_results_", ref_set,"_", gold_standard_type, "_reference_driver_level.csv"))
  }else{
    message("SL Prediction stats file already exists. Reading in previous results.")
    message(paste0("To stop this, remove file: ", "results/CCLE_",network_choice,"/full_results_", ref_set,"_", gold_standard_type, "_reference_driver_level.csv"))
    prediction_stats <- read_csv(paste0("results/CCLE_",network_choice,"/full_results_", ref_set,"_", gold_standard_type, "_reference_driver_level.csv"))
  }
  
}

#run_evaluation("CGC","specific")
run_evaluation("CGC","all")
run_evaluation("all","specific")









make_evaluation_plot <- function(ref_set,gold_standard_type,algs){
  
  prediction_stats <- read_csv(paste0("results/CCLE_",network_choice,"/full_results_", ref_set,"_", gold_standard_type, "_reference_driver_level.csv"))
  
  # Dynamic N-Max
  
  gs_lengths <- prediction_stats %>% dplyr::select(cell_ID,n_gs) %>% unique() %>% pull(n_gs)
    
  ## Calculates 2*median number of sensitive genes for cell lines with >3 sensitive genes
  N_max <- gs_lengths[which(gs_lengths  > 3 )] %>% median()*2
  
  
  
  prediction_stats_summarised <- prediction_stats %>%
    group_by(algorithm,n) %>%
    filter(n()>5) %>%
    summarise(
      mean_TP = mean(TP),
      mean_precision = mean(precision),
      mean_recall = mean(recall),
      mean_F1 = median(F1),
    ) %>%
    pivot_longer(cols = -c(algorithm,n), names_to = "measure", values_to = "value") %>%
    mutate(measure = factor(measure, levels = c("mean_TP","mean_recall","mean_precision","mean_F1")))
  
  
  ggplot(prediction_stats_summarised %>% filter(algorithm %in% algs) %>% filter(n<=N_max), aes(x = n, y = value, color = algorithm)) +
    geom_line() +
    ylab("Average Value") +
    xlab("Number of Predicted Sensitive Genes") +
    theme_light() +
    facet_wrap(~measure, nrow = 1, scales = "free")
  
  ggsave(paste0("results/CCLE_",network_choice,"/",ref_set,"_",gold_standard_type,"_reference_driver_level_stats.png"),width = 50, height = 20, units = "cm", dpi = 300)
  
}


alg_of_interest <- c(
  #"consensus_top_mean",
  #"consensus_top_median",
  "consensus_topklists_geomean",
  #"consensus_top_l2norm",
  #"consensus_BiG",
  
  "DawnRank",
  "OncoImpact",
  "PNC",
  "PRODIGY",
  "PersonaDrive",
  "SCS",
  "sysSVM2",
  "PhenoDriverR",
  
  "CSN_DFVS",
  "CSN_MDS",
  "CSN_MMS",
  "CSN_NCUA",
  "LIONESS_DFVS",
  "LIONESS_MDS",
  "LIONESS_MMS",
  "LIONESS_NCUA",
  "SPCC_DFVS",
  "SPCC_MDS",
  "SPCC_MMS",
  "SPCC_NCUA",
  "SSN_DFVS",
  "SSN_MDS",
  "SSN_MMS",
  "SSN_NCUA"
)

make_evaluation_plot("CGC","all",alg_of_interest)

make_evaluation_plot("all","specific",alg_of_interest)















test <- read_csv("results/CCLE_own/full_results_CGC_all_reference_driver_level.csv") %>%
  group_by(algorithm,n) %>%
  filter(n()>5) %>%
  summarise(
    mean_TP = mean(TP),
    mean_precision = mean(precision),
    mean_recall = mean(recall),
    mean_F1 = median(F1),
  ) %>%
  pivot_longer(cols = -c(algorithm,n), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure, levels = c("mean_TP","mean_recall","mean_precision","mean_F1"))) %>%
  filter(n==10, measure=="mean_F1")


test2 <- read_csv("results/CCLE_own/full_results_all_specific_reference_driver_level.csv") %>%
  group_by(algorithm,n) %>%
  filter(n()>5) %>%
  summarise(
    mean_TP = mean(TP),
    mean_precision = mean(precision),
    mean_recall = mean(recall),
    mean_F1 = median(F1),
  ) %>%
  pivot_longer(cols = -c(algorithm,n), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure, levels = c("mean_TP","mean_recall","mean_precision","mean_F1"))) %>%
  filter(n==10, measure=="mean_F1")

