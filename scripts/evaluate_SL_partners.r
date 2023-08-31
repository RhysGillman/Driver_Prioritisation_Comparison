# evaluate_SL_partners.r

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(foreach, quietly = T))
suppressPackageStartupMessages (library(doParallel, quietly = T))






# Handling input arguments
option_list = list(
  make_option(c("-s", "--synleth_partners"), type="character", default=5, 
              help="Maximum number of SL-partners to use for each gene", metavar ="Synthetic Lethal Partners"),
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network being used", metavar ="Network"),
  make_option(c("-a", "--algorithms"), type="character", default="ALL", 
              help="algorithms to include in comparison separated by spaces, or 'ALL' (Default), or lead with '-' to exclude", metavar ="Algorithms"),
  make_option(c("-t", "--threads"), type="integer", default=4, 
              help="Number of threads to use if running in parallel", metavar ="Rank Aggregation Method")
  
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

max_SL <- opt$synleth_partners
network_choice <- opt$network
algorithms <- opt$algorithms
threads <- opt$threads

if(threads>1){
  #registerDoParallel(cores=threads)
  cl <- makeCluster(threads, outfile = "log/evaluate_SL_partners.log")
  registerDoParallel(cl)
}

#############################
# Sample Info
#############################

sample_info <- read_csv(paste0("validation_data/CCLE_",network_choice,"/sample_info.csv"))
samples <- sample_info$cell_ID %>% sort()


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

#aggregated_results <- aggregated_results %>% filter(lineage=="Liver")


#############################
# Consensus Results
#############################

consensus_drivers <- foreach(result=list.files(paste0("results/CCLE_",network_choice,"/consensus/")), .combine = "rbind") %do% {
  read_csv(paste0("results/CCLE_",network_choice,"/consensus/",result))
}

aggregated_results <- aggregated_results %>% rbind(consensus_drivers)

#aggregated_results <- aggregated_results %>% filter(str_detect(algorithm,"consensus"))

#############################
# Read In Gold Standards
#############################


gold_standard <- read_csv(paste0("validation_data/CCLE_",network_choice,"/gold_standards.csv")) %>%
  dplyr::select(cell_ID,gene_ID) %>%
  unique() %>%
  group_by(cell_ID) %>%
  summarise(sensitive_genes = list(gene_ID)) %>%
  deframe()


#############################
# Read In Synthetic Lethality Predictions
#############################

SL <- read_csv(paste0("validation_data/CCLE_",network_choice,"/SL_partners_max_",max_SL,".csv")) %>%
  dplyr::select(gene_ID=gene1,partner=gene2,SL_rank=rank)


#############################
# Read In LOF/GOF Predictions
#############################

LOF_GOF <- read_csv("data/LOF_GOF_annotations.csv")

#############################
# N Max
#############################


# Calculates 2*median number of sensitive genes for cell lines with >3 sensitive genes
N_max <- gold_standard[lapply(gold_standard, length) > 3] %>% lapply(length) %>% unlist() %>% median()*2



####################################
# Driver-SL-Partner-Level Evaluation
####################################

SL_partner_results <- aggregated_results %>%
  dplyr::rename(driver_rank=rank) %>%
  left_join(LOF_GOF, by = c("driver"="gene_ID","cell_ID")) %>%
  # assume that missing annotations are LOF
  mutate(annotation=ifelse(is.na(annotation),"LOF",annotation)) %>%
  # If mutation is GOF, add the gene itself to the target list
  mutate(target = ifelse(annotation %in% c("GOF"), driver, NA)) %>%
  # Add SL partners
  left_join(SL, by = c("driver"="gene_ID"), multiple="all") %>%
  # Remove SL partners for GOF mutations
  mutate(partner = ifelse(annotation == "GOF", NA, partner),
         SL_rank = ifelse(annotation == "GOF", NA, SL_rank)) %>%
  unique() %>%
  # Make the SL partner the target for LOF and both
  mutate(target = ifelse(annotation %in% c("LOF", "both"), partner, target))

# Add the original gene_ID as an additional target for genes with both GOF and LOF mutations
SL_partner_results <- SL_partner_results %>%
  rbind(SL_partner_results %>% 
          filter(annotation == "both") %>% 
          mutate(target = driver, SL_rank = 0)) %>%
  arrange(driver_rank,SL_rank) %>%
  #filter(!is.na(target)) %>%
  group_by(cell_ID, algorithm) %>%
  # For each cell, after arranging by driver rank, then SL_rank, remove duplicated target genes
  filter(!duplicated(target)) %>%
  # Finally remove missing targets (LOF mutations without an SL partner)
  filter(!is.na(target)) %>%
  mutate(final_rank = row_number()) %>%
  ungroup() %>%
  arrange(lineage,algorithm, cell_ID, final_rank)

stats <- SL_partner_results %>%
  group_by(lineage,cell_ID,algorithm) %>%
  summarise(number_of_targets=n())



tmp_samples <- intersect(SL_partner_results$cell_ID, names(gold_standard))

if(!file.exists(paste0("results/CCLE_",network_choice,"/full_results_SL_partner_level.csv"))){
  
SL_prediction_stats <- foreach(sample=tmp_samples, .combine = "rbind", .packages = c("tidyverse","foreach")) %dopar% {
  
  lineage <- sample_info %>% filter(cell_ID==sample) %>% pull("lineage")
  gs <- gold_standard[sample] %>% unlist()
  tmp1 <- SL_partner_results %>% filter(cell_ID == sample)
  
  foreach(alg=unique(tmp1$algorithm), .combine = "rbind", .packages = c("tidyverse","foreach")) %do% {
    
    tmp2 <- tmp1 %>% filter(algorithm == alg)
    if(nrow(tmp2)==0){break}
    max_drivers <- max(tmp2 %>% pull(final_rank))
    
    print(paste0("Calculating stats for ", sample, "(",which(tmp_samples == sample),"/",length(tmp_samples),")", " targets using ", alg))
    
    foreach(n=1:max_drivers, .combine = "rbind", .packages = c("tidyverse")) %do% {
      
      predicted <- tmp2 %>% filter(final_rank <= n) %>% pull(target)
      correct <- predicted[which(predicted %in% gs)]
      wrong <- predicted[which(!predicted %in% gs)]
      TP <- length(correct)
      FP <- length(wrong)
      PPV <- TP/(TP+FP)
      precision <- TP/n
      recall <- TP/length(gs)
      F1 <- 2*((precision*recall)/(precision + recall))
      if(is.nan(F1)){
        F1 <- 0
      }
      
      data.frame(lineage=lineage,
                 cell_ID = sample, 
                 algorithm = alg, 
                 n = n,
                 correct = paste0(correct, collapse = ";"), 
                 TP = TP,
                 FP = FP,
                 PPV = PPV,
                 precision = precision,
                 recall = recall,
                 F1 = F1)
      
      
    }
    
    
    
  }
  
}
write_csv(SL_prediction_stats, paste0("results/CCLE_",network_choice,"/full_results_SL_partner_level.csv"))
}else{
  message("SL Prediction stats file already exists. Reading in previous results.")
  message(paste0("To stop this, remove file: ", "results/CCLE_",network_choice,"/full_results_SL_partner_level.csv"))
  SL_prediction_stats <- read_csv(paste0("results/CCLE_",network_choice,"/full_results_SL_partner_level.csv"))
}




alg_of_interest <- c(
  #"consensus_top_mean",
  #"consensus_top_median",
  "consensus_top_geomean",
  #"consensus_top_l2norm"
  
  "DawnRank",
  "OncoImpact",
  "PNC",
  "PRODIGY",
  "PersonaDrive",
  "SCS"
  
  #"CSN_DFVS",
  #"CSN_MDS",
  #"CSN_MMS",
  #"CSN_NCUA",
  #"LIONESS_DFVS",
  #"LIONESS_MDS",
  #"LIONESS_MMS",
  #"LIONESS_NCUA",
  #"SPCC_DFVS",
  #"SPCC_MDS",
  #"SPCC_MMS",
  #"SPCC_NCUA",
  #"SSN_DFVS",
  #"SSN_MDS",
  #"SSN_MMS",
  #"SSN_NCUA"
)





SL_prediction_stats_summarised <- SL_prediction_stats %>%
  group_by(algorithm,n) %>%
  filter(n()>10) %>%
  summarise(
    mean_TP = mean(TP),
    mean_PPV = mean(PPV),
    mean_precision = mean(precision),
    mean_recall = mean(recall),
    mean_F1 = median(F1),
    #median_TP = median(TP),
    #median_PPV = median(PPV),
    #median_precision = median(precision),
    #median_recall = median(recall),
    #median_F1 = median(F1)
  ) %>%
  pivot_longer(cols = -c(algorithm,n), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure, levels = c("mean_TP","mean_PPV","mean_recall","mean_precision","mean_F1")))

ggplot(SL_prediction_stats_summarised %>% filter(algorithm %in% alg_of_interest) %>% filter(n<=100), aes(x = n, y = value, color = algorithm)) +
  geom_line() +
  ylab("Average Value") +
  xlab("Number of Predicted Sensitive Genes") +
  theme_light() +
  facet_wrap(~measure, nrow = 1, scales = "free")

ggsave(paste0("results/CCLE_",network_choice,"/SL_partner_level_stats.png"),width = 50, height = 20, units = "cm", dpi = 300)


SL_prediction_stats_summarised <- SL_prediction_stats %>%
  group_by(lineage,algorithm,n) %>%
  summarise(
    mean_TP = mean(TP),
    mean_PPV = mean(PPV),
    mean_precision = mean(precision),
    mean_recall = mean(recall),
    mean_F1 = median(F1)
  ) %>%
  pivot_longer(cols = -c(lineage,algorithm,n), names_to = "measure", values_to = "value")

ggplot(SL_prediction_stats_summarised %>% 
         filter(algorithm %in% alg_of_interest) %>% 
         filter(n<=100), 
       aes(x = n, y = value, color = algorithm)) +
  geom_line() +
  ylab("Average Value") +
  xlab("Number of Predicted Sensitive Genes") +
  theme_light() +
  facet_wrap(lineage ~ measure, scales = "free", ncol = 5)
  #facet_grid(rows = vars(lineage), cols = vars(measure), scales = "free")

ggsave(paste0("results/CCLE_",network_choice,"/SL_partner_level_stats_by_lineage.png"),width = 30, height = 100, units = "cm", dpi = 300)




























### Original

#N_max <- max(stats$number_of_targets)


suppressWarnings(rm(SL_prediction_stats))

tmp_samples <- intersect(SL_partner_results$cell_ID, names(gold_standard))

for(sample in tmp_samples){
  lineage <- sample_info %>% filter(cell_ID==sample) %>% pull("lineage")
  gs <- gold_standard[sample] %>% unlist()
  tmp1 <- SL_partner_results %>% filter(cell_ID == sample)
  
  for(alg in unique(SL_partner_results$algorithm)){
    
    tmp2 <- tmp1 %>% filter(algorithm == alg)
    
    if(nrow(tmp2)==0){break}
    
    max_drivers <- max(tmp2 %>% pull(final_rank))
    
    print(paste0("Calculating stats for ", sample, "(",which(tmp_samples == sample),"/",length(tmp_samples),")", " targets using ", alg))
    
    for(n in 1:max_drivers){
      
      
      predicted <- tmp2 %>% filter(final_rank <= n) %>% pull(target)
      correct <- predicted[which(predicted %in% gs)]
      wrong <- predicted[which(!predicted %in% gs)]
      TP <- length(correct)
      FP <- length(wrong)
      PPV <- TP/(TP+FP)
      precision <- TP/n
      recall <- TP/length(gs)
      F1 <- 2*((precision*recall)/(precision + recall))
      if(is.nan(F1)){
        F1 <- 0
      }
      indiv_result <- data.frame(lineage=lineage,
                                 cell_ID = sample, 
                                 algorithm = alg, 
                                 n = n,
                                 correct = paste0(correct, collapse = ";"), 
                                 TP = TP,
                                 FP = FP,
                                 PPV = PPV,
                                 precision = precision,
                                 recall = recall,
                                 F1 = F1)
      if(!exists("SL_prediction_stats")){
        SL_prediction_stats <- indiv_result
      }else{
        SL_prediction_stats <- rbind(SL_prediction_stats, indiv_result)
      }
    }
    
  }
  
}


write_csv(SL_prediction_stats, paste0("results/CCLE_",network_choice,"/full_results_SL_partner_level.csv"))


SL_prediction_stats <- read_csv(paste0("results/CCLE_",network_choice,"/full_results_SL_partner_level.csv"))