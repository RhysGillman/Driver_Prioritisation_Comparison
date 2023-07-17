#evaluate_compare_outputs.r

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(UpSetR, quietly = T))


# Handling input arguments
option_list = list(
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network to use", metavar ="Network"),
  make_option(c("-a", "--algorithms"), type="character", default="ALL", 
              help="algorithms to include in comparison separated by spaces, or 'ALL' (Default), or lead with '-' to exclude", metavar ="Algorithms")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

network_choice <- opt$network
algorithms <- opt$algorithms

algorithms <- str_split(algorithms, " ") %>% unlist()


#############################
# Sample Info
#############################

sample_info <- read_csv(paste0("validation_data/CCLE_",network_choice,"/sample_info.csv"))
samples <- sample_info$cell_ID %>% sort()


#############################
# Read In Results
#############################

algorithms <- "-combined_de_novo_methods"

suppressWarnings(rm(aggregated_results))
for(alg in list.dirs(paste0("results/CCLE_",network_choice), recursive = F)){
  # First, get names of algorithms that have been run
  alg = str_extract(alg,"(?<=/)[^/]*$")
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

aggregated_results <- aggregated_results %>% filter(lineage=="Liver")

#############################
# Cosine Similarity
#############################
top_ns <- c(5)
alg_pairs <- as.list(as.data.frame(combn(unique(aggregated_results$algorithm), 2)))

suppressWarnings(rm(combined_cos_sim))
for(alg_pair in alg_pairs){
  
  alg1 <- alg_pair[1]
  alg2 <- alg_pair[2]
  
  tmp_results_1 <- aggregated_results %>% filter(algorithm %in% c(alg1, alg2))
  
  comp_cells <- intersect(tmp_results_1 %>% filter(algorithm==alg1) %>% pull(cell_ID),
                          tmp_results_1 %>% filter(algorithm==alg2) %>% pull(cell_ID)
  )
  
  
  for(cell in comp_cells){
    
    tmp_results_2 <- tmp_results_1 %>% filter(cell_ID == cell)
    
    print(paste("Calculating cosine similarity for", cell, "in", alg1, "vs", alg2, sep = " "))
    
    for(n in top_ns){
      
      tmp_results_3 <- tmp_results_2 %>% filter(rank <= n) %>% 
        ungroup() %>% 
        dplyr::select(driver, algorithm) %>%
        table() %>%
        as.data.frame() %>%
        pivot_wider(names_from = algorithm, values_from = Freq)
      
      cos_sim <- (
        # This calculates the sum of the products of the two vectors
        as.numeric(t(pull(tmp_results_3,2))%*%pull(tmp_results_3,3))
        # Divide
      ) /
        # This multiplies the square root of the sum squares of each vector
        (sqrt(sum(pull(tmp_results_3,2)^2))*sqrt(sum(pull(tmp_results_3,3)^2)))
      
      indiv_result <- data.frame(cell_ID = cell, 
                                 algorithm_1 = alg1, 
                                 algorithm_2 = alg2,
                                 top_n = n,
                                 cos_sim = cos_sim
      )
      
      if(!exists("combined_cos_sim")){
        combined_cos_sim <- indiv_result
      }else{
        combined_cos_sim <- rbind(combined_cos_sim, indiv_result)
        
      }
      
    }
  }
}

summarised_cos_sim <- combined_cos_sim %>%
  group_by(algorithm_1,algorithm_2) %>%
  summarise(mean_cos_sim = mean(cos_sim)) %>%
  mutate(mean_cos_dist = 1-mean_cos_sim)

dist_matrix <- summarised_cos_sim %>%
  dplyr::select(algorithm_1,algorithm_2,mean_cos_dist) %>%
  pivot_wider(names_from = algorithm_2, values_from = mean_cos_dist)


