#get_consensus_drivers.r

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(TopKLists, quietly = T))
suppressPackageStartupMessages (library(foreach, quietly = T))
suppressPackageStartupMessages (library(doParallel, quietly = T))
suppressPackageStartupMessages (library(truncnorm, quietly = T))


# Code sources

# BiG is not longer available as an R package on CRAN
# Using code from https://github.com/baillielab/comparison_of_RA_methods/tree/main/algorithms
# BiG_code_platform_changed.R
# useBiG.R

source("scripts/BiG_code_platform_changed.R")


# Handling input arguments
option_list = list(
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network being used", metavar ="Network"),
  make_option(c("-a", "--algorithms"), type="character", default="ALL", 
              help="algorithms to include in consensus separated by spaces, or 'ALL' (Default), or lead with '-' to exclude", metavar ="Algorithms"),
  make_option(c("-m", "--method"), type="character", default="topklists", 
              help="Rank aggregation method to be used", metavar ="Rank Aggregation Method"),
  make_option(c("-t", "--threads"), type="integer", default=4, 
              help="Number of threads to use if running in parallel", metavar ="Rank Aggregation Method")
  
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

network_choice <- opt$network
algorithms <- opt$algorithms
RAmethod <- opt$method
threads <- opt$threads

if(threads>1){
  cl <- makeCluster(threads, outfile = "log/get_consensus_drivers.log")
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

algorithms <- "-combined_de_novo_methods"
#algorithms <- "ALL"

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
# Rank Aggregation
#############################

if(RAmethod=="topklists"){
  
  algs <- sort(unique(aggregated_results$algorithm))
  
  # For each cell, get a list of the ranked drivers from each algorithm
  
  start <- Sys.time()
  
  consensus_drivers <- foreach(cell=unique(aggregated_results$cell_ID), .packages = c("tidyverse","TopKLists","foreach"), .combine = "rbind") %dopar% {
    
    
    
    #message(paste0("Getting consensus drivers for ",cell))
    
    indiv_ranks <- foreach(alg=algs, .packages = "tidyverse") %do% {
      x <- aggregated_results %>% 
        filter(cell_ID==cell) %>%
        filter(algorithm==alg) %>%
        pull(driver)
      
      # NA indicates there were no results from that algorithm for that cell_ID. In this case, replace with empty vector
      if(any(is.na(x))){
        x <- c()
      }
      
      x
      
    }
    
    names(indiv_ranks) <- algs
    
    
    lineage <- sample_info %>% filter(cell_ID==cell) %>% pull(lineage)
    
    topk <- TopKLists::Borda(indiv_ranks)[[1]]
    
    topk %>%
      mutate(lineage=lineage,cell_ID=cell, rank=row_number()) %>%
      pivot_longer(cols = c(mean,median,geo.mean,l2norm), names_to = "algorithm", values_to = "driver") %>%
      mutate(algorithm = paste0("consensus_topklists_",gsub("\\.","",algorithm))) %>%
      dplyr::select(lineage,cell_ID,driver,rank,algorithm)
    
  }
  
  end <- Sys.time()
  
  message("Rank aggregation took ", round(end-start,2), " seconds")
  
}

suppressWarnings(dir.create(paste0("results/CCLE_",network_choice,"/consensus")))

write_csv(consensus_drivers,paste0("results/CCLE_",network_choice,"/consensus/",RAmethod,".csv"))



if(RAmethod=="BiG"){
  
  algs <- sort(unique(aggregated_results$algorithm))
  
  # For each cell, get a list of the ranked drivers from each algorithm
  
  start <- Sys.time()
  
  consensus_drivers <- foreach(cell=unique(aggregated_results$cell_ID), .packages = c("tidyverse","TopKLists","foreach"), .combine = "rbind") %dopar% {
    
    #message(paste0("Getting consensus drivers for ",cell))
    
    indiv_ranks <- aggregated_results %>%
      filter(cell_ID==cell) %>%
      pivot_wider(names_from = "algorithm", values_from = "driver", values_fill = NA) %>%
      arrange(rank) %>%
      dplyr::select(-c(cell_ID,lineage,rank)) %>%
      as.matrix()
      
    
    lineage <- sample_info %>% filter(cell_ID==cell) %>% pull(lineage)
    
    lengths <- colSums(!is.na(indiv_ranks))
    
    BiG <- BiG_diffuse(r=indiv_ranks,n_T = lengths,n_p1=20000, M=2000, burnin=1000, prior="IG")
    
    topk %>%
      mutate(lineage=lineage,cell_ID=cell, rank=row_number()) %>%
      pivot_longer(cols = c(mean,median,geo.mean,l2norm), names_to = "algorithm", values_to = "driver") %>%
      mutate(algorithm = paste0("consensus_topklists_",gsub("\\.","",algorithm))) %>%
      dplyr::select(lineage,cell_ID,driver,rank,algorithm)
    
  }
  
  end <- Sys.time()
  
  message("Rank aggregation took ", round(end-start,2), " seconds")
  
}



