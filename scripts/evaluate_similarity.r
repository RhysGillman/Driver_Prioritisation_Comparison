#evaluate_compare_outputs.r

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(MASS, quietly = T))
suppressPackageStartupMessages (library(ggrepel, quietly = T))


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
algorithms <- "ALL"

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
# Similarity of Results
#############################
#top_ns <- c(1,5,10,50,100,150,200)
#top_ns <- append(c(1,5),seq(10,max(aggregated_results$rank), by = 10))
top_ns <- c(1,2,3,4,5,6,7,8,9,10,50,100)
alg_pairs <- as.list(as.data.frame(combn(unique(aggregated_results$algorithm), 2)))

suppressWarnings(rm(combined_cos_sim))
for(alg_pair in alg_pairs){
  
  alg1 <- alg_pair[1]
  alg2 <- alg_pair[2]
  
  tmp_results_1 <- aggregated_results %>% filter(algorithm %in% c(alg1, alg2))
  
  comp_cells <- intersect(tmp_results_1 %>% filter(algorithm==alg1) %>% pull(cell_ID),
                          tmp_results_1 %>% filter(algorithm==alg2) %>% pull(cell_ID)
  )
  
  message(paste("Calculating cosine and jaccard similarity for", alg1, "vs", alg2, sep = " "))
  
  for(cell in comp_cells){
    
    tmp_results_2 <- tmp_results_1 %>% filter(cell_ID == cell)
    
    
    
    for(n in top_ns){
      
      # Calculate the cosine similarity of results from each cell
      
      cos_sim_matrix <- tmp_results_2 %>% filter(rank <= n) %>% 
        ungroup() %>% 
        dplyr::select(driver, algorithm) %>%
        table() %>%
        as.data.frame() %>%
        pivot_wider(names_from = algorithm, values_from = Freq)
      
      cos_sim <- (
        # This calculates the sum of the products of the two vectors
        as.numeric(t(pull(cos_sim_matrix,2))%*%pull(cos_sim_matrix,3))
        # Divide
      ) /
        # This multiplies the square root of the sum squares of each vector
        (sqrt(sum(pull(cos_sim_matrix,2)^2))*sqrt(sum(pull(cos_sim_matrix,3)^2)))
      
      
      # Calculate the jaccard similarity of results from each cell
      
      jac_sim <- (
        as.numeric(length(intersect(tmp_results_2 %>% filter(rank <= n) %>% filter(algorithm==alg1) %>% pull(driver),
                                    tmp_results_2 %>% filter(rank <= n) %>% filter(algorithm==alg2) %>% pull(driver)))) / 
        as.numeric(length(union(tmp_results_2 %>% filter(rank <= n) %>% filter(algorithm==alg1) %>% pull(driver),
                                tmp_results_2 %>% filter(rank <= n) %>% filter(algorithm==alg2) %>% pull(driver))))
      )
      
      
      
      indiv_result <- data.frame(cell_ID = cell, 
                                 algorithm_1 = alg1, 
                                 algorithm_2 = alg2,
                                 top_n = n,
                                 cos_sim = cos_sim,
                                 jac_sim = jac_sim
      )
      
      if(!exists("combined_sim")){
        combined_sim <- indiv_result
      }else{
        combined_sim <- rbind(combined_sim, indiv_result)
        
      }
      
    }
  }
}

#############################
# Cosine Similarity
#############################

summarised_cos_sim <- combined_sim %>%
  group_by(algorithm_1,algorithm_2,top_n) %>%
  summarise(mean_cos_sim = mean(cos_sim, na.rm = T)) %>%
  mutate(mean_cos_dist = 1-mean_cos_sim)

#ggplot(summarised_cos_sim, aes(x=top_n,y=mean_cos_dist)) +
#  geom_point()


suppressWarnings(rm(cos_similarity_plot))
for(n in top_ns){
  dist_matrix <- summarised_cos_sim %>%
    filter(top_n==n) %>%
    dplyr::select(algorithm_1,algorithm_2,mean_cos_dist)
  dist_matrix <- dist_matrix %>%
    # Need to duplicate values to make full symmetric matrix
    rbind(dist_matrix %>% dplyr::select(algorithm_1 = algorithm_2, algorithm_2 = algorithm_1, mean_cos_dist)) %>%
    # Make distance matrix
    pivot_wider(names_from = algorithm_2, values_from = mean_cos_dist, values_fill = 0) %>%
    column_to_rownames("algorithm_1") %>%
    as.matrix()
  dist_matrix <- dist_matrix[sort(unique(aggregated_results$algorithm)),sort(unique(aggregated_results$algorithm))]
  
  mds <- isoMDS(dist_matrix, k = 1, tol = 0.1, p = 2)
  
  message(paste0("model stress for n = ",n," is ", mds$stress))
  
  tmp_similarity_plot <- mds$points %>%
    as.data.frame() %>%
    dplyr::rename(pos=V1) %>%
    rownames_to_column("algorithm") %>%
    mutate(n_drivers = n)
  
  if(!exists("cos_similarity_plot")){
    cos_similarity_plot <- tmp_similarity_plot
  }else{
    cos_similarity_plot <- rbind(cos_similarity_plot, tmp_similarity_plot)
    
  }
  
}

# Invert values if order of algorithms becomes reversed
for(n in top_ns){
  
  if(n!=top_ns[1]){
    current_order <- cos_similarity_plot %>% 
      filter(n_drivers==n) %>%
      arrange(desc(pos)) %>%
      pull(algorithm)
    
    current_match <- sum(current_order==prev_order)
    reverse_match <- sum(current_order==rev(prev_order))
    
    #print(paste0(current_match," vs ",reverse_match," reverse = ", reverse_match>current_match))
    
    if(reverse_match>current_match){
      cos_similarity_plot <- cos_similarity_plot %>%
        mutate(pos=ifelse(n_drivers==n,pos/-1,pos))
    }
    
    
  }
  
  prev_order <- cos_similarity_plot %>% 
    filter(n_drivers==n) %>%
    arrange(desc(pos)) %>%
    pull(algorithm)
  
}

#cos_similarity_plot <- cos_similarity_plot %>%
#  mutate(log_pos = log10(pos + abs(min(cos_similarity_plot$pos))+1))



ggplot(cos_similarity_plot %>% mutate(start_label = if_else(n_drivers == min(n_drivers), as.character(algorithm), NA_character_),
                                  end_label = if_else(n_drivers == max(n_drivers), as.character(algorithm), NA_character_)), 
       aes(x=as.factor(n_drivers),
           y=pos, colour=algorithm, group=algorithm)) +
  geom_point() +
  geom_line() +
  ylab("Cosine Similarity (Dimensionality Reduced)") +
  xlab("Top N Drivers") +
  geom_label_repel(aes(label=start_label),
                   nudge_x = -5,
                   na.rm = T) +
  geom_label_repel(aes(label=end_label),
                   nudge_x = 5,
                   na.rm = T) +
  theme_classic() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")

ggsave(paste0("results/CCLE_",network_choice,"/compare_results_cosine.png"), width = 5000, height = 3000 ,units = "px")

#############################
# Jaccard Distance
#############################

summarised_jac_sim <- combined_sim %>%
  group_by(algorithm_1,algorithm_2,top_n) %>%
  summarise(mean_jac_sim = mean(jac_sim, na.rm = T)) %>%
  mutate(mean_jac_dist = 1-mean_jac_sim)

#ggplot(summarised_cos_sim, aes(x=top_n,y=mean_cos_dist)) +
#  geom_point()


suppressWarnings(rm(jac_similarity_plot))
for(n in top_ns){
  dist_matrix <- summarised_jac_sim %>%
    filter(top_n==n) %>%
    dplyr::select(algorithm_1,algorithm_2,mean_jac_dist)
  dist_matrix <- dist_matrix %>%
    # Need to duplicate values to make full symmetric matrix
    rbind(dist_matrix %>% dplyr::select(algorithm_1 = algorithm_2, algorithm_2 = algorithm_1, mean_jac_dist)) %>%
    # Make distance matrix
    pivot_wider(names_from = algorithm_2, values_from = mean_jac_dist, values_fill = 0) %>%
    column_to_rownames("algorithm_1") %>%
    as.matrix()
  dist_matrix <- dist_matrix[sort(unique(aggregated_results$algorithm)),sort(unique(aggregated_results$algorithm))]
  
  mds <- isoMDS(dist_matrix, k = 1, tol = 0.1)
  
  message(paste0("model stress for n = ",n," is ", mds$stress))
  
  tmp_similarity_plot <- mds$points %>%
    as.data.frame() %>%
    dplyr::rename(pos=V1) %>%
    rownames_to_column("algorithm") %>%
    mutate(n_drivers = n)
  
  if(!exists("jac_similarity_plot")){
    jac_similarity_plot <- tmp_similarity_plot
  }else{
    jac_similarity_plot <- rbind(jac_similarity_plot, tmp_similarity_plot)
    
  }
  
}

# Invert values if order of algorithms becomes reversed
for(n in top_ns){
  
  if(n!=top_ns[1]){
    current_order <- jac_similarity_plot %>% 
      filter(n_drivers==n) %>%
      arrange(desc(pos)) %>%
      pull(algorithm)
    
    current_match <- sum(current_order==prev_order)
    reverse_match <- sum(current_order==rev(prev_order))
    
    #print(paste0(current_match," vs ",reverse_match," reverse = ", reverse_match>current_match))
    
    if(reverse_match>current_match){
      jac_similarity_plot <- jac_similarity_plot %>%
        mutate(pos=ifelse(n_drivers==n,pos/-1,pos))
    }
    
    
  }
  
  prev_order <- jac_similarity_plot %>% 
    filter(n_drivers==n) %>%
    arrange(desc(pos)) %>%
    pull(algorithm)
  
}



ggplot(jac_similarity_plot %>% mutate(start_label = if_else(n_drivers == min(n_drivers), as.character(algorithm), NA_character_),
                                  end_label = if_else(n_drivers == max(n_drivers), as.character(algorithm), NA_character_)), 
       aes(x=as.factor(n_drivers),
           y=pos, colour=algorithm, group=algorithm)) +
  geom_point() +
  geom_line() +
  ylab("Jaccard Similarity (Dimensionality Reduced)") +
  xlab("Top N Drivers") +
  geom_label_repel(aes(label=start_label),
                   nudge_x = -5,
                   na.rm = T) +
  geom_label_repel(aes(label=end_label),
                   nudge_x = 5,
                   na.rm = T) +
  theme_classic() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")

ggsave(paste0("results/CCLE_",network_choice,"/compare_results_jaccard.png"), width = 5000, height = 3000 ,units = "px")
