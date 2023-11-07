#evaluate_compare_outputs.r

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(MASS, quietly = T))
suppressPackageStartupMessages (library(ggrepel, quietly = T))
suppressPackageStartupMessages (library(foreach, quietly = T))
suppressPackageStartupMessages (library(doParallel, quietly = T))


# Handling input arguments
option_list = list(
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network to use", metavar ="Network"),
  make_option(c("-a", "--algorithms"), type="character", default="ALL", 
              help="algorithms to include in comparison separated by semicolons (;), or 'ALL' (Default), or lead with '-' to exclude", metavar ="Algorithms"),
  make_option(c("-t", "--threads"), type="integer", default=4, 
              help="Number of threads to use if running in parallel", metavar ="Rank Aggregation Method"),
  make_option(c("-p", "--plot"), type="character", default="-consensus_BiG;-consensus_topklists_mean;-consensus_topklists_median;-consensus_topklists_l2norm", 
              help="Algorithms to include in the final plot, separated by semicolons (;), or 'ALL' (Default), or lead with '-' to exclude ", metavar ="Rank Aggregation Method")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

network_choice <- opt$network
algorithms <- opt$algorithms
threads <- opt$threads
plot_algorithms <- opt$plot

algorithms <- str_split(algorithms, ";") %>% unlist()

if(toupper(plot_algorithms[1]) != "ALL"){
  plot_algorithms <- str_split(plot_algorithms, ";") %>% unlist()
  plot_keep <- plot_algorithms[str_detect(plot_algorithms,"^[^-]")]
  plot_exclude <- plot_algorithms[str_detect(plot_algorithms,"^-")]
  plot_exclude <- gsub("^-","",plot_exclude)
}



if(length(plot_keep)>0 & length(plot_exclude)>0){
  warning("Warning: Can't combine included and excluded algorithms for -p option")
}

threads <- 10

if(threads>1){
  #registerDoParallel(cores=threads)
  cl <- makeCluster(threads, outfile = "log/evaluate_similarity.log")
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
#algorithms <- "ALL"

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

#aggregated_results <- aggregated_results %>% filter(lineage=="Liver")

#############################
# Similarity of Results
#############################
#top_ns <- c(1,5,10,50,100,150,200)
#top_ns <- append(c(1,5),seq(10,max(aggregated_results$rank), by = 10))
top_ns <- c(1,2,3,4,5,6,7,8,9,10,50,100)
alg_pairs <- as.list(as.data.frame(combn(unique(aggregated_results$algorithm), 2)))

suppressWarnings(rm(combined_cos_sim))


if(!file.exists(paste0("results/CCLE_",network_choice,"/full_results_similarity.csv"))){
  

combined_sim <- foreach(alg_pair=alg_pairs, .combine = "rbind",  .packages = c("tidyverse","foreach"), .export = c("aggregated_results", "top_ns")) %dopar% {
  
  alg1 <- alg_pair[1]
  alg2 <- alg_pair[2]
  
  tmp_results_1 <- aggregated_results %>% filter(algorithm %in% c(alg1, alg2))
  
  comp_cells <- intersect(tmp_results_1 %>% filter(algorithm==alg1) %>% pull(cell_ID),
                          tmp_results_1 %>% filter(algorithm==alg2) %>% pull(cell_ID)
  )
  
  message(paste("Calculating cosine and jaccard similarity for", alg1, "vs", alg2, sep = " "))
  
  
  foreach(cell=comp_cells, .combine = "rbind", .packages = c("tidyverse","foreach"), .export = c("tmp_results_1","cell","alg1","alg2")) %do% {
    
    tmp_results_2 <- tmp_results_1 %>% filter(cell_ID == cell)
    
    foreach(n=top_ns, .combine = "rbind", .packages = c("tidyverse","foreach"), .export = c("tmp_results_2","cell","alg1","alg2")) %do% {
      
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
      
      
      
      data.frame(cell_ID = cell,
                 algorithm_1 = alg1, 
                 algorithm_2 = alg2,
                 top_n = n,
                 cos_sim = cos_sim,
                 jac_sim = jac_sim
      )
    }
  }
  
  
}

write_csv(combined_sim, paste0("results/CCLE_",network_choice,"/full_results_similarity.csv"))
}else{
  message("Similarity stats file already exists. Reading in previous results.")
  message(paste0("To stop this, remove file: ", "results/CCLE_",network_choice,"/full_results_similarity.csv"))
  combined_sim <- read_csv(paste0("results/CCLE_",network_choice,"/full_results_similarity.csv"))
}


#############################
# Cosine Similarity
#############################

summarised_cos_sim <- combined_sim %>%
  group_by(algorithm_1,algorithm_2,top_n) %>%
  summarise(mean_cos_sim = mean(cos_sim, na.rm = T)) %>%
  mutate(mean_cos_dist = 1-mean_cos_sim)

if(toupper(plot_algorithms[1]) != "ALL"){
  if(length(plot_keep)>0){
    summarised_cos_sim <- summarised_cos_sim %>% filter(algorithm_1 %in% plot_keep & algorithm_2 %in% plot_keep)
  }else if(length(plot_exclude)>0){
    summarised_cos_sim <- summarised_cos_sim %>% filter(!algorithm_1 %in% plot_exclude & !algorithm_2 %in% plot_exclude)
  }
    
  
}

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
    unique() %>%
    # Make distance matrix
    pivot_wider(names_from = algorithm_2, values_from = mean_cos_dist, values_fill = 0) %>%
    column_to_rownames("algorithm_1") %>%
    as.matrix()
  dist_matrix <- dist_matrix[sort(unique(append(summarised_cos_sim$algorithm_1,summarised_cos_sim$algorithm_2))),sort(unique(append(summarised_cos_sim$algorithm_1,summarised_cos_sim$algorithm_2)))]
  
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

alg_colours <- read_csv("data/algorithm_colours.csv") %>% deframe()

cos_similarity_plot <- cos_similarity_plot %>%
  mutate(algorithm=factor(algorithm, levels = names(alg_colours)))




ggplot(cos_similarity_plot %>% mutate(start_label = if_else(n_drivers == min(n_drivers), as.character(algorithm), NA_character_),
                                  end_label = if_else(n_drivers == max(n_drivers), as.character(algorithm), NA_character_)), 
       aes(x=as.factor(n_drivers),
           y=pos, 
           colour=algorithm, 
           group=algorithm
           )) +
  #geom_point() +
  geom_line() +
  scale_colour_manual(values = alg_colours) +
  ylab("Cosine Similarity (Dimensionality Reduced)") +
  xlab("Top N Drivers") +
  geom_label_repel(aes(label=start_label),
                   nudge_x = -5,
                   na.rm = T, direction = "y", arrow = arrow(type = "open", angle = 10, length = unit(0.3,"cm")),
                   max.overlaps = 100) +
  #geom_label_repel(aes(label=end_label),
  #                 nudge_x = 5,
  #                 na.rm = T) +
  theme_classic() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #legend.position = "none",
        #panel.background = element_rect(fill="#1B2631")
        ) +
  guides(colour="none") +
  geom_vline(xintercept = 1, colour = "black", linetype = "dashed")

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
