# find_best_consensus.r

# evaluate_SL_partners.r

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(foreach, quietly = T))
suppressPackageStartupMessages (library(doParallel, quietly = T))
suppressPackageStartupMessages (library(TopKLists, quietly = T))






# Handling input arguments
option_list = list(
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network being used", metavar ="Network"),
  make_option(c("-a", "--algorithms"), type="character", default="ALL", 
              help="algorithms to include in comparison separated by spaces, or 'ALL' (Default), or lead with '-' to exclude", metavar ="Algorithms"),
  make_option(c("-t", "--threads"), type="integer", default=4, 
              help="Number of threads to use if running in parallel", metavar ="Threads"),
  make_option(c("-c", "--consensusmethod"), type="character", default="topklists", 
              help="Consensus method to use", metavar ="Consensus Method"),
  make_option(c("-s", "--samplesize"), type="integer", default=50, 
              help="Number of randomly selected samples to use", metavar ="Sample Size")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


network_choice <- opt$network
algorithms <- opt$algorithms
threads <- opt$threads
consensus_method <- opt$consensusmethod
sample_size <- opt$samplesize

threads = 10


if(threads>1){
  cl <- makeCluster(threads, outfile = "log/find_best_consensus.log")
  registerDoParallel(cl)
}



# Best De Novo Method
SL_level <- fread("results/CCLE_STRINGv11/full_results_SL_partner_level.csv") %>% 
  filter(n==20) %>%
  # Only keeping results for samples with > 10 gold standards
  filter(length_gs >= 10) %>%
  group_by(algorithm,n) %>%
  # Only keep measurements where more than 10 cells are available to calculate mean
  filter(n()>=10) %>%
  summarise(
    mean_TP = mean(TP),
    mean_precision = mean(precision),
    mean_recall = mean(recall),
    mean_F1 = median(F1)
  )

ref_level <- fread("results/CCLE_STRINGv11/full_results_CGC_all_reference_driver_level.csv") %>% 
  filter(n==20) %>%
  # Only keeping results for samples with > 10 gold standards
  filter(n_gs >= 10) %>%
  group_by(algorithm,n) %>%
  # Only keep measurements where more than 10 cells are available to calculate mean
  filter(n()>=10) %>%
  summarise(
    mean_TP = mean(TP),
    mean_precision = mean(precision),
    mean_recall = mean(recall),
    mean_F1 = median(F1)
  )


# Top Picks (hand chosen based on results above)

algorithms <- c(
  "CSN_NCUA",
  "PersonaDrive",
  "SCS",
  "DawnRank",
  "PRODIGY",
  "OncoImpact",
  "sysSVM2",
  "PNC",
  "PhenoDriverR"
)

topn <- c(5,10,50)


#############################
# Sample Info
#############################

sample_info <- read_csv(paste0("validation_data/CCLE_",network_choice,"/sample_info.csv"))
samples <- sample_info$cell_ID %>% sort()

#############################
# Data
#############################

# Gold Standards

gold_standard <- read_csv(paste0("validation_data/CCLE_",network_choice,"/gold_standards.csv")) %>%
  dplyr::select(cell_ID,gene_ID) %>%
  unique() %>%
  group_by(cell_ID) %>%
  summarise(sensitive_genes = list(gene_ID)) %>%
  deframe()


rare_gold_standard <- read_csv(paste0("validation_data/CCLE_",network_choice,"/rare_gold_standards.csv")) %>%
  dplyr::select(cell_ID,gene_ID) %>%
  unique() %>%
  group_by(cell_ID) %>%
  summarise(sensitive_genes = list(gene_ID)) %>%
  deframe()

# SL
max_SL=5

SL <- read_csv(paste0("validation_data/CCLE_",network_choice,"/SL_partners_max_",max_SL,".csv")) %>%
  dplyr::select(gene_ID=gene1,partner=gene2,SL_rank=rank)


#############################
# Read In LOF/GOF Predictions
#############################

LOF_GOF <- read_csv("data/LOF_GOF_annotations.csv")





#############################
# Read In Results
#############################

de_novo_collection <- c("CSN_DFVS",     "CSN_MDS",     "CSN_MMS",      "CSN_NCUA",
                        "LIONESS_DFVS", "LIONESS_MDS", "LIONESS_MMS",  "LIONESS_NCUA",
                        "SPCC_DFVS",    "SPCC_MDS",    "SPCC_MMS",     "SPCC_NCUA",
                        "SSN_DFVS",     "SSN_MDS",     "SSN_MMS",      "SSN_NCUA")

suppressWarnings(rm(aggregated_results,specific_de_novo))
for(alg in list.dirs(paste0("results/CCLE_",network_choice), recursive = F)){
  # First, get names of algorithms that have been run
  alg = str_extract(alg,"(?<=/)[^/]*$")
  # Skip consensus results. These will be added separately
  if(alg=="consensus"){next}
  
  
  # Check which algorithms are to be included
  
  # First check whether all are to be included, if so, then skip ahead to reading in files
  if(algorithms[1]!="ALL"){
    
    # Check if all of the names begin with "-" to indicate they are to be excluded
    if(all(substr(algorithms,1,1)=="-")){
      
      # If so, then skip any algorithms listed
      if(alg %in% gsub("-","",algorithms)){next}
      
    }else if(any(substr(algorithms,1,1)=="-")){
      warning("Invalid input for --algorithms. Cannot mix inclusions and exclusions in one statement.")
      stop() 
      # If the requested algorithm isn't present then skip it
    }else if(!alg %in% algorithms){
      
      # Unless it is a specific de novo method
      if(any(algorithms %in% de_novo_collection)){
        specific_de_novo <- algorithms[algorithms %in% de_novo_collection]
        alg <- "combined_de_novo_methods"
      }else{
        next
      }
    }
  }
  for(result_file in list.files(paste0("results/CCLE_", network_choice,"/",alg), pattern = "*.csv")){
    message(paste0("Reading result file for ", alg, "-", result_file))
    indiv_result <- data.table::fread(paste0("results/CCLE_", network_choice,"/",alg,"/",result_file))
    if(!"algorithm" %in% colnames(indiv_result)){
      indiv_result <- indiv_result %>% mutate(algorithm = alg)
    }
    if(alg == "combined_de_novo_methods" & exists("specific_de_novo")){
      indiv_result <- indiv_result %>% filter(algorithm %in% specific_de_novo)
    }
    if(!exists("aggregated_results")){
      aggregated_results <- indiv_result
    }else{
      aggregated_results <- rbind(aggregated_results,indiv_result)
    }
  }
}

#############################
# Downsampling
#############################
# set seed for reproducible results
set.seed(1000)

# Select 10 random cells from (sample size/10) random lineages
selected_lineages <- sample(unique(sample_info$lineage),sample_size/10, replace = FALSE)

# Only sample from cells with gold standards available and of a decent length
selected_samples <- foreach(lin=selected_lineages, .combine = "c") %do% {
  sample(intersect(sample_info %>% filter(lineage==lin) %>% pull(cell_ID) %>% unique(), names(which(lapply(gold_standard, length)>10))), 10, replace = FALSE)
  }


results_sample <- aggregated_results %>%
  filter(cell_ID %in% selected_samples)


#############################
# Finding All Combinations
#############################

if(algorithms[1]=="ALL"){
  algs <- unique(results_sample$algorithm)
}else{
  algs <- algorithms
}

for(i in seq(2,length(algs))){
  indiv_comb <- as.list(as.data.frame(combn(algs,i)))
  if(!exists("alg_combinations")){
    alg_combinations <- indiv_comb
  }else{
    alg_combinations <- append(alg_combinations,indiv_comb)
  }
}

#############################
# Rank Aggregation
#############################

if(consensus_method=="topklists"){
  
  if(!file.exists(paste0("results/CCLE_",network_choice,"/sampled_topklist_consensus_combinations.csv"))){

  consensus_drivers <- foreach(cell=selected_samples, .packages = c("tidyverse","TopKLists","foreach"), .combine = "rbind") %dopar% {
    
    message(paste0("Getting consensus drivers for ",cell, " (", which(selected_samples==cell),"/",length(selected_samples),")"))
    
    alg_set_results <- foreach(index=1:length(alg_combinations), .combine = "rbind") %do% {
      
      if((index/length(alg_combinations)) %in% c(round(length(alg_combinations)/4,0)/length(alg_combinations), #~ 25%
                                                 round(length(alg_combinations)/2,0)/length(alg_combinations), #~ 50%
                                                 round(length(alg_combinations)*3/4,0)/length(alg_combinations), #~ 75%
                                                 1) ){
        message(paste0(round(index/length(alg_combinations)*100,0),"%"))
      }
      
      alg_set <- alg_combinations[[index]]
      
      indiv_ranks <- foreach(alg=alg_set) %do% {
      
      x <- results_sample %>% 
        filter(cell_ID==cell) %>%
        filter(algorithm==alg) %>%
        pull(driver)
      
      # NA indicates there were no results from that algorithm for that cell_ID. In this case, replace with empty vector
      if(any(is.na(x))){
        x <- c()
      }
      
      x
      
      }
      
      names(indiv_ranks) <- alg_set
      
      topk <- TopKLists::Borda(indiv_ranks)[[1]]
      
      topk %>%
        mutate(consensus_algorithms=paste(alg_set, collapse = ","),cell_ID=cell, rank=row_number()) %>%
        pivot_longer(cols = c(mean,median,geo.mean,l2norm), names_to = "type", values_to = "driver") %>%
        dplyr::select(consensus_algorithms,type,cell_ID,driver,rank)
      
      
    }
    
    
    
  }
  
  write_csv(consensus_drivers, paste0("results/CCLE_",network_choice,"/sampled_topklist_consensus_combinations.csv"))
  
  }else{
  
    consensus_drivers <- fread(paste0("results/CCLE_",network_choice,"/sampled_topklist_consensus_combinations.csv"))
}


}



#############################
# Evaluate SL Partners
#############################



# Driver-SL-Partner-Level Evaluation


SL_partner_results <- consensus_drivers %>%
  mutate(algorithm=paste0(consensus_algorithms,"_",type)) %>%
  dplyr::select(algorithm, cell_ID, driver, rank) %>%
  dplyr::rename(driver_rank=rank) %>%
  # Have to filter results down to help compuation
  filter(driver_rank <= 100) %>%
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
  arrange(algorithm, cell_ID, final_rank)



tmp_samples <- intersect(SL_partner_results$cell_ID, names(gold_standard))

if(!file.exists(paste0("results/CCLE_",network_choice,"/sampled_topklist_consensus_combinations_SL_results.csv"))){
  
  SL_prediction_stats <- foreach(sample=tmp_samples, .combine = "rbind", .packages = c("tidyverse","foreach")) %dopar% {
    print(paste0("Calculating stats for ", sample, "(",which(tmp_samples == sample),"/",length(tmp_samples),")"))
    gs <- gold_standard[sample] %>% unlist()
    tmp1 <- SL_partner_results %>% filter(cell_ID == sample)
    
    foreach(alg=unique(tmp1$algorithm), .combine = "rbind", .packages = c("tidyverse","foreach")) %do% {
      
      tmp2 <- tmp1 %>% filter(algorithm == alg)
      if(nrow(tmp2)==0){break}
      
      foreach(n=topn, .combine = "rbind", .packages = c("tidyverse")) %do% {
      
        
      predicted <- tmp2 %>% filter(final_rank <= n) %>% pull(target)
      correct <- predicted[which(predicted %in% gs)]
      wrong <- predicted[which(!predicted %in% gs)]
      TP <- length(correct)
      FP <- length(wrong)
      precision <- TP/n
      recall <- TP/length(gs)
      F1 <- 2*((precision*recall)/(precision + recall))
      if(is.nan(F1)){
        F1 <- 0
      }
      
      data.frame(cell_ID = sample, 
                 algorithm = alg, 
                 n = n,
                 correct = paste0(correct, collapse = ";"), 
                 TP = TP,
                 FP = FP,
                 precision = precision,
                 recall = recall,
                 F1 = F1,
                 length_gs = length(gs))
        
        
      }
      
      
      
    }
    
  }
  write_csv(SL_prediction_stats, paste0("results/CCLE_",network_choice,"/sampled_topklist_consensus_combinations_SL_results.csv"))
}else{
  message("SL Prediction stats file already exists. Reading in previous results.")
  message(paste0("To stop this, remove file: ", "results/CCLE_",network_choice,"/sampled_topklist_consensus_combinations_SL_results.csv"))
  SL_prediction_stats <- read_csv(paste0("results/CCLE_",network_choice,"/sampled_topklist_consensus_combinations_SL_results.csv"))
}


SL_prediction_stats_summarised <- SL_prediction_stats %>%
  group_by(algorithm,n) %>%
  summarise(
    mean_TP = mean(TP),
    mean_precision = mean(precision),
    mean_recall = mean(recall),
    mean_F1 = median(F1),
    sample_size = n()
  ) %>%
  pivot_longer(cols = -c(algorithm,n,sample_size), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure, levels = c("mean_TP","mean_recall","mean_precision","mean_F1")))


# Find best overall

suppressWarnings(rm(ranked_methods))
for(index in topn){
  for(type in unique(SL_prediction_stats_summarised$measure)){
    rank <- SL_prediction_stats_summarised %>% filter(n==index, measure==type) %>% arrange(desc(value)) %>% pull(algorithm)
    rank <- list(rank)
    names(rank) <- paste0(type,"_",index)
    if(!exists("ranked_methods")){
      ranked_methods <- rank
    }else{
      ranked_methods <- append(ranked_methods,rank)
    }
  }
}



top_methods <- TopKLists::Borda(ranked_methods)[[1]] %>%
  mutate(rank=row_number()) %>%
  pivot_longer(cols = c(mean,median,geo.mean,l2norm), names_to = "measure", values_to = "consensus_algorithm") %>%
  filter(measure=="geo.mean") %>%
  arrange(rank)

#top_method <- "CSN_NCUA,sysSVM2,PhenoDriverR_mean"
#top_method <- "SCS,sysSVM2,PNC,PhenoDriverR_median"
#top_method <- "CSN_NCUA,OncoImpact,sysSVM2,PhenoDriverR_geo.mean"
top_method <- "CSN_NCUA,OncoImpact,sysSVM2_geo.mean"
#top_method <- top_methods %>% head(1) %>% pull(consensus_algorithm)

#alg_colours <- c("red",rep("black",length(unique(SL_prediction_stats_summarised$algorithm))-1))
#names(alg_colours) <- c(top_method, unique(SL_prediction_stats_summarised$algorithm)[which(unique(SL_prediction_stats_summarised$algorithm)!=top_method)] )

#alg_widths <- alg_colours
#alg_widths[alg_widths=="red"] <- 2
#alg_widths[alg_widths=="black"] <- 0.5
#alg_widths <- as.numeric(alg_widths)
#names(alg_widths) <- names(alg_colours)

plot_data <- SL_prediction_stats_summarised %>%
  mutate(plot_colour = ifelse(algorithm==top_method, "red", "black"),
         plot_linewidth = ifelse(algorithm==top_method, 2, 0.5),
         plot_linetype = ifelse(algorithm==top_method, "dashed", "solid"))


ggplot(mapping=aes(x = n, y = value, color = algorithm, linewidth=algorithm, linetype=algorithm)) +
  geom_line(data=plot_data %>% filter(algorithm != top_method)) +
  geom_line(data=plot_data %>% filter(algorithm == top_method)) +
  scale_color_manual(breaks = plot_data$algorithm, values = plot_data$plot_colour) +
  scale_linewidth_manual(breaks = plot_data$algorithm, values = plot_data$plot_linewidth) +
  scale_linetype_manual(breaks = plot_data$algorithm, values = plot_data$plot_linetype) +
  ylab("Average Value") +
  xlab("Number of Predicted Sensitive Genes") +
  ggtitle(top_method) +
  theme_light() +
  theme(legend.position = "none") +
  facet_wrap(~measure, nrow = 1, scales = "free")

ggsave(paste0("results/CCLE_",network_choice,"/choosing_top_consensus_method_SL_partner_level.png"),width = 40, height = 30, units = "cm", dpi = 300)
  
















########################
# Rare SL Partner Level
########################


global_freq_thresh <- 0.5
lineage_freq_thresh <- 0.25


rare_SL_partner_predictions <- SL_partner_results %>%
  left_join(sample_info %>% dplyr::select(cell_ID,lineage), by = "cell_ID") %>%
  dplyr::select(lineage, cell_ID, algorithm, target, final_rank) %>%
  ungroup() %>%
  # Because some algorithms give very long lists, they will be greatly affected by trying to find only rare drivers
  # need to filter to only high ranks first
  filter(final_rank <= 100) %>%
  # Get the total # of cells
  mutate(total_n = length(unique(cell_ID))) %>%
  # Get the global frequency of predicted essential gene for each algorithm
  group_by(algorithm,target) %>% mutate(global_essentiality_count = length(unique(cell_ID))) %>% ungroup() %>%
  # Get the total # of cells per lineage
  group_by(lineage) %>% mutate(lineage_n = length(unique(cell_ID))) %>% ungroup() %>%
  # Get the lineage frequency of essentiality
  group_by(algorithm,lineage,target) %>% mutate(lineage_essentiality_count = length(unique(cell_ID))) %>% ungroup() %>%
  mutate(global_essentiality_frequency = global_essentiality_count / total_n,
         lineage_essentiality_frequency = lineage_essentiality_count / lineage_n) %>%
  # Filter frequency thresholds
  filter(global_essentiality_frequency < global_freq_thresh, lineage_essentiality_frequency < lineage_freq_thresh) %>%
  dplyr::select(lineage, cell_ID, algorithm, target, final_rank, global_essentiality_frequency, lineage_essentiality_frequency) %>%
  # update the ranks
  group_by(algorithm,cell_ID) %>%
  arrange(final_rank) %>%
  mutate(final_rank = row_number()) %>%
  ungroup() %>%
  arrange(algorithm,lineage,cell_ID,final_rank)


tmp_samples <- intersect(rare_SL_partner_predictions$cell_ID, names(rare_gold_standard))

if(!file.exists(paste0("results/CCLE_",network_choice,"/sampled_topklist_consensus_combinations_rare_SL_results.csv"))){
  
  rare_SL_prediction_stats <- foreach(sample=tmp_samples, .combine = "rbind", .packages = c("tidyverse","foreach")) %dopar% {
    print(paste0("Calculating stats for ", sample, "(",which(tmp_samples == sample),"/",length(tmp_samples),")"))
    gs <- rare_gold_standard[sample] %>% unlist()
    tmp1 <- rare_SL_partner_predictions %>% filter(cell_ID == sample)
    
    foreach(alg=unique(tmp1$algorithm), .combine = "rbind", .packages = c("tidyverse","foreach")) %do% {
      
      tmp2 <- tmp1 %>% filter(algorithm == alg)
      if(nrow(tmp2)==0){break}
      
      foreach(n=topn, .combine = "rbind", .packages = c("tidyverse")) %do% {
        
        
        predicted <- tmp2 %>% filter(final_rank <= n) %>% pull(target)
        correct <- predicted[which(predicted %in% gs)]
        wrong <- predicted[which(!predicted %in% gs)]
        TP <- length(correct)
        FP <- length(wrong)
        precision <- TP/n
        recall <- TP/length(gs)
        F1 <- 2*((precision*recall)/(precision + recall))
        if(is.nan(F1)){
          F1 <- 0
        }
        
        data.frame(cell_ID = sample, 
                   algorithm = alg, 
                   n = n,
                   correct = paste0(correct, collapse = ";"), 
                   TP = TP,
                   FP = FP,
                   precision = precision,
                   recall = recall,
                   F1 = F1,
                   length_gs = length(gs))
        
        
      }
      
      
      
    }
    
  }
  write_csv(rare_SL_prediction_stats, paste0("results/CCLE_",network_choice,"/sampled_topklist_consensus_combinations_rare_SL_results.csv"))
}else{
  message("SL Prediction stats file already exists. Reading in previous results.")
  message(paste0("To stop this, remove file: ", "results/CCLE_",network_choice,"/sampled_topklist_consensus_combinations_rare_SL_results.csv"))
  rare_SL_prediction_stats <- read_csv(paste0("results/CCLE_",network_choice,"/sampled_topklist_consensus_combinations_rare_SL_results.csv"))
}

rare_SL_prediction_stats_summarised <- rare_SL_prediction_stats %>%
  group_by(algorithm,n) %>%
  summarise(
    mean_TP = mean(TP),
    mean_precision = mean(precision),
    mean_recall = mean(recall),
    mean_F1 = median(F1),
    sample_size = n()
  ) %>%
  pivot_longer(cols = -c(algorithm,n,sample_size), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure, levels = c("mean_TP","mean_recall","mean_precision","mean_F1")))


# Find best overall

suppressWarnings(rm(ranked_methods))
for(index in topn[1:2]){
  for(type in unique(rare_SL_prediction_stats_summarised$measure)){
    rank <- rare_SL_prediction_stats_summarised %>% filter(n==index, measure==type) %>% arrange(desc(value)) %>% pull(algorithm)
    rank <- list(rank)
    names(rank) <- paste0(type,"_",index)
    if(!exists("ranked_methods")){
      ranked_methods <- rank
    }else{
      ranked_methods <- append(ranked_methods,rank)
    }
  }
}



top_methods <- TopKLists::Borda(ranked_methods)[[1]] %>%
  mutate(rank=row_number()) %>%
  pivot_longer(cols = c(mean,median,geo.mean,l2norm), names_to = "measure", values_to = "consensus_algorithm") %>%
  filter(measure=="geo.mean") %>%
  arrange(rank)

top_method <- "SCS,DawnRank,sysSVM2_median"
#top_method <- top_methods %>% head(1) %>% pull(consensus_algorithm)

#alg_colours <- c("red",rep("black",length(unique(SL_prediction_stats_summarised$algorithm))-1))
#names(alg_colours) <- c(top_method, unique(SL_prediction_stats_summarised$algorithm)[which(unique(SL_prediction_stats_summarised$algorithm)!=top_method)] )

#alg_widths <- alg_colours
#alg_widths[alg_widths=="red"] <- 2
#alg_widths[alg_widths=="black"] <- 0.5
#alg_widths <- as.numeric(alg_widths)
#names(alg_widths) <- names(alg_colours)

plot_data <- rare_SL_prediction_stats_summarised %>%
  mutate(plot_colour = ifelse(algorithm==top_method, "red", "black"),
         plot_linewidth = ifelse(algorithm==top_method, 2, 0.5),
         plot_linetype = ifelse(algorithm==top_method, "dashed", "solid"))


ggplot(mapping=aes(x = n, y = value, color = algorithm, linewidth=algorithm, linetype=algorithm)) +
  geom_line(data=plot_data %>% filter(algorithm != top_method)) +
  geom_line(data=plot_data %>% filter(algorithm == top_method)) +
  scale_color_manual(breaks = plot_data$algorithm, values = plot_data$plot_colour) +
  scale_linewidth_manual(breaks = plot_data$algorithm, values = plot_data$plot_linewidth) +
  scale_linetype_manual(breaks = plot_data$algorithm, values = plot_data$plot_linetype) +
  ylab("Average Value") +
  xlab("Number of Predicted Sensitive Genes") +
  ggtitle(top_method) +
  theme_light() +
  theme(legend.position = "none") +
  facet_wrap(~measure, nrow = 1, scales = "free")

ggsave(paste0("results/CCLE_",network_choice,"/choosing_top_consensus_method_rare_SL_partner.png"),width = 40, height = 30, units = "cm", dpi = 300)


















#############################
# Evaluate Ref Drivers
#############################

# Get all altered genes


mutation <- fread(paste0("validation_data/CCLE_",network_choice,"/mutations.csv"), select = c("gene_ID",selected_samples)) %>% 
  arrange(gene_ID) %>%
  column_to_rownames("gene_ID")

cnv <- fread(paste0("validation_data/CCLE_",network_choice,"/cnv.csv"), select = c("gene_ID",selected_samples)) %>%
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

# Reference Drivers

CGC_all <- read_csv("data/cancer_reference_genes/CGC/cancer_gene_census.csv") %>%
  dplyr::select(gene_ID = `Gene Symbol`) %>%
  unique()


gs <- CGC_all %>% pull(gene_ID) %>% unique()

# A minimum number of gold standards allowed to each cell. If less are available, the cell is removed from calculations
min_gs <- 10


if(!file.exists(paste0("results/CCLE_",network_choice,"/sampled_topklist_consensus_combinations_ref_driver_results.csv"))){
  
  ref_prediction_stats <- foreach(sample=selected_samples, .combine = "rbind", .packages = c("tidyverse","foreach"), .export = c("consensus_drivers","altered_genes","min_gs")) %dopar% {
    
    alt_genes <- altered_genes %>% filter(cell_ID==sample) %>% pull(gene_ID)
    
    gs_specific <- gs[which(gs %in% alt_genes)]
    
    if(length(gs_specific)<min_gs){message("Skipping a cell due to too few gold standards")}else{
      
      
      
      tmp1 <- consensus_drivers %>% filter(cell_ID == sample) %>% mutate(algorithm=paste0(consensus_algorithms,"_",type)) %>% dplyr::select(algorithm, cell_ID, driver, rank)
      
      foreach(alg=unique(tmp1$algorithm), .combine = "rbind", .packages = c("tidyverse","foreach"), .export = c("run_evaluation","tmp1")) %do% {
        
        tmp2 <- tmp1 %>% filter(algorithm == alg)
        if(nrow(tmp2)==0){break}else{
          max_drivers <- max(tmp2 %>% pull(rank))
          
          print(paste0("Calculating stats for ", sample, "(",which(samples == sample),"/",length(samples),")", " targets using ", alg))
          
          foreach(n=topn, .combine = "rbind", .packages = c("tidyverse"), .export = "run_evaluation") %do% {
            
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
            
            data.frame(cell_ID = sample, 
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
    
  }
  write_csv(ref_prediction_stats, paste0("results/CCLE_",network_choice,"/sampled_topklist_consensus_combinations_ref_driver_results.csv"))
}else{
  message("stats file already exists. Reading in previous results.")
  message(paste0("To stop this, remove file: ", "results/CCLE_",network_choice,"/sampled_topklist_consensus_combinations_ref_driver_results.csv"))
  ref_prediction_stats <- read_csv(paste0("results/CCLE_",network_choice,"/sampled_topklist_consensus_combinations_ref_driver_results.csv"))
}


ref_prediction_stats_summarised <- ref_prediction_stats %>%
  group_by(algorithm,n) %>%
  summarise(
    mean_TP = mean(TP),
    mean_precision = mean(precision),
    mean_recall = mean(recall),
    mean_F1 = median(F1),
    sample_size = n()
  ) %>%
  pivot_longer(cols = -c(algorithm,n,sample_size), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure, levels = c("mean_TP","mean_recall","mean_precision","mean_F1")))


# Find best overall

suppressWarnings(rm(ranked_methods))
for(index in topn[1:2]){
  for(type in unique(ref_prediction_stats_summarised$measure)){
    rank <- ref_prediction_stats_summarised %>% filter(n==index, measure==type) %>% arrange(desc(value)) %>% pull(algorithm)
    rank <- list(rank)
    names(rank) <- paste0(type,"_",index)
    if(!exists("ranked_methods")){
      ranked_methods <- rank
    }else{
      ranked_methods <- append(ranked_methods,rank)
    }
  }
}



top_methods <- TopKLists::Borda(ranked_methods)[[1]] %>%
  mutate(rank=row_number()) %>%
  pivot_longer(cols = c(mean,median,geo.mean,l2norm), names_to = "measure", values_to = "consensus_algorithm") %>%
  filter(measure=="geo.mean") %>%
  arrange(rank)

#top_method <- "CSN_NCUA,sysSVM2,PhenoDriverR_mean"
top_method <- "CSN_NCUA,PersonaDrive,PRODIGY,OncoImpact,sysSVM2,PhenoDriverR_geo.mean"
#top_method <- top_methods %>% head(1) %>% pull(consensus_algorithm)

#alg_colours <- c("red",rep("black",length(unique(SL_prediction_stats_summarised$algorithm))-1))
#names(alg_colours) <- c(top_method, unique(SL_prediction_stats_summarised$algorithm)[which(unique(SL_prediction_stats_summarised$algorithm)!=top_method)] )

#alg_widths <- alg_colours
#alg_widths[alg_widths=="red"] <- 2
#alg_widths[alg_widths=="black"] <- 0.5
#alg_widths <- as.numeric(alg_widths)
#names(alg_widths) <- names(alg_colours)

plot_data <- ref_prediction_stats_summarised %>%
  mutate(plot_colour = ifelse(algorithm==top_method, "red", "black"),
         plot_linewidth = ifelse(algorithm==top_method, 2, 0.5),
         plot_linetype = ifelse(algorithm==top_method, "dashed", "solid"))


ggplot(mapping=aes(x = n, y = value, color = algorithm, linewidth=algorithm, linetype=algorithm)) +
  geom_line(data=plot_data %>% filter(algorithm != top_method)) +
  geom_line(data=plot_data %>% filter(algorithm == top_method)) +
  scale_color_manual(breaks = plot_data$algorithm, values = plot_data$plot_colour) +
  scale_linewidth_manual(breaks = plot_data$algorithm, values = plot_data$plot_linewidth) +
  scale_linetype_manual(breaks = plot_data$algorithm, values = plot_data$plot_linetype) +
  ylab("Average Value") +
  xlab("Number of Predicted Driver Genes") +
  ggtitle(top_method) +
  theme_light() +
  theme(legend.position = "none") +
  facet_wrap(~measure, nrow = 1, scales = "free")

ggsave(paste0("results/CCLE_",network_choice,"/choosing_top_consensus_method_reference_drivers.png"),width = 40, height = 30, units = "cm", dpi = 300)
