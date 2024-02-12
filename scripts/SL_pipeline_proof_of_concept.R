#SL_pipeline_proof_of_concept.R



# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(foreach, quietly = T))
suppressPackageStartupMessages (library(doParallel, quietly = T))
suppressPackageStartupMessages (library(ggpubr, quietly = T))
suppressPackageStartupMessages (library(cowplot, quietly = T))




# Handling input arguments
option_list = list(
  make_option(c("-s", "--synleth_partners"), type="integer", default=5, 
              help="Maximum number of SL-partners to use for each gene", metavar ="Synthetic Lethal Partners"),
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network being used", metavar ="Network"),
  make_option(c("-a", "--algorithms"), type="character", default="sysSVM2", 
              help="algorithms to include in comparison separated by spaces, or 'ALL' (Default), or lead with '-' to exclude", metavar ="Algorithms"),
  make_option(c("-t", "--threads"), type="integer", default=4, 
              help="Number of threads to use if running in parallel", metavar ="Threads"),
  make_option(c("-b", "--badDriver"), type="character", default="yes", 
              help="Include badDriver predictions (yes/no)", metavar ="badDriver")
  
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

max_SL <- opt$synleth_partners
network_choice <- opt$network
algorithms <- opt$algorithms
threads <- opt$threads
include_badDriver <- tolower(opt$badDriver)

SL_th <- 0.5

threads <- 8



if(threads>1){
  #registerDoParallel(cores=threads)
  cl <- makeCluster(threads, outfile = "log/proof_of_concept.log")
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
# Read In Results
#############################

#algorithms <- "-combined_de_novo_methods"
#algorithms <- "ALL"
#algorithms <- "sysSVM2"

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
      #warning("Invalid input for --algorithms. Cannot mix inclusions and exclusions in one statement.")
      #stop()
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
# badDriver Results
#############################

suppressWarnings(rm(badDriver_results))
if(include_badDriver == "yes" | include_badDriver == "y"){
  for(ct in list.dirs(paste0("bad_driver_simulations/CCLE_",network_choice,"/"), recursive = F)){
    ct <- gsub(".*/","",ct)
    message("badDriver results for ", ct)
    
    for(run in list.files(paste0("bad_driver_simulations/CCLE_",network_choice,"/",ct))) {
      alg <- gsub("_[0-9]+.csv","",run)
      run <- str_match(run,"_ref_([0-9]+).csv")[2]
      
      indiv_result <- read_csv(paste0("bad_driver_simulations/CCLE_",network_choice,"/",ct,"/", alg, "_", run, ".csv"), col_types = cols()) %>%
        mutate(algorithm=paste0(alg,"_",run))
      
      if(!exists("badDriver_results")){
        badDriver_results <- indiv_result
      }else{
        badDriver_results <- rbind(badDriver_results,indiv_result)
      }
      
    }
    
  }
  
  badDriver_results <- badDriver_results %>% filter(rank <= 100)
  
  aggregated_results <- aggregated_results %>% rbind(badDriver_results)
}




#############################
# Read In Gold Standards
#############################


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
#############################
# Read In Synthetic Lethality Predictions
#############################

SL <- read_csv(paste0("validation_data/CCLE_",network_choice,"/SL_partners_max_",max_SL,".csv")) %>%
  filter(score>=SL_th) %>%
  dplyr::select(gene_ID=gene1,partner=gene2,SL_rank=rank)


#############################
# Read In LOF/GOF Predictions
#############################

LOF_GOF <- read_csv("data/LOF_GOF_annotations.csv")





get_results <- function(ranked_predictions, sample_list, gold_standard_list){
  results <- foreach(sample=sample_list, .combine = "rbind", .packages = c("tidyverse","foreach"), .export = "sample_info") %dopar% {
    lineage <- sample_info %>% filter(cell_ID==sample) %>% pull("lineage")
    gs <- gold_standard_list[sample] %>% unlist()
    tmp1 <- ranked_predictions %>% filter(cell_ID == sample)
    
    foreach(alg=unique(tmp1$algorithm), .combine = "rbind", .packages = c("tidyverse","foreach")) %do% {
      
      tmp2 <- tmp1 %>% filter(algorithm == alg)
      if(nrow(tmp2)==0){break}
      max_drivers <- max(tmp2 %>% pull(final_rank))
      
      print(paste0("Calculating stats for ", sample, "(",which(sample_list == sample),"/",length(sample_list),")", " targets using ", alg))
      
      foreach(n=1:max_drivers, .combine = "rbind", .packages = c("tidyverse")) %do% {
        
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
        
        data.frame(lineage=lineage,
                   cell_ID = sample, 
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
  return(results)
}

get_rare_results <- function(results){
  
  all_algs <- unique(results$algorithm)  
  # Loop through each algorithm
  rare_results <- foreach(alg=all_algs, .combine = "rbind") %do% {
    
    # Get all predicted essential genes for the algorithm
    tmp1 <- results %>% filter(algorithm==alg)
    
    # Loop through each sample
    sample_result <- foreach(sample=unique(tmp1$cell_ID), .combine = "rbind",.packages = c("tidyverse","coin", "foreach")) %dopar% {
      
      print(paste0("Mann-Whitney U Test to find rare predicted essential genes. Algorithm: ", 
                   "(",which(all_algs == alg),"/",length(all_algs),"). Sample: "
                   , "(",which(unique(tmp1$cell_ID) == sample),"/",length(unique(tmp1$cell_ID)),")"))
      
      # Get the predicted essential genes for an individual sample
      sample_predictions <- tmp1 %>% filter(cell_ID==sample) 
      
      # Get the lineage of the sample
      lin <- sample_predictions %>% pull(lineage) %>% head(1)
      
      sample_predictions <- sample_predictions %>% pull(target)
      
      # Get the maximum number of predicted essential genes by that algorithm in that lineage
      max_pred <- tmp1 %>% filter(lineage==lin) %>% group_by(cell_ID) %>% summarise(count=n()) %>% pull(count) %>% max()
      
      # Get ranks of each predicted essential gene from all samples in the lineage
      lineage_ranks <- tmp1 %>% 
        filter(lineage==lin) %>%
        dplyr::select(cell_ID,target,final_rank) %>%
        pivot_wider(names_from = "target", values_from = "final_rank", values_fill = max_pred+1)
      
      
      gene_result <- foreach(gene=sample_predictions, .combine = "rbind") %do% {
        
        # Get the rank of a predicted essential gene in the sample compared with others in the lineage
        sample_rank <- lineage_ranks[c("cell_ID",gene)] %>%
          mutate(cell_ID=ifelse(cell_ID!=sample, "other", "sample")) %>%
          mutate(cell_ID=factor(cell_ID, levels = c("sample","other"))) %>%
          rename(rank=2) %>%
          arrange(cell_ID)
        
        
        stat <- wilcox_test(rank~cell_ID, data=sample_rank, conf.level=0.95,ties.method="mid-ranks", alternative = "less", paired=FALSE)
        
        data.frame(
          algorithm=alg,
          lineage=lin,
          cell_ID=sample,
          target=gene,
          final_rank=sample_rank[1,"rank"],
          pvalue=pvalue(stat),
          effect_size=as.numeric(stat@statistic@linearstatistic),
          expectation=as.numeric(expectation(stat))
        )
        
        
      }
      
    }
    
    
    sample_result %>% mutate(padj=p.adjust(sample_result$pvalue, method = "BH"))
    
    
  }
  
  rare_results <- rare_results %>%
    dplyr::rename(final_rank=rank) %>%
    filter(padj < 0.05) %>%
    group_by(algorithm,lineage,cell_ID) %>%
    arrange(final_rank) %>%
    mutate(final_rank=row_number()) %>%
    ungroup() %>%
    arrange(algorithm,lineage,cell_ID,final_rank)
  
  return(rare_results)
}


#############################
# vs Gold Standards
#############################

available_samples <- intersect(aggregated_results$cell_ID, names(gold_standard))


# Driver Level Only

driver_vs_gs <- get_results(aggregated_results %>% dplyr::rename(final_rank=rank, target=driver) %>% filter(final_rank <= 50), available_samples, gold_standard)


# Only Consider GOF Mutations

GOF_results <- aggregated_results %>%
  dplyr::rename(target=driver) %>%
  left_join(LOF_GOF, by = c("target"="gene_ID","cell_ID")) %>%
  filter(annotation == "GOF") %>%
  group_by(cell_ID, algorithm) %>%
  filter(!duplicated(target)) %>%
  filter(!is.na(target)) %>%
  mutate(final_rank = row_number()) %>%
  filter(final_rank <= 50)
  
gof_driver_vs_gs <- get_results(GOF_results, available_samples, gold_standard)


# All SL Partners

SL_aggregated_results <- aggregated_results %>%
  dplyr::rename(driver_rank=rank) %>%
  # Add SL partners
  left_join(SL, by = c("driver"="gene_ID"), multiple="all") %>%
  # Make SL partner the target
  dplyr::rename(target=partner) %>%
  arrange(driver_rank,SL_rank) %>%
  group_by(cell_ID, algorithm) %>%
  # For each cell, after arranging by driver rank, then SL_rank, remove duplicated target genes
  filter(!duplicated(target)) %>%
  # Finally remove missing targets (LOF mutations without an SL partner)
  filter(!is.na(target)) %>%
  mutate(final_rank = row_number()) %>%
  ungroup() %>%
  arrange(lineage,algorithm, cell_ID, final_rank) %>%
  filter(final_rank <= 50)
  
  
SL_vs_gs <- get_results(SL_aggregated_results, available_samples, gold_standard)

# Lof SL Partners Only

LOF_SL_aggregated_results <- aggregated_results %>%
  dplyr::rename(driver_rank=rank) %>%
  left_join(LOF_GOF, by = c("driver"="gene_ID","cell_ID")) %>%
  # assume that missing annotations are LOF
  mutate(annotation=ifelse(is.na(annotation),"LOF",annotation)) %>%
  filter(annotation=="LOF") %>%
  # Add SL partners
  left_join(SL, by = c("driver"="gene_ID"), multiple="all") %>%
  mutate(target =  partner) %>%
  arrange(driver_rank,SL_rank) %>%
  group_by(cell_ID, algorithm) %>%
  # For each cell, after arranging by driver rank, then SL_rank, remove duplicated target genes
  filter(!duplicated(target)) %>%
  # Finally remove missing targets (LOF mutations without an SL partner)
  filter(!is.na(target)) %>%
  mutate(final_rank = row_number()) %>%
  ungroup() %>%
  arrange(lineage,algorithm, cell_ID, final_rank) %>%
  filter(final_rank <= 50)

LOF_SL_vs_gs <- get_results(LOF_SL_aggregated_results, available_samples, gold_standard)


# Full Pipeline

full_pipeline_aggregated_results <- aggregated_results %>%
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
full_pipeline_aggregated_results <- full_pipeline_aggregated_results %>%
  rbind(full_pipeline_aggregated_results %>% 
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
  arrange(lineage,algorithm, cell_ID, final_rank) #%>%
  #filter(final_rank <= 50)

full_pipeline_aggregated_results_rare <- get_rare_results(full_pipeline_aggregated_results) %>%
  filter(final_rank < 50)

full_pipeline_aggregated_results <- full_pipeline_aggregated_results %>%
  filter(final_rank <= 50)


full_pipeline_vs_gs <- get_results(full_pipeline_aggregated_results, available_samples, gold_standard)



all_stats <- rbind(
  driver_vs_gs %>% mutate(pipeline="Drivers Only"),
  gof_driver_vs_gs %>% mutate(pipeline = "GOF Drivers Only"),
  SL_vs_gs %>% mutate(pipeline="All SL Partners"),
  LOF_SL_vs_gs %>% mutate(pipeline="LOF SL Partners Only"),
  full_pipeline_vs_gs %>% mutate(pipeline="Full Pipeline")
)

write_csv(all_stats, paste0("results/CCLE_",network_choice,"/proof_of_concept_results_all_gold_standards.csv"))


pipeline_colours <- c("Drivers Only"="black","GOF Drivers Only"="red","All SL Partners"="blue","LOF SL Partners Only"="purple","Full Pipeline"="green")

all_stats_summarised <- all_stats %>%
  # Combining all badDriver simulations to one mean
  mutate(algorithm=ifelse(str_detect(algorithm,"bad_sim"), "badDriver", algorithm)) %>%
  # Only keeping results for samples with > 10 gold standards
  filter(length_gs >= 10) %>%
  group_by(pipeline,algorithm,n) %>%
  # Only keep measurements where more than 10 cells are available to calculate mean
  filter(n()>=10) %>%
  summarise(
    mean_TP = mean(TP),
    mean_precision = mean(precision),
    mean_recall = mean(recall),
    mean_F1 = median(F1),
    sample_size = n()
  ) %>%
  pivot_longer(cols = -c(pipeline,algorithm,n,sample_size), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure, levels = c("mean_TP","mean_recall","mean_precision","mean_F1")),
         algorithm = factor(algorithm, levels = c(algorithms,"badDriver")),
         pipeline = factor(pipeline, levels = names(pipeline_colours)),
         sample_size = ifelse(algorithm=="badDriver",length(available_samples),sample_size))

max_predictions<-50

# Dynamic N_max
gs_lengths <- all_stats %>% dplyr::select(cell_ID,length_gs) %>% unique() %>% pull(length_gs)


if(max_predictions=="median gold standards"){
  ## Calculates 2*median number of sensitive genes for cell lines with >3 sensitive genes
  N_max <- gs_lengths[which(gs_lengths  > 3 )] %>% median()*2
} else {
  N_max <- max_predictions
}




ggplot(all_stats_summarised %>% filter(n<=N_max), mapping=aes(x = n, y = value, color = pipeline, linetype = algorithm)) +
  geom_line() +
  scale_color_manual(breaks = names(pipeline_colours),values = pipeline_colours) +
  scale_linetype_manual(breaks = c(algorithms,"badDriver"), values = c(rep("solid",length(algorithms)), "dashed")) +
  ylab("Average Value") +
  xlab("Number of Predicted Sensitive Genes") +
  guides(linetype=guide_legend(title="Algorithm (Linetype)"),colour=guide_legend(title="Pipeline (Colour)")) +
  theme_light() +
  facet_wrap(~measure, nrow = 1, scales = "free")

ggsave(paste0("results/CCLE_",network_choice,"/pipeline_proof_of_concept_all_gold_standards.png"),width = 50, height = 20, units = "cm", dpi = 300)

#############################
# vs Rare Gold Standards
#############################

driver_vs_rare_gs <- get_results(aggregated_results %>% dplyr::rename(final_rank=rank, target=driver) %>% filter(final_rank <= 50), available_samples, rare_gold_standard)
gof_driver_vs_rare_gs <- get_results(GOF_results, available_samples, rare_gold_standard)
SL_vs_rare_gs <- get_results(SL_aggregated_results, available_samples, rare_gold_standard)
LOF_SL_vs_rare_gs <- get_results(LOF_SL_aggregated_results, available_samples, rare_gold_standard)
full_pipeline_vs_rare_gs <- get_results(full_pipeline_aggregated_results, available_samples, rare_gold_standard)
full_pipeline_rare_vs_rare_gs <- get_results(full_pipeline_aggregated_results_rare, available_samples, rare_gold_standard)


all_stats_rare <- rbind(
  driver_vs_rare_gs %>% mutate(pipeline="Drivers Only"),
  gof_driver_vs_rare_gs %>% mutate(pipeline = "GOF Drivers Only"),
  SL_vs_rare_gs %>% mutate(pipeline="All SL Partners"),
  LOF_SL_vs_rare_gs %>% mutate(pipeline="LOF SL Partners Only"),
  full_pipeline_vs_rare_gs %>% mutate(pipeline="Full Pipeline All"),
  full_pipeline_rare_vs_rare_gs %>% mutate(pipeline="Full Pipeline Rare")
)

write_csv(all_stats_rare, paste0("results/CCLE_",network_choice,"/proof_of_concept_results_rare_gold_standards.csv"))


pipeline_colours <- c("Drivers Only"="black","GOF Drivers Only"="red","All SL Partners"="blue","LOF SL Partners Only"="purple","Full Pipeline All"="green", "Full Pipeline Rare"="orange")

all_stats_rare_summarised <- all_stats_rare %>%
  # Combining all badDriver simulations to one mean
  mutate(algorithm=ifelse(str_detect(algorithm,"bad_sim"), "badDriver", algorithm)) %>%
  # Only keeping results for samples with > 10 gold standards
  filter(length_gs >= 10) %>%
  group_by(pipeline,algorithm,n) %>%
  # Only keep measurements where more than 10 cells are available to calculate mean
  filter(n()>=10) %>%
  summarise(
    mean_TP = mean(TP),
    mean_precision = mean(precision),
    mean_recall = mean(recall),
    mean_F1 = median(F1),
    sample_size = n()
  ) %>%
  pivot_longer(cols = -c(pipeline,algorithm,n,sample_size), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure, levels = c("mean_TP","mean_recall","mean_precision","mean_F1")),
         algorithm = factor(algorithm, levels = c(algorithms,"badDriver")),
         pipeline = factor(pipeline, levels = names(pipeline_colours)),
         sample_size = ifelse(algorithm=="badDriver",length(available_samples),sample_size))

max_predictions<-50

# Dynamic N_max
gs_lengths <- all_stats_rare %>% dplyr::select(cell_ID,length_gs) %>% unique() %>% pull(length_gs)


if(max_predictions=="median gold standards"){
  ## Calculates 2*median number of sensitive genes for cell lines with >3 sensitive genes
  N_max <- gs_lengths[which(gs_lengths  > 3 )] %>% median()*2
} else {
  N_max <- max_predictions
}




ggplot(all_stats_rare_summarised %>% filter(n<=N_max), mapping=aes(x = n, y = value, color = pipeline, linetype = algorithm)) +
  geom_line() +
  scale_color_manual(breaks = names(pipeline_colours),values = pipeline_colours) +
  scale_linetype_manual(breaks = c(algorithms,"badDriver"), values = c(rep("solid",length(algorithms)), "dashed")) +
  ylab("Average Value") +
  xlab("Number of Predicted Sensitive Genes") +
  guides(linetype=guide_legend(title="Algorithm (Linetype)"),colour=guide_legend(title="Pipeline (Colour)")) +
  theme_light() +
  facet_wrap(~measure, nrow = 1, scales = "free")

ggsave(paste0("results/CCLE_",network_choice,"/proof_of_concept_results_rare_gold_standards.png"),width = 50, height = 20, units = "cm", dpi = 300)




#########################################
# Rare Predictions vs Rare Gold Standards
#########################################

available_samples <- intersect(aggregated_results$cell_ID, names(rare_gold_standard))





rare_drivers <- get_rare_results(aggregated_results %>% dplyr::rename(final_rank=rank, target=driver))


#global_freq_thresh <- 0.5
#lineage_freq_thresh <- 0.25

#freq_n_max <- 50

#rare_driver_vs_rare_gs_results <- aggregated_results %>% 
#  dplyr::rename(final_rank=rank, target=driver) %>% 
  #dplyr::select(lineage, cell_ID, algorithm, target, final_rank) %>%
#  ungroup() %>%
  # Because some algorithms give very long lists, they will be greatly affected by trying to find only rare drivers
  # need to filter to only high ranks first
#  filter(final_rank <= freq_n_max) %>%
  # Get the total # of cells
#  mutate(total_n = length(unique(cell_ID))) %>%
  # Get the global frequency of predicted essential gene for each algorithm
#  group_by(algorithm,target) %>% mutate(global_essentiality_count = length(unique(cell_ID))) %>% ungroup() %>%
  # Get the total # of cells per lineage
#  group_by(lineage) %>% mutate(lineage_n = length(unique(cell_ID))) %>% ungroup() %>%
  # Get the lineage frequency of essentiality
#  group_by(algorithm,lineage,target) %>% mutate(lineage_essentiality_count = length(unique(cell_ID))) %>% ungroup() %>%
#  mutate(global_essentiality_frequency = global_essentiality_count / total_n,
#         lineage_essentiality_frequency = lineage_essentiality_count / lineage_n) %>%
  # Filter frequency thresholds
#  filter(global_essentiality_frequency < global_freq_thresh, lineage_essentiality_frequency < lineage_freq_thresh) %>%
#  dplyr::select(lineage, cell_ID, algorithm, target, final_rank, global_essentiality_frequency, lineage_essentiality_frequency) %>%
  # update the ranks
#  group_by(algorithm,cell_ID) %>%
#  arrange(final_rank) %>%
#  mutate(final_rank = row_number()) %>%
#  ungroup() %>%
#  arrange(algorithm,lineage,cell_ID,final_rank) %>%
#  filter(final_rank <= 50)

rare_driver_vs_rare_gs <- get_results(rare_drivers, available_samples, rare_gold_standard)



rare_gof_drivers <- get_rare_results(aggregated_results  %>%
                                       dplyr::rename(target=driver) %>%
                                       left_join(LOF_GOF, by = c("target"="gene_ID","cell_ID")) %>%
                                       filter(annotation == "GOF") %>%
                                       group_by(cell_ID, algorithm) %>%
                                       filter(!duplicated(target)) %>%
                                       filter(!is.na(target)) %>%
                                       mutate(final_rank = row_number()) %>%
                                       ungroup()
                                     )


#rare_gof_driver_vs_rare_gs_results <- aggregated_results %>%
#  dplyr::rename(target=driver) %>%
#  left_join(LOF_GOF, by = c("target"="gene_ID","cell_ID")) %>%
#  filter(annotation == "GOF") %>%
#  group_by(cell_ID, algorithm) %>%
#  filter(!duplicated(target)) %>%
#  filter(!is.na(target)) %>%
#  mutate(final_rank = row_number()) %>%
#  ungroup() %>%
  # Because some algorithms give very long lists, they will be greatly affected by trying to find only rare drivers
  # need to filter to only high ranks first
#  filter(final_rank <= freq_n_max) %>%
  # Get the total # of cells
#  mutate(total_n = length(unique(cell_ID))) %>%
  # Get the global frequency of predicted essential gene for each algorithm
#  group_by(algorithm,target) %>% mutate(global_essentiality_count = length(unique(cell_ID))) %>% ungroup() %>%
  # Get the total # of cells per lineage
#  group_by(lineage) %>% mutate(lineage_n = length(unique(cell_ID))) %>% ungroup() %>%
  # Get the lineage frequency of essentiality
#  group_by(algorithm,lineage,target) %>% mutate(lineage_essentiality_count = length(unique(cell_ID))) %>% ungroup() %>%
#  mutate(global_essentiality_frequency = global_essentiality_count / total_n,
#         lineage_essentiality_frequency = lineage_essentiality_count / lineage_n) %>%
  # Filter frequency thresholds
#  filter(global_essentiality_frequency < global_freq_thresh, lineage_essentiality_frequency < lineage_freq_thresh) %>%
#  dplyr::select(lineage, cell_ID, algorithm, target, final_rank, global_essentiality_frequency, lineage_essentiality_frequency) %>%
  # update the ranks
#  group_by(algorithm,cell_ID) %>%
#  arrange(final_rank) %>%
#  mutate(final_rank = row_number()) %>%
#  ungroup() %>%
#  arrange(algorithm,lineage,cell_ID,final_rank) %>%
#  filter(final_rank <= 50)


rare_gof_driver_vs_rare_gs <- get_results(rare_gof_drivers %>% dplyr::rename(final_rank=rank), available_samples, rare_gold_standard)




#rare_SL_vs_rare_gs_results <- aggregated_results %>%
#  dplyr::rename(driver_rank=rank) %>%
  # Add SL partners
#  left_join(SL, by = c("driver"="gene_ID"), multiple="all") %>%
  # Make SL partner the target
#  dplyr::rename(target=partner) %>%
#  arrange(driver_rank,SL_rank) %>%
#  group_by(cell_ID, algorithm) %>%
  # For each cell, after arranging by driver rank, then SL_rank, remove duplicated target genes
#  filter(!duplicated(target)) %>%
  # Finally remove missing targets (LOF mutations without an SL partner)
#  filter(!is.na(target)) %>%
#  mutate(final_rank = row_number()) %>%
#  ungroup() %>%
#  arrange(lineage,algorithm, cell_ID, final_rank) %>%
#  left_join(LOF_GOF, by = c("target"="gene_ID","cell_ID")) %>%
#  filter(annotation == "GOF") %>%
#  group_by(cell_ID, algorithm) %>%
#  filter(!duplicated(target)) %>%
#  filter(!is.na(target)) %>%
#  mutate(final_rank = row_number()) %>%
#  ungroup() %>%
  # Because some algorithms give very long lists, they will be greatly affected by trying to find only rare drivers
  # need to filter to only high ranks first
#  filter(final_rank <= freq_n_max) %>%
  # Get the total # of cells
#  mutate(total_n = length(unique(cell_ID))) %>%
  # Get the global frequency of predicted essential gene for each algorithm
#  group_by(algorithm,target) %>% mutate(global_essentiality_count = length(unique(cell_ID))) %>% ungroup() %>%
  # Get the total # of cells per lineage
#  group_by(lineage) %>% mutate(lineage_n = length(unique(cell_ID))) %>% ungroup() %>%
  # Get the lineage frequency of essentiality
#  group_by(algorithm,lineage,target) %>% mutate(lineage_essentiality_count = length(unique(cell_ID))) %>% ungroup() %>%
#  mutate(global_essentiality_frequency = global_essentiality_count / total_n,
#         lineage_essentiality_frequency = lineage_essentiality_count / lineage_n) %>%
  # Filter frequency thresholds
#  filter(global_essentiality_frequency < global_freq_thresh, lineage_essentiality_frequency < lineage_freq_thresh) %>%
#  dplyr::select(lineage, cell_ID, algorithm, target, final_rank, global_essentiality_frequency, lineage_essentiality_frequency) %>%
  # update the ranks
#  group_by(algorithm,cell_ID) %>%
#  arrange(final_rank) %>%
#  mutate(final_rank = row_number()) %>%
#  ungroup() %>%
#  arrange(algorithm,lineage,cell_ID,final_rank) %>%
#  filter(final_rank <= 50)


rare_SL_vs_rare_gs <- get_results(rare_SL_vs_rare_gs_results, available_samples, rare_gold_standard)



rare_LOF_SL_vs_rare_gs_results <- aggregated_results %>%
  dplyr::rename(driver_rank=rank) %>%
  left_join(LOF_GOF, by = c("driver"="gene_ID","cell_ID")) %>%
  # assume that missing annotations are LOF
  mutate(annotation=ifelse(is.na(annotation),"LOF",annotation)) %>%
  filter(annotation=="LOF") %>%
  # Add SL partners
  left_join(SL, by = c("driver"="gene_ID"), multiple="all") %>%
  mutate(target =  partner) %>%
  arrange(driver_rank,SL_rank) %>%
  group_by(cell_ID, algorithm) %>%
  # For each cell, after arranging by driver rank, then SL_rank, remove duplicated target genes
  filter(!duplicated(target)) %>%
  # Finally remove missing targets (LOF mutations without an SL partner)
  filter(!is.na(target)) %>%
  mutate(final_rank = row_number()) %>%
  ungroup() %>%
  arrange(lineage,algorithm, cell_ID, final_rank) %>%
  # Because some algorithms give very long lists, they will be greatly affected by trying to find only rare drivers
  # need to filter to only high ranks first
  filter(final_rank <= freq_n_max) %>%
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
  arrange(algorithm,lineage,cell_ID,final_rank) %>%
  filter(final_rank <= 50)
  

rare_LOF_SL_vs_rare_gs <- get_results(rare_LOF_SL_vs_rare_gs_results, available_samples, rare_gold_standard)




rare_full_pipeline_vs_rare_gs_results <- aggregated_results %>%
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
rare_full_pipeline_vs_rare_gs_results <- rare_full_pipeline_vs_rare_gs_results %>%
  rbind(rare_full_pipeline_vs_rare_gs_results %>% 
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
  arrange(lineage,algorithm, cell_ID, final_rank) %>%
  # Because some algorithms give very long lists, they will be greatly affected by trying to find only rare drivers
  # need to filter to only high ranks first
  filter(final_rank <= freq_n_max) %>%
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
  arrange(algorithm,lineage,cell_ID,final_rank) %>%
  filter(final_rank <= 50)


rare_full_pipeline_vs_rare_gs <- get_results(rare_full_pipeline_vs_rare_gs_results, available_samples, rare_gold_standard)



all_stats_rare_vs_rare <- rbind(
  rare_driver_vs_rare_gs %>% mutate(pipeline="Drivers Only"),
  rare_gof_driver_vs_rare_gs %>% mutate(pipeline = "GOF Drivers Only"),
  rare_SL_vs_rare_gs %>% mutate(pipeline="All SL Partners"),
  rare_LOF_SL_vs_rare_gs %>% mutate(pipeline="LOF SL Partners Only"),
  rare_full_pipeline_vs_rare_gs %>% mutate(pipeline="Full Pipeline")
)

write_csv(all_stats_rare_vs_rare, paste0("results/CCLE_",network_choice,"/proof_of_concept_rare_results_rare_gold_standards.csv"))


pipeline_colours <- c("Drivers Only"="black","GOF Drivers Only"="red","All SL Partners"="blue","LOF SL Partners Only"="purple","Full Pipeline"="green")

all_stats_rare_vs_rare_summarised <- all_stats_rare_vs_rare %>%
  # Combining all badDriver simulations to one mean
  mutate(algorithm=ifelse(str_detect(algorithm,"bad_sim"), "badDriver", algorithm)) %>%
  # Only keeping results for samples with > 10 gold standards
  filter(length_gs >= 10) %>%
  group_by(pipeline,algorithm,n) %>%
  # Only keep measurements where more than 10 cells are available to calculate mean
  filter(n()>=10) %>%
  summarise(
    mean_TP = mean(TP),
    mean_precision = mean(precision),
    mean_recall = mean(recall),
    mean_F1 = median(F1),
    sample_size = n()
  ) %>%
  pivot_longer(cols = -c(pipeline,algorithm,n,sample_size), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure, levels = c("mean_TP","mean_recall","mean_precision","mean_F1")),
         algorithm = factor(algorithm, levels = c(algorithms,"badDriver")),
         pipeline = factor(pipeline, levels = names(pipeline_colours)),
         sample_size = ifelse(algorithm=="badDriver",length(available_samples),sample_size))

max_predictions<-50

# Dynamic N_max
gs_lengths <- all_stats_rare_vs_rare %>% dplyr::select(cell_ID,length_gs) %>% unique() %>% pull(length_gs)


if(max_predictions=="median gold standards"){
  ## Calculates 2*median number of sensitive genes for cell lines with >3 sensitive genes
  N_max <- gs_lengths[which(gs_lengths  > 3 )] %>% median()*2
} else {
  N_max <- max_predictions
}




ggplot(all_stats_rare_vs_rare_summarised %>% filter(n<=N_max), mapping=aes(x = n, y = value, color = pipeline, linetype = algorithm)) +
  geom_line() +
  scale_color_manual(breaks = names(pipeline_colours),values = pipeline_colours) +
  scale_linetype_manual(breaks = c(algorithms,"badDriver"), values = c(rep("solid",length(algorithms)), "dashed")) +
  ylab("Average Value") +
  xlab("Number of Predicted Sensitive Genes") +
  guides(linetype=guide_legend(title="Algorithm (Linetype)"),colour=guide_legend(title="Pipeline (Colour)")) +
  theme_light() +
  facet_wrap(~measure, nrow = 1, scales = "free")

ggsave(paste0("results/CCLE_",network_choice,"/proof_of_concept_rare_results_rare_gold_standards.png"),width = 50, height = 20, units = "cm", dpi = 300)


















############################
# Gene Effect Density Plots
############################


# Altered Genes


## Mutations

CGC <- read_csv("data/cancer_reference_genes/CGC/cancer_gene_census.csv") %>%
  dplyr::select(gene_ID = `Gene Symbol`) %>%
  pull(gene_ID) %>%
  unique()



mutation <- fread(paste0("validation_data/CCLE_",network_choice,"/mutations.csv"), select = c("gene_ID",samples)) %>% 
  filter(gene_ID %in% genes) %>%
  arrange(gene_ID) %>%
  pivot_longer(-gene_ID, names_to = "cell_ID", values_to = "mutated") %>%
  filter(mutated!=0) %>%
  dplyr::select(-mutated)

cnv <- fread(paste0("validation_data/CCLE_",network_choice,"/cnv.csv"), select = c("gene_ID",samples)) %>%
  filter(gene_ID %in% genes) %>%
  arrange(gene_ID) %>%
  pivot_longer(-gene_ID, names_to = "cell_ID", values_to = "mutated") %>%
  filter(mutated!=0) %>%
  dplyr::select(-mutated)

altered_genes <- rbind(mutation,cnv) %>% 
  unique() %>%
  filter(gene_ID %in% CGC)


LOF_GOF <- fread("data/LOF_GOF_annotations_all_evidence.csv") %>%
  #filter(LOF_total_score >=5 | GOF_total_score >= 5) %>%
  dplyr::select(cell_ID,gene_ID,annotation=combined_prediction) %>%
  unique() %>%
  ungroup() %>%
  # Need to indicate if a single gene has both GOF and LOF mutations
  group_by(cell_ID,gene_ID) %>%
  summarise(annotation=paste0(annotation, collapse = "and")) %>%
  ungroup() %>%
  mutate(annotation=ifelse(str_detect(annotation,"and"),"both",annotation))


picklesv3 <- fread("data/PICKLESv3/master_table_Avana_Score_TKOv3_BF_Zscore_Expr_LOF_GOF_18959genes_1162cells_lineage_disease.txt") %>%
  dplyr::select(-c(V1,LOF,GOF,primary_disease,lineage,Expr,Cell_Line)) %>%
  dplyr::rename(cell_ID=stripped_cell_line_name, gene_ID=Gene) %>%
  relocate(cell_ID) %>%
  unique() %>%
  filter(cell_ID %in% samples, gene_ID %in% genes)

#ggplot(LOF_GOF) +
#  geom_density(aes(x=LOF_total_score), colour="red") +
#  geom_density(aes(x=GOF_total_score), colour="blue")



gene_effect <- fread("data/CCLE/CRISPRGeneEffect.csv") %>%
  as.data.frame() %>%
  dplyr::rename(DepMap_ID=V1) %>%
  inner_join(sample_info[c("DepMap_ID","cell_ID","lineage")], by = "DepMap_ID") %>%
  dplyr::select(-DepMap_ID) %>%
  pivot_longer(-c(cell_ID,lineage), names_to = "gene_ID", values_to = "gene_effect") %>%
  mutate(gene_ID = gsub(" \\(.*","", gene_ID)) %>%
  filter(gene_ID %in% genes, cell_ID %in% samples) %>%
  group_by(cell_ID,gene_ID) %>%
  summarise(gene_effect=min(gene_effect)) %>%
  ungroup()



density_all_mutations <- inner_join(altered_genes,gene_effect, by = c("cell_ID", "gene_ID")) %>% 
  mutate(group="All Mutations")

density_GOF <- altered_genes %>%
  left_join(LOF_GOF, by = c("gene_ID","cell_ID")) %>%
  filter(annotation == "GOF") %>%
  dplyr::select(-annotation) %>% 
  inner_join(gene_effect, by = c("cell_ID", "gene_ID")) %>% 
  mutate(group="GOF Mutations")

density_LOF_SL <- altered_genes %>%
  left_join(LOF_GOF, by = c("gene_ID","cell_ID")) %>%
  # assume that missing annotations are LOF
  mutate(annotation=ifelse(is.na(annotation),"LOF",annotation)) %>%
  filter(annotation=="LOF") %>%
  # Add SL partners
  left_join(SL, by = c("gene_ID"), multiple="all") %>%
  dplyr::select(cell_ID, gene_ID = partner) %>%
  unique() %>% 
  inner_join(gene_effect, by = c("cell_ID", "gene_ID")) %>% 
  mutate(group="LOF SL Partners")

density_full_pipeline <- altered_genes %>%
  left_join(LOF_GOF, by = c("gene_ID","cell_ID")) %>%
  # assume that missing annotations are LOF
  mutate(annotation=ifelse(is.na(annotation),"LOF",annotation)) %>%
  # If mutation is GOF, add the gene itself to the target list
  mutate(target = ifelse(annotation %in% c("GOF"), gene_ID, NA)) %>%
  # Add SL partners
  left_join(SL, by = c("gene_ID"), multiple="all") %>%
  # Remove SL partners for GOF mutations
  mutate(partner = ifelse(annotation == "GOF", NA, partner),
         SL_rank = ifelse(annotation == "GOF", NA, SL_rank)) %>%
  unique() %>%
  # Make the SL partner the target for LOF and both
  mutate(target = ifelse(annotation %in% c("LOF", "both"), partner, target))

# Add the original gene_ID as an additional target for genes with both GOF and LOF mutations
density_full_pipeline <- density_full_pipeline %>%
  rbind(density_full_pipeline %>% 
          filter(annotation == "both") %>% 
          mutate(target = gene_ID, SL_rank = 0, partner = NA) %>%
          unique()) %>%
  dplyr::select(cell_ID, gene_ID = target) %>%
  unique() %>% 
  inner_join(gene_effect, by = c("cell_ID", "gene_ID")) %>% 
  mutate(group="Full Pipeline")
  

pipeline_density_plot_data <- rbind(
  density_all_mutations,
  density_GOF,
  density_LOF_SL,
  density_full_pipeline
) %>%
  left_join(picklesv3, by = c("cell_ID", "gene_ID"))

pipeline_medians <- pipeline_density_plot_data %>%
  group_by(group) %>%
  summarise(across(is.numeric, ~ median(.x, na.rm=T)))


ggplot(pipeline_density_plot_data, aes(x=Avana_BF, colour = group)) +
  geom_density(linewidth=1) +
  geom_vline(data = pipeline_medians, aes(xintercept=Avana_BF, colour=group), linewidth=1, linetype="dashed") +
  xlim(-100, 100)


ggsave("plots/pipeline_essentiality.png")









LOF_GOF_Density <- altered_genes %>%
  inner_join(LOF_GOF, by = c("gene_ID","cell_ID")) %>%
  #inner_join(gene_effect, by = c("cell_ID", "gene_ID"))
  inner_join(picklesv3, by = c("cell_ID", "gene_ID"))

LOF_GOF_medians <- LOF_GOF_Density %>%
  group_by(annotation) %>%
  summarise(across(is.numeric, ~ median(.x, na.rm=T)))

ggplot(LOF_GOF_Density, aes(x=Avana_BF, colour=annotation)) +
  geom_density(linewidth=1) +
  geom_vline(data = LOF_GOF_medians, aes(xintercept=Avana_BF, colour=annotation), linewidth=1, linetype="dashed") +
  xlim(-100, 100)

ggsave("plots/LOF_GOF_essentiality.png")





















###############################################
# Using Tier 1 CGC Drivers for Proof of Concept
###############################################

#1 Gather cell lines with mutations in tier 1 CGC drivers

## CGC

CGC <- read_csv("data/cancer_reference_genes/CGC/cancer_gene_census.csv") %>%
  filter(Tier==1) %>%
  dplyr::select(gene_ID = `Gene Symbol`) %>%
  pull(gene_ID) %>%
  unique()


## Altered Genes

mutation <- fread(paste0("validation_data/CCLE_",network_choice,"/mutations.csv"), select = c("gene_ID",samples)) %>% 
  #filter(gene_ID %in% genes) %>%
  arrange(gene_ID) %>%
  pivot_longer(-gene_ID, names_to = "cell_ID", values_to = "mutated") %>%
  filter(mutated!=0) %>%
  dplyr::select(-mutated)

cnv <- fread(paste0("validation_data/CCLE_",network_choice,"/cnv.csv"), select = c("gene_ID",samples)) %>%
  #filter(gene_ID %in% genes) %>%
  arrange(gene_ID) %>%
  pivot_longer(-gene_ID, names_to = "cell_ID", values_to = "mutated") %>%
  filter(mutated!=0) %>%
  dplyr::select(-mutated)

altered_genes <- rbind(mutation,cnv) %>% 
  unique() %>%
  filter(gene_ID %in% CGC)


SL_th <- 0.5

#for(SL_th in c(0,0.1,0.25,0.5,0.75,0.9)){

# LOF_GOF
LOF_GOF <- read_csv("data/LOF_GOF_annotations.csv")
# SL
SL <- read_csv(paste0("validation_data/CCLE_",network_choice,"/SL_partners_max_",max_SL,".csv")) %>%
  filter(score >= SL_th) %>%
  dplyr::select(gene_ID=gene1,partner=gene2)


#2 Convert these to predicted essential genes

predicted_essential_genes <- altered_genes %>%
  left_join(LOF_GOF, by = c("cell_ID","gene_ID")) %>%
  dplyr::rename(driver=gene_ID) %>%
  # assume that missing annotations are LOF
  #mutate(annotation=ifelse(is.na(annotation),"LOF",annotation)) %>%
  # If mutation is GOF, add the gene itself to the target list
  mutate(target = ifelse(annotation %in% c("GOF"), driver, NA)) %>%
  # Add SL partners
  left_join(SL, by = c("driver"="gene_ID"), multiple="all") %>%
  # Remove SL partners for GOF mutations
  mutate(partner = ifelse(annotation == "GOF", NA, partner),
         ) %>%
  unique() %>%
  # Make the SL partner the target for LOF and both
  mutate(target = ifelse(annotation %in% c("LOF", "both"), partner, target))

# Add the original gene_ID as an additional target for genes with both GOF and LOF mutations
predicted_essential_genes <- predicted_essential_genes %>%
  rbind(predicted_essential_genes %>% 
          filter(annotation == "both") %>% 
          mutate(annotation = "GOF") %>%
          mutate(target = driver, partner = NA) %>%
          unique()) %>%
  mutate(annotation=ifelse(annotation=="both","LOF",annotation)) %>%
  filter(!is.na(target)) %>%
  mutate(predicted_essential_gene="All Predicted Essential Genes") %>%
  dplyr::select(cell_ID,gene_ID=target,annotation,predicted_essential_gene) %>%
  pivot_longer(cols = c(annotation,predicted_essential_gene), values_to = "annotation", names_to = "type") %>%
  dplyr::select(-type) %>%
  unique()
  
  

# Compare Essentiality

essentiality <- fread("data/PICKLESv3/master_table_Avana_Score_TKOv3_BF_Zscore_Expr_LOF_GOF_18959genes_1162cells_lineage_disease.txt") %>%
  dplyr::select(-c(V1,LOF,GOF,primary_disease,lineage,Expr,Cell_Line,TKOv3_Z,TKOv3_BF)) %>%
  dplyr::rename(cell_ID=stripped_cell_line_name, gene_ID=Gene) %>%
  relocate(cell_ID) %>%
  left_join(predicted_essential_genes, by = c("cell_ID","gene_ID")) %>%
  mutate(annotation = ifelse(is.na(annotation), "No Annotation", annotation))

# Keep all rows for predicted essential genes
keep_rows <- which(essentiality$annotation!="No Annotation")
# Also keep a random sample of 10000 other rows
keep_rows <- append(keep_rows, sample(which(!1:nrow(essentiality) %in% keep_rows), 50000))

#keep_rows <- 1:nrow(essentiality)

essentiality_plot <- essentiality[keep_rows,] %>%
  pivot_longer(Avana_Z:Score_BF, names_to = "measure", values_to = "value") %>%
  filter(measure=="Avana_BF") %>%
  mutate(annotation = ifelse(annotation=="LOF","LOF Drivers - SL Partner", annotation)) %>%
  mutate(annotation = ifelse(annotation=="GOF","GOF Drivers", annotation)) %>%
  mutate(annotation = factor(annotation, levels = 
    c("No Annotation","All Predicted Essential Genes","GOF Drivers","LOF Drivers - SL Partner")
  ))

test <- essentiality_plot %>% 
  dplyr::select(annotation,value) %>%
  arrange(annotation) %>%
  na.omit()

write_csv(test, "test.csv")


my_comparisons <- list( c("No Annotation", "All Predicted Essential Genes"), 
                        c("No Annotation", "GOF Drivers"),
                        c("No Annotation", "LOF Drivers - SL Partner")
                        )

pairwise_stats <- compare_means(value~annotation,data=essentiality_plot,p.adjust.method = "fdr",method = "wilcox.test", comparisons=my_comparisons) %>%
  filter(group1=="No Annotation"|group2=="No Annotation") %>%
  mutate(custom=ifelse(p.signif=="****", "<0.0001", ""))

all_stats <- compare_means(value~annotation,data=essentiality_plot,p.adjust.method = "fdr",method = "kruskal.test") %>%
  mutate(custom=ifelse(p.signif=="****", "<0.0001", ""))


ggplot(essentiality_plot, aes(x=annotation, y=value)) +
  geom_boxplot() +
  #facet_wrap(~measure, nrow = 1, scales = "free") +
  stat_pvalue_manual(pairwise_stats,label="p = {custom}", y.position = c(240,255,270)) +
  #stat_pvalue_manual(all_stats,label="Kruskal-Wallis p = {custom}") +
  ylab("Avana Bagel2 LogBF Score (Essentiality)") +
  xlab("") +
  theme_bw()

ggsave("plots/paper_fig_CGC_tier_1_pipeline_proof_of_concept.png", units = "px", dpi = 300)
#ggsave(paste0("plots/pipeline_poc_CGCtier1_sl_threshold_",SL_th,".png"))




#3 Convert these to predicted drug sensitivities

drug_targets <- fread(paste0("validation_data/CCLE_",network_choice,"/drug_targets.csv")) %>%
  dplyr::select(drug_ID,gene_ID,source) %>%
  unique()

drug_sensitivity <- fread("validation_data/drug_sensitivity.csv")

# foreach loop through sources

drugs_plot <- foreach(src=c("GDSC1","GDSC2","PRISM"), .combine = "rbind") %do% {
  ## get cell-specific predicted drug sensitivty for each source (predicted_essential_genes join with drug targets, copy code from eval_drug_sens)
  predicted_drugs <- predicted_essential_genes %>%
    dplyr::select(cell_ID,gene_ID,annotation) %>%
    left_join(drug_targets, by = "gene_ID", relationship="many-to-many") %>%
    filter(!is.na(drug_ID)) %>%
    dplyr::select(-gene_ID) %>%
    #mutate(predicted_sensitive_drug=T) %>%
    unique()
  
  drug_sensitivity <- fread("validation_data/drug_sensitivity.csv") %>%
    filter(sensitivity_source==src) %>%
    dplyr::select(cell_ID,drug_ID,sensitivity_value,source=sensitivity_source) %>%
    left_join(predicted_drugs, by = c("cell_ID","drug_ID","source")) %>%
    mutate(annotation = ifelse(is.na(annotation), "No Annotation", annotation)) %>%
    mutate(annotation = ifelse(annotation=="LOF","LOF Drivers - SL Partner", annotation)) %>%
    mutate(annotation = ifelse(annotation=="GOF","GOF Drivers", annotation)) %>%
    mutate(annotation = factor(annotation, levels = 
                                 c("No Annotation","All Predicted Essential Genes","GOF Drivers","LOF Drivers - SL Partner")
    ))
  
  drug_sensitivity
  
  ## New "keep_rows" for all predicted senstive drugs + x others
  
  # Keep all rows for predicted essential genes
  #keep_rows <- which(drug_sensitivity$predicted_sensitive_drug)
  # Also keep a random sample of 10000 other rows
  #keep_rows <- append(keep_rows, sample(which(!1:nrow(drug_sensitivity) %in% keep_rows), 10000))
  
  #drug_sensitivity[keep_rows,]
  
}



my_comparisons <- list( c("No Annotation", "All Predicted Essential Genes"), 
                        c("No Annotation", "GOF Drivers"),
                        c("No Annotation", "LOF Drivers - SL Partner")
)

# GDSC1

pairwise_stats <- compare_means(sensitivity_value~annotation,data=drugs_plot %>% filter(source=="GDSC1"),p.adjust.method = "fdr",method = "wilcox.test", comparisons=my_comparisons) %>%
  filter(group1=="No Annotation"|group2=="No Annotation") %>%
  mutate(custom=ifelse(p.signif=="****", "<0.0001", ""))


GDSC1_p <- ggplot(drugs_plot %>% filter(source=="GDSC1"), aes(x=annotation, y=sensitivity_value)) +
  geom_boxplot() +
  ylab("GDSC1 lnIC50") +
  xlab("") +
  stat_pvalue_manual(pairwise_stats,label="p = {custom}", y.position = c(8,9,10)) +
  theme_bw()

#ggsave(paste0("plots/pipeline_poc_CGCtier1_sl_threshold_",SL_th,"_drugs.png"))


# GDSC2

pairwise_stats <- compare_means(sensitivity_value~annotation,data=drugs_plot %>% filter(source=="GDSC2"),p.adjust.method = "fdr",method = "wilcox.test", comparisons=my_comparisons) %>%
  filter(group1=="No Annotation"|group2=="No Annotation") %>%
  mutate(custom=ifelse(p.signif=="****", "<0.0001", ""))


GDSC2_p <- ggplot(drugs_plot %>% filter(source=="GDSC2"), aes(x=annotation, y=sensitivity_value)) +
  geom_boxplot() +
  ylab("GDSC2 lnIC50") +
  xlab("") +
  stat_pvalue_manual(pairwise_stats,label="p = {custom}", y.position = c(11,12,13)) +
  theme_bw()


# PRISM

pairwise_stats <- compare_means(sensitivity_value~annotation,data=drugs_plot %>% filter(source=="PRISM"),p.adjust.method = "fdr",method = "wilcox.test", comparisons=my_comparisons) %>%
  filter(group1=="No Annotation"|group2=="No Annotation") %>%
  mutate(custom=ifelse(p.signif=="****", "<0.0001", ""))


PRISM_p <- ggplot(drugs_plot %>% filter(source=="PRISM"), aes(x=annotation, y=sensitivity_value)) +
  geom_boxplot() +
  ylab("PRISM LFC Viability") +
  xlab("") +
  stat_pvalue_manual(pairwise_stats,label="p = {custom}", y.position = c(4,5,6)) +
  theme_bw()





cowplot::plot_grid(GDSC1_p,GDSC2_p,PRISM_p)


ggsave("plots/paper_fig_CGC_tier_1_pipeline_proof_of_concept_drugs.png", units = "px", dpi = 300, width = 4000, height = 4000)







####################################
# Gold Standards LOF/GOF Proportions
####################################


# Altered Genes


## Mutations

mutation <- fread(paste0("validation_data/CCLE_",network_choice,"/mutations.csv"), select = c("gene_ID",samples)) %>% 
  filter(gene_ID %in% genes) %>%
  arrange(gene_ID) %>%
  pivot_longer(-gene_ID, names_to = "cell_ID", values_to = "mutated") %>%
  filter(mutated!=0) %>%
  dplyr::select(-mutated) %>%
  unique()

cnv <- fread(paste0("validation_data/CCLE_",network_choice,"/cnv.csv"), select = c("gene_ID",samples)) %>%
  filter(gene_ID %in% genes) %>%
  arrange(gene_ID) %>%
  pivot_longer(-gene_ID, names_to = "cell_ID", values_to = "mutated") %>%
  filter(mutated!=0) %>%
  dplyr::select(-mutated) %>%
  unique()

altered_genes <- rbind(mutation,cnv) %>% 
  unique()
  

LOF_GOF <- fread("data/LOF_GOF_annotations_all_evidence.csv") %>%
  #filter(LOF_total_score >=5 | GOF_total_score >= 5) %>%
  dplyr::select(cell_ID,gene_ID,annotation=combined_prediction) %>%
  unique() %>%
  ungroup() %>%
  # Need to indicate if a single gene has both GOF and LOF mutations
  group_by(cell_ID,gene_ID) %>%
  summarise(annotation=paste0(annotation, collapse = "and")) %>%
  ungroup() %>%
  mutate(annotation=ifelse(str_detect(annotation,"and"),"both",annotation))

SL <- read_csv(paste0("validation_data/CCLE_",network_choice,"/SL_partners_max_",max_SL,".csv")) %>%
  filter(score>=SL_th) %>%
  dplyr::select(gene_ID=gene1,partner=gene2,SL_rank=rank)

# Gold Standards


gold_standard_df <- read_csv(paste0("validation_data/CCLE_",network_choice,"/gold_standards.csv")) %>%
  dplyr::select(cell_ID,gene_ID) %>%
  unique()

rare_gold_standard_df <- read_csv(paste0("validation_data/CCLE_",network_choice,"/rare_gold_standards.csv")) %>%
  dplyr::select(cell_ID,gene_ID) %>%
  unique()




plot_data <- rbind(
  altered_genes %>% mutate(group="Genetically Altered"),
  LOF_GOF %>% dplyr::rename(group=annotation),
  LOF_GOF %>% filter(annotation %in% c("LOF","both")) %>% left_join(SL, by="gene_ID") %>% dplyr::select(cell_ID,gene_ID=partner) %>% mutate(group="LOF SL Partner")
) %>%
  full_join(gold_standard_df %>% mutate(gold_standard=T), by=c("gene_ID","cell_ID")) %>%
  full_join(rare_gold_standard_df %>% mutate(rare_gold_standard=T), by=c("gene_ID","cell_ID")) %>%
  filter(gold_standard|rare_gold_standard) %>%
  mutate(group=ifelse(is.na(group), "No Annotation", group)) %>%
  mutate(group=ifelse(group=="both","LOF + GOF",group))

ggplot(plot_data %>% filter(gold_standard
                            ,group != "No Annotation"
                            ) %>% group_by(group) %>% summarise(count=n()), aes(x=group,y=count)) +
  geom_col()


