# evaluate_SL_partners.r

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(foreach, quietly = T))
suppressPackageStartupMessages (library(doParallel, quietly = T))
suppressPackageStartupMessages (library(ggpubr, quietly = T))
suppressPackageStartupMessages (library(coin, quietly = T))
suppressPackageStartupMessages (library(ggpointdensity, quietly = T))



# Handling input arguments
option_list = list(
  make_option(c("-s", "--synleth_partners"), type="integer", default=5, 
              help="Maximum number of SL-partners to use for each gene", metavar ="Synthetic Lethal Partners"),
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network being used", metavar ="Network"),
  make_option(c("-a", "--algorithms"), type="character", default="ALL", 
              help="algorithms to include in comparison separated by spaces, or 'ALL' (Default), or lead with '-' to exclude", metavar ="Algorithms"),
  make_option(c("-t", "--threads"), type="integer", default=4, 
              help="Number of threads to use if running in parallel", metavar ="Threads"),
  make_option(c("-b", "--badDriver"), type="character", default="yes", 
              help="Include badDriver predictions (yes/no)", metavar ="badDriver"),
  make_option(c("-c", "--celltype"), type="character", default="ALL", 
              help="Cell types to include in the analysis separated by semicolons, or 'ALL' (Default)", metavar ="badDriver")
  
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

max_SL <- opt$synleth_partners
network_choice <- opt$network
algorithms <- opt$algorithms
threads <- opt$threads
include_badDriver <- tolower(opt$badDriver)
celltypes <- opt$celltype


SL_th <- 0.5

#celltypes <- "Liver;Biliary_Tract"

if(celltypes != "ALL"){
  celltypes <- str_split_1(celltypes, pattern = ";")
}


threads <- 4

#network_choice <- "own"
#include_badDriver <- "no"

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


genes <- fread(paste0("validation_data/CCLE_",network_choice,"/tpm.csv"), select = "gene_ID") %>% pull(gene_ID) %>% unique()




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




if(celltypes[1] != "ALL"){
  aggregated_results <- aggregated_results %>% filter(lineage %in% celltypes)
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
##########################################
# Read In Synthetic Lethality Predictions
##########################################

SL <- read_csv(paste0("validation_data/CCLE_",network_choice,"/SL_partners_max_",max_SL,".csv")) %>%
  filter(score>=SL_th) %>%
  dplyr::select(gene_ID=gene1,partner=gene2,SL_rank=rank)


#############################
# Read In LOF/GOF Predictions
#############################

LOF_GOF <- read_csv("data/LOF_GOF_annotations.csv")

#############################
# N Max
#############################


# Calculates 2*median number of sensitive genes for cell lines with >3 sensitive genes
#N_max <- gold_standard[lapply(gold_standard, length) > 3] %>% lapply(length) %>% unlist() %>% median()*2



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
          mutate(target = driver, SL_rank = 0, partner = NA) %>%
          unique()) %>%
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


write_csv(SL_partner_results, paste0("results/CCLE_",network_choice,"/aggregated_predicted_essential_genes.csv"))



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
write_csv(SL_prediction_stats, paste0("results/CCLE_",network_choice,"/full_results_SL_partner_level.csv"))
}else{
  message("SL Prediction stats file already exists. Reading in previous results.")
  message(paste0("To stop this, remove file: ", "results/CCLE_",network_choice,"/full_results_SL_partner_level.csv"))
  SL_prediction_stats <- fread(paste0("results/CCLE_",network_choice,"/full_results_SL_partner_level.csv"))
}


if(celltypes[1] != "ALL"){
  SL_prediction_stats <- SL_prediction_stats %>% filter(lineage %in% celltypes)
}



#old <- read_csv(paste0("results/CCLE_",network_choice,"/full_results_SL_partner_level.csv"))
#new <- rbind(SL_prediction_stats, old) %>% unique()
#write_csv(new, paste0("results/CCLE_",network_choice,"/full_results_SL_partner_level.csv"))
#SL_prediction_stats <- new

alg_of_interest <- c(
  #"consensus_topklists_geomean_CSN_NCUA_PersonaDrive_PRODIGY_OncoImpact_sysSVM2_PhenoDriverR", # best for STRINGv11 reference drivers
  #"consensus_topklists_median_SCS_DawnRank_sysSVM2",
  #"consensus_topklists_geomean_CSN_NCUA_OncoImpact_sysSVM2",
  
  #"consensus_topklists_geomean_CSN_NCUA_OncoImpact_sysSVM2", # Best for STRINGv11 SL partners
  #"consensus_topklists_median_SCS_DawnRank_sysSVM2", # Best for STRINGv11 rare SL partners
  "consensus_topklists_median_CSN_NCUA_CSN_MDS_OncoImpact_sysSVM2",
  "consensus_topklists_geomean_DawnRank_OncoImpact_PersonaDrive",
  "consensus_topklists_geomean_ALL",
  "post_consensus_rare",
  
  #"consensus_topklists_geomean_sysSVM2_PNC", # Best for own SL partners
  #"consensus_topklists_geomean_PRODIGY_sysSVM2_PhenoDriverR", # Best for own rare SL partners
  
  
  "DawnRank",
  "OncoImpact",
  "PNC",
  "PRODIGY",
  "PersonaDrive",
  "SCS",
  "sysSVM2",
  "PhenoDriverR",
  
  #"CSN_DFVS",
  "CSN_MDS",
  #"CSN_MMS",
  "CSN_NCUA",
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
  #"SSN_NCUA",
  
  "badDriver"
)


alg_colours <- read_csv("data/algorithm_colours.csv") %>% deframe()


SL_prediction_stats_summarised <- SL_prediction_stats %>%
  # Combining all badDriver simulations to one mean
  mutate(algorithm=ifelse(str_detect(algorithm,"bad_sim"), "badDriver", algorithm)) %>%
  filter(algorithm %in% alg_of_interest) %>%
  # Only keeping results for samples with > 10 gold standards
  filter(length_gs >= 10) %>%
  group_by(algorithm,n) %>%
  # Only keep measurements where more than 10 cells are available to calculate mean
  filter(n()>=10) %>%
  summarise(
    mean_TP = mean(TP),
    mean_precision = mean(precision),
    mean_recall = mean(recall),
    mean_F1 = median(F1),
    sample_size = n()
  ) %>%
  pivot_longer(cols = -c(algorithm,n,sample_size), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure, levels = c("mean_TP","mean_recall","mean_precision","mean_F1")),
         algorithm = factor(algorithm, levels = names(alg_colours)),
         sample_size = ifelse(algorithm=="badDriver",length(tmp_samples),sample_size))

max_predictions<-50

# Dynamic N_max
gs_lengths <- SL_prediction_stats %>% dplyr::select(cell_ID,length_gs) %>% unique() %>% pull(length_gs)


if(max_predictions=="median gold standards"){
  ## Calculates 2*median number of sensitive genes for cell lines with >3 sensitive genes
  N_max <- gs_lengths[which(gs_lengths  > 3 )] %>% median()*2
} else {
  N_max <- max_predictions
}




ggplot(mapping=aes(x = n, y = value, color = algorithm)) +
  geom_line(data=SL_prediction_stats_summarised %>% filter(algorithm %in% alg_of_interest) %>% filter(n<=N_max) %>% filter(!str_detect(algorithm, "consensus|badDriver")), mapping=aes(alpha = sample_size)) +
  geom_line(data=SL_prediction_stats_summarised %>% filter(algorithm %in% alg_of_interest) %>% filter(n<=N_max) %>% filter(str_detect(algorithm, "consensus|badDriver")), linetype = "dashed") +
  scale_color_manual(values = alg_colours) +
  #scale_linetype_manual(
  #  breaks = names(alg_colours),
  #  values = ifelse(str_detect(names(alg_colours),"consensus|badDriver"),"dashed","solid")
  #) +
  ylab("Average Value") +
  xlab("Number of Predicted Sensitive Genes") +
  guides(colour=guide_legend(title="Algorithm (Colour)"),alpha=guide_legend(title = "Sample Size Remaining (Opacity)")) +
  theme_light() +
  facet_wrap(~measure, nrow = 1, scales = "free")

ggsave(paste0("results/CCLE_",network_choice,"/SL_partner_level_stats_max_",max_predictions,"_",paste(celltypes, collapse = "-"),".png"),width = 50, height = 20, units = "cm", dpi = 300)

# Stats Plot

stats_top_n <- 5
plot_measure="precision"

stats_plot <- SL_prediction_stats %>%
  ungroup() %>%
  filter(n==stats_top_n) %>%
  # Combining all badDriver simulations to one mean
  mutate(algorithm=ifelse(str_detect(algorithm,"bad_sim"), "badDriver", algorithm)) %>%
  filter(algorithm %in% alg_of_interest) %>%
  # Only keeping results for samples with > 10 gold standards
  filter(length_gs >= 10) %>%
  pivot_longer(cols = c(TP,FP,precision,recall,F1), names_to = "measure", values_to = "value") %>%
  mutate(algorithm = factor(algorithm, levels = names(alg_colours)))

# Precision

my_comparisons <- data.frame(a=rep("badDriver",length(alg_of_interest)-1),b=alg_of_interest[which(alg_of_interest!="badDriver")]) %>%
  transpose() %>%
  as.list()

ggplot(stats_plot %>% filter(measure==plot_measure), aes(x=algorithm, y=value, fill=algorithm)) +
  geom_violin() +
  scale_fill_manual(breaks = names(alg_colours), values = alg_colours) +
  guides(fill="none") +
  ggtitle(paste0(plot_measure," of Top ",stats_top_n," Predicted Essential Genes")) +
  labs(y=plot_measure, x="Algorithm" ) +
  stat_compare_means(method = "wilcox.test",ref.group = "badDriver", label = "p.signif", hide.ns = F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave(paste0("results/CCLE_",network_choice,"/SL_partner_level_stats_F1_top_",stats_top_n,".png"),width = 50, height = 20, units = "cm", dpi = 300)



# All Predictions vs Rare Essentials


tmp_samples <- intersect(SL_partner_results$cell_ID, names(rare_gold_standard))

if(!file.exists(paste0("results/CCLE_",network_choice,"/full_results_SL_partner_level_all_vs_rare.csv"))){
  
  SL_prediction_stats_all_vs_rare <- foreach(sample=tmp_samples, .combine = "rbind", .packages = c("tidyverse","foreach")) %dopar% {
    
    lineage <- sample_info %>% filter(cell_ID==sample) %>% pull("lineage")
    gs <- rare_gold_standard[sample] %>% unlist()
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
  write_csv(SL_prediction_stats_all_vs_rare, paste0("results/CCLE_",network_choice,"/full_results_SL_partner_level_all_vs_rare.csv"))
}else{
  message("SL Prediction stats file already exists. Reading in previous results.")
  message(paste0("To stop this, remove file: ", "results/CCLE_",network_choice,"/full_results_SL_partner_level_all_vs_rare.csv"))
  SL_prediction_stats_all_vs_rare <- fread(paste0("results/CCLE_",network_choice,"/full_results_SL_partner_level_all_vs_rare.csv"))
}


if(celltypes[1] != "ALL"){
  SL_prediction_stats_all_vs_rare <- SL_prediction_stats_all_vs_rare %>% filter(lineage %in% celltypes)
}



SL_prediction_all_vs_rare_stats_summarised <- SL_prediction_stats_all_vs_rare %>%
  # Combining all badDriver simulations to one mean
  mutate(algorithm=ifelse(str_detect(algorithm,"bad_sim"), "badDriver", algorithm)) %>%
  filter(algorithm %in% alg_of_interest) %>%
  # Only keeping results for samples with > 10 gold standards
  filter(length_gs >= 10) %>%
  group_by(algorithm,n) %>%
  # Only keep measurements where more than 10 cells are available to calculate mean
  filter(n()>=10) %>%
  summarise(
    mean_TP = mean(TP),
    mean_precision = mean(precision),
    mean_recall = mean(recall),
    mean_F1 = median(F1),
    sample_size = n()
  ) %>%
  pivot_longer(cols = -c(algorithm,n,sample_size), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure, levels = c("mean_TP","mean_recall","mean_precision","mean_F1")),
         algorithm = factor(algorithm, levels = names(alg_colours)),
         sample_size = ifelse(algorithm=="badDriver",length(tmp_samples),sample_size))

max_predictions<-10

# Dynamic N_max
gs_lengths <- SL_prediction_stats_all_vs_rare %>% dplyr::select(cell_ID,length_gs) %>% unique() %>% pull(length_gs)


if(max_predictions=="median gold standards"){
  ## Calculates 2*median number of sensitive genes for cell lines with >3 sensitive genes
  N_max <- gs_lengths[which(gs_lengths  > 3 )] %>% median()*2
} else {
  N_max <- max_predictions
}




ggplot(mapping=aes(x = n, y = value, color = algorithm)) +
  geom_line(data=SL_prediction_all_vs_rare_stats_summarised %>% filter(algorithm %in% alg_of_interest) %>% filter(n<=N_max) %>% filter(!str_detect(algorithm, "consensus|badDriver")), mapping=aes(alpha = sample_size)) +
  geom_line(data=SL_prediction_all_vs_rare_stats_summarised %>% filter(algorithm %in% alg_of_interest) %>% filter(n<=N_max) %>% filter(str_detect(algorithm, "consensus|badDriver")), linetype = "dashed") +
  scale_color_manual(values = alg_colours) +
  #scale_linetype_manual(
  #  breaks = names(alg_colours),
  #  values = ifelse(str_detect(names(alg_colours),"consensus|badDriver"),"dashed","solid")
  #) +
  ylab("Average Value") +
  xlab("Number of Predicted Sensitive Genes") +
  guides(colour=guide_legend(title="Algorithm (Colour)"),alpha=guide_legend(title = "Sample Size Remaining (Opacity)")) +
  theme_light() +
  facet_wrap(~measure, nrow = 1, scales = "free")

ggsave(paste0("results/CCLE_",network_choice,"/SL_partner_level_stats_all_vs_rare_max_",max_predictions,"_",paste(celltypes, collapse = "-"),".png"),width = 50, height = 20, units = "cm", dpi = 300)







################################
# Predicting Rare Sensitive Genes
################################

#main_algs <- c("DawnRank","OncoImpact","PersonaDrive","PhenoDriverR","PNC","PRODIGY","SCS","sysSVM2","CSN_NCUA")

global_freq_thresh <- 0.5
lineage_freq_thresh <- 0.25
#freq_n_max <- SL_partner_results %>% 
#  filter(algorithm %in% main_algs) %>% 
#  group_by(cell_ID,algorithm) %>% 
#  summarise(count = n()) %>% 
#  group_by(algorithm) %>% 
#  summarise(mean_count=mean(count)) %>% 
#  pull(mean_count) %>%
#  median()
freq_n_max <- 50

rare_SL_partner_predictions <- SL_partner_results %>%
  dplyr::select(lineage, cell_ID, algorithm, target, final_rank) %>%
  ungroup() %>%
  # Because some algorithms give very long lists, they will be greatly affected by trying to find only rare drivers
  # need to filter to only high ranks first
  filter(final_rank <= freq_n_max) %>%
  # Get the total # of cells analysed for each algorithm
  group_by(algorithm) %>% mutate(total_n = length(unique(cell_ID))) %>% ungroup() %>%
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

rare_prediction_counts <- rare_SL_partner_predictions %>%
  group_by(algorithm, cell_ID) %>%
  summarise(count = n()) %>%
  group_by(algorithm) %>%
  summarise(mean_count = mean(count))


tmp_samples <- intersect(rare_SL_partner_predictions$cell_ID, names(rare_gold_standard))

if(!file.exists(paste0("results/CCLE_",network_choice,"/full_results_SL_partner_level_rare.csv"))){
  
  rare_SL_prediction_stats <- foreach(sample=tmp_samples, .combine = "rbind", .packages = c("tidyverse","foreach")) %dopar% {
    
    lineage <- sample_info %>% filter(cell_ID==sample) %>% pull("lineage")
    gs <- rare_gold_standard[sample] %>% unlist()
    tmp1 <- rare_SL_partner_predictions %>% filter(cell_ID == sample)
    
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
  write_csv(rare_SL_prediction_stats, paste0("results/CCLE_",network_choice,"/full_results_SL_partner_level_rare.csv"))
}else{
  message("SL Prediction stats file already exists. Reading in previous results.")
  message(paste0("To stop this, remove file: ", "results/CCLE_",network_choice,"/full_results_SL_partner_level_rare.csv"))
  rare_SL_prediction_stats <- read_csv(paste0("results/CCLE_",network_choice,"/full_results_SL_partner_level_rare.csv"))
}


alg_colours <- read_csv("data/algorithm_colours.csv") %>% deframe()


rare_SL_prediction_stats_summarised <- rare_SL_prediction_stats %>%
  # Combining all badDriver simulations to one mean
  mutate(algorithm=ifelse(str_detect(algorithm,"bad_sim"), "badDriver", algorithm)) %>%
  filter(algorithm %in% alg_of_interest) %>%
  # Only keeping results for samples with > 10 gold standards
  filter(length_gs >= 10) %>%
  group_by(algorithm,n) %>%
  # Only keep measurements where more than 10 cells are available to calculate mean
  filter(n()>=10) %>%
  summarise(
    mean_TP = mean(TP),
    #mean_PPV = mean(PPV),
    mean_precision = mean(precision),
    mean_recall = mean(recall),
    mean_F1 = median(F1),
    #median_TP = median(TP),
    #median_PPV = median(PPV),
    #median_precision = median(precision),
    #median_recall = median(recall),
    #median_F1 = median(F1)
    sample_size = n()
  ) %>%
  pivot_longer(cols = -c(algorithm,n,sample_size), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure, levels = c("mean_TP","mean_recall","mean_precision","mean_F1")),
         algorithm = factor(algorithm, levels = names(alg_colours)),
         sample_size = ifelse(algorithm=="badDriver",length(samples),sample_size))

max_predictions<-50

# Dynamic N_max
gs_lengths <- rare_SL_prediction_stats %>% dplyr::select(cell_ID,length_gs) %>% unique() %>% pull(length_gs)


if(max_predictions=="median gold standards"){
  ## Calculates 2*median number of sensitive genes for cell lines with >3 sensitive genes
  N_max <- gs_lengths[which(gs_lengths  > 3 )] %>% median()*2
} else {
  N_max <- max_predictions
}




ggplot(mapping=aes(x = n, y = value, color = algorithm)) +
  geom_line(data=rare_SL_prediction_stats_summarised %>% filter(algorithm %in% alg_of_interest) %>% filter(n<=N_max) %>% filter(!str_detect(algorithm, "consensus|badDriver")), mapping=aes(alpha = sample_size)) +
  geom_line(data=rare_SL_prediction_stats_summarised %>% filter(algorithm %in% alg_of_interest) %>% filter(n<=N_max) %>% filter(str_detect(algorithm, "consensus|badDriver")), linetype = "dashed") +
  scale_color_manual(values = alg_colours) +
  #scale_linetype_manual(
  #  breaks = names(alg_colours),
  #  values = ifelse(str_detect(names(alg_colours),"consensus|badDriver"),"dashed","solid")
  #) +
  ylab("Average Value") +
  xlab("Number of Predicted Sensitive Genes") +
  guides(colour=guide_legend(title="Algorithm (Colour)"),alpha=guide_legend(title = "Sample Size Remaining (Opacity)")) +
  theme_light() +
  facet_wrap(~measure, nrow = 1, scales = "free")

ggsave(paste0("results/CCLE_",network_choice,"/rare_SL_partner_level_stats_max_",max_predictions,".png"),width = 50, height = 20, units = "cm", dpi = 300)



alg_colours <- read_csv("data/algorithm_colours.csv") %>% deframe()

lineage_plot_measure <- "mean_precision"


lineage_rare_SL_prediction_stats_summarised <- rare_SL_prediction_stats %>%
  # Combining all badDriver simulations to one mean
  mutate(algorithm=ifelse(str_detect(algorithm,"bad_sim"), "badDriver", algorithm)) %>%
  filter(algorithm %in% alg_of_interest) %>%
  # Only keeping results for samples with > 10 gold standards
  filter(length_gs >= 10) %>%
  group_by(lineage,algorithm,n) %>%
  # Only keep measurements where more than 10 cells are available to calculate mean
  filter(n()>=10) %>%
  summarise(
    mean_TP = mean(TP),
    #mean_PPV = mean(PPV),
    mean_precision = mean(precision),
    mean_recall = mean(recall),
    mean_F1 = median(F1),
    #median_TP = median(TP),
    #median_PPV = median(PPV),
    #median_precision = median(precision),
    #median_recall = median(recall),
    #median_F1 = median(F1)
    sample_size = n()
  ) %>%
  pivot_longer(cols = -c(lineage,algorithm,n,sample_size), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure, levels = c("mean_TP","mean_recall","mean_precision","mean_F1")),
         algorithm = factor(algorithm, levels = names(alg_colours)),
         sample_size = ifelse(algorithm=="badDriver",length(samples),sample_size)) %>%
  filter(measure==lineage_plot_measure)

N_max<-50


ggplot(mapping=aes(x = n, y = value, color = algorithm)) +
  geom_line(data=lineage_rare_SL_prediction_stats_summarised %>% filter(algorithm %in% alg_of_interest) %>% filter(n<=N_max) %>% filter(!str_detect(algorithm, "consensus|badDriver")), mapping=aes(alpha = sample_size)) +
  geom_line(data=lineage_rare_SL_prediction_stats_summarised %>% filter(algorithm %in% alg_of_interest) %>% filter(n<=N_max) %>% filter(str_detect(algorithm, "consensus|badDriver")), linetype = "dashed") +
  scale_color_manual(values = alg_colours) +
  #scale_linetype_manual(
  #  breaks = names(alg_colours),
  #  values = ifelse(str_detect(names(alg_colours),"consensus|badDriver"),"dashed","solid")
  #) +
  ylab("Average Value") +
  xlab("Number of Predicted Sensitive Genes") +
  guides(colour=guide_legend(title="Algorithm (Colour)"),alpha=guide_legend(title = "Sample Size Remaining (Opacity)")) +
  theme_light() +
  facet_wrap(~lineage, scales = "fixed")

ggsave(paste0("results/CCLE_",network_choice,"/rare_SL_partner_level_stats_max_",N_max,"_by_lineage.png"),width = 60, height = 40, units = "cm", dpi = 300)



























# Using Mann-Whitney U Test (New Idea 16.11.23)

if(!file.exists(paste0("results/CCLE_",network_choice,"/rare_essential_genes_mann_whitney.csv"))){

all_algs <- unique(SL_partner_results$algorithm)
  
# Loop through each algorithm
rare_SL_partner_predictions <- foreach(alg=all_algs, .combine = "rbind", .packages = c("tidyverse","coin", "foreach")) %dopar% {
  
  # Get all predicted essential genes for the algorithm
  tmp1 <- SL_partner_results %>% filter(algorithm==alg)
  
  # Loop through each sample
  sample_result <- foreach(sample=unique(tmp1$cell_ID), .combine = "rbind") %do% {
    
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
write_csv(rare_SL_partner_predictions, paste0("results/CCLE_",network_choice,"/rare_essential_genes_mann_whitney.csv"))
}else{
  message("Rare essential genes file already exists. Reading in previous results.")
  message(paste0("To stop this, remove file: ", "results/CCLE_",network_choice,"/rare_essential_genes_mann_whitney.csv"))
  rare_SL_partner_predictions <- read_csv(paste0("results/CCLE_",network_choice,"/rare_essential_genes_mann_whitney.csv"))
}

# Idea 1

# get p-value-based ranks
# borda count aggregate driver-based ranking and rareness ranking

#combined_ranking_rare_predictions <- rare_SL_partner_predictions %>%
  # Get p-value based rankings
#  group_by(algorithm,lineage,cell_ID) %>%
#  arrange(padj) %>%
#  mutate(p_rank=row_number()) %>%
  # Perform Borda count
#  mutate(driver_rank_borda_score=(length(target)-rank)^2,
#         p_rank_borda_score=(length(target)-p_rank)
#         ) %>%
#  mutate(driver_rank_borda_score=(driver_rank_borda_score - min(driver_rank_borda_score)) / (max(driver_rank_borda_score) - min(driver_rank_borda_score)),
#         p_rank_borda_score=(p_rank_borda_score - min(p_rank_borda_score)) / (max(p_rank_borda_score) - min(p_rank_borda_score))
#  ) %>%
  #rowwise() %>%
#  mutate(borda_score=driver_rank_borda_score + p_rank_borda_score) %>%
#  group_by(algorithm,lineage,cell_ID) %>%
#  arrange(desc(borda_score)) %>%
#  mutate(borda_rank=row_number()) %>%
#  ungroup() %>%
#  arrange(algorithm,lineage,cell_ID,borda_rank)


# Idea 2: Filter for significant increase
combined_ranking_rare_predictions <- rare_SL_partner_predictions %>%
  filter(padj < 0.05) %>%
  group_by(algorithm,lineage,cell_ID) %>%
  arrange(rank) %>%
  mutate(rank=row_number()) %>%
  ungroup() %>%
  arrange(algorithm,lineage,cell_ID,rank)



post_consensus_rare <- get_consensus(df=combined_ranking_rare_predictions, algs=c("DawnRank","OncoImpact","PersonaDrive"), title="post_consensus_rare")

combined_ranking_rare_predictions <- combined_ranking_rare_predictions %>%
  rbind(post_consensus_rare %>% mutate(pvalue=NA,effect_size=NA,expectation=NA,padj=NA))


tmp_samples <- intersect(combined_ranking_rare_predictions$cell_ID, names(rare_gold_standard))

if(!file.exists(paste0("results/CCLE_",network_choice,"/full_results_mann_whitney_rankings.csv"))){
  
  rare_SL_prediction_stats <- foreach(sample=tmp_samples, .combine = "rbind", .packages = c("tidyverse","foreach")) %dopar% {
    
    lineage <- sample_info %>% filter(cell_ID==sample) %>% pull("lineage")
    gs <- rare_gold_standard[sample] %>% unlist()
    tmp1 <- combined_ranking_rare_predictions %>% filter(cell_ID == sample)
    
    foreach(alg=unique(tmp1$algorithm), .combine = "rbind", .packages = c("tidyverse","foreach")) %do% {
      
      tmp2 <- tmp1 %>% filter(algorithm == alg)
      if(nrow(tmp2)==0){break}
      max_drivers <- max(tmp2 %>% pull(rank))
      
      print(paste0("Calculating stats for ", sample, "(",which(tmp_samples == sample),"/",length(tmp_samples),")", " targets using ", alg))
      
      foreach(n=1:max_drivers, .combine = "rbind", .packages = c("tidyverse")) %do% {
        
        predicted <- tmp2 %>% filter(rank <= n) %>% pull(target)
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
  write_csv(rare_SL_prediction_stats, paste0("results/CCLE_",network_choice,"/full_results_mann_whitney_rankings.csv"))
}else{
  message("SL Prediction stats file already exists. Reading in previous results.")
  message(paste0("To stop this, remove file: ", "results/CCLE_",network_choice,"/full_results_mann_whitney_rankings.csv"))
  rare_SL_prediction_stats <- read_csv(paste0("results/CCLE_",network_choice,"/full_results_mann_whitney_rankings.csv"))
}

if(celltypes[1] != "ALL"){
  rare_SL_prediction_stats <- rare_SL_prediction_stats %>% filter(lineage %in% celltypes)
}


alg_colours <- read_csv("data/algorithm_colours.csv") %>% deframe()


rare_SL_prediction_stats_summarised <- rare_SL_prediction_stats %>%
  # Combining all badDriver simulations to one mean
  mutate(algorithm=ifelse(str_detect(algorithm,"bad_sim"), "badDriver", algorithm)) %>%
  filter(algorithm %in% alg_of_interest) %>%
  # Only keeping results for samples with > 10 gold standards
  filter(length_gs >= 10) %>%
  group_by(algorithm,n) %>%
  # Only keep measurements where more than 10 cells are available to calculate mean
  filter(n()>=10) %>%
  summarise(
    mean_TP = mean(TP),
    #mean_PPV = mean(PPV),
    mean_precision = mean(precision),
    mean_recall = mean(recall),
    mean_F1 = median(F1),
    #median_TP = median(TP),
    #median_PPV = median(PPV),
    #median_precision = median(precision),
    #median_recall = median(recall),
    #median_F1 = median(F1)
    sample_size = n()
  ) %>%
  pivot_longer(cols = -c(algorithm,n,sample_size), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure, levels = c("mean_TP","mean_recall","mean_precision","mean_F1")),
         algorithm = factor(algorithm, levels = names(alg_colours)),
         sample_size = ifelse(algorithm=="badDriver",length(samples),sample_size))

max_predictions<-50

# Dynamic N_max
gs_lengths <- rare_SL_prediction_stats %>% dplyr::select(cell_ID,length_gs) %>% unique() %>% pull(length_gs)


if(max_predictions=="median gold standards"){
  ## Calculates 2*median number of sensitive genes for cell lines with >3 sensitive genes
  N_max <- gs_lengths[which(gs_lengths  > 3 )] %>% median()*2
} else {
  N_max <- max_predictions
}




ggplot(mapping=aes(x = n, y = value, color = algorithm)) +
  geom_line(data=rare_SL_prediction_stats_summarised %>% filter(algorithm %in% alg_of_interest) %>% filter(n<=N_max) %>% filter(!str_detect(algorithm, "consensus|badDriver")), mapping=aes(alpha = sample_size)) +
  geom_line(data=rare_SL_prediction_stats_summarised %>% filter(algorithm %in% alg_of_interest) %>% filter(n<=N_max) %>% filter(str_detect(algorithm, "consensus|badDriver")), linetype = "dashed") +
  scale_color_manual(values = alg_colours) +
  #scale_linetype_manual(
  #  breaks = names(alg_colours),
  #  values = ifelse(str_detect(names(alg_colours),"consensus|badDriver"),"dashed","solid")
  #) +
  ylab("Average Value") +
  xlab("Number of Predicted Sensitive Genes") +
  guides(colour=guide_legend(title="Algorithm (Colour)"),alpha=guide_legend(title = "Sample Size Remaining (Opacity)")) +
  theme_light() +
  facet_wrap(~measure, nrow = 1, scales = "free")

ggsave(paste0("results/CCLE_",network_choice,"/rare_SL_partner_level_stats_max_",max_predictions,"_",paste(celltypes, collapse = "-"),"_mann_whitney_rank.png"),width = 50, height = 20, units = "cm", dpi = 300)










































################################
#Gene Effect Quantitative Plots
################################

# For the quantitative plots, still onyl plot samples with a decent number of rare-gold standards
quant_plot_samples <- names(rare_gold_standard[lapply(rare_gold_standard, length)>=10])



top_colour = 25

gene_effect_z_scores <- fread("validation_data/gene_effect_z_scores.csv")

if(!file.exists(paste0("results/CCLE_",network_choice,"/full_results_SL_partner_level_gene_effect.csv"))){

gene_effect_results <- foreach(alg=unique(SL_partner_results$algorithm), .combine = "rbind", .export = c("gene_effect_z_scores", "SL_partner_results"), .packages = c("tidyverse","foreach")) %dopar% {
  
  # Get SL-partner level results for specific algorithm
  alg_results <- SL_partner_results %>% 
    filter(algorithm==alg) %>% 
    filter(final_rank <= 100) %>%
    # Just keep the necessary info
    dplyr::select(algorithm,lineage,cell_ID,target,final_rank) %>%
    # Only keep data with > 10 replicates
    group_by(lineage,algorithm,final_rank) %>%
    filter(n()>=10) %>%
    ungroup()
  
  # Combine this information with the gene effect z-scores
  
  alg_results <- alg_results %>%
    left_join(gene_effect_z_scores,
              by = c("target"="gene_ID","cell_ID","lineage"), 
              relationship = "many-to-one")
  
  # Summarise the information down to an average z-score per n for each lineage
  
  #alg_results <- alg_results %>%
  #  group_by(lineage,algorithm,final_rank) %>%
  #  summarise(mean_local_z = mean(local_z_score, na.rm = T), mean_global_z = mean(global_z_score, na.rm = T)) %>%
  #  ungroup()
  
  #alg_results
  
  
  foreach(rank_n=seq(1,100), .combine = "rbind") %do% {
    # Get cumulative mean z-scores for each rank
    
    message(paste0("Getting cumulative gene effect means for ", alg, " top ", rank_n))
    
    rank_n_results <- alg_results %>%
      filter(final_rank <= rank_n) %>%
      group_by(lineage,algorithm) %>%
      summarise(mean_local_z = mean(local_z_score, na.rm = T), 
                mean_global_z = mean(global_z_score, na.rm = T), 
                mean_weighted_z = mean(weighted_average, na.rm = T), 
                .groups = "drop" ) %>%
      ungroup() %>%
      mutate(final_rank=rank_n)
    
    rank_n_results
    
    
  }
  
  
}
write_csv(gene_effect_results, paste0("results/CCLE_",network_choice,"/full_results_SL_partner_level_gene_effect.csv"))
}else{
  gene_effect_results <- read_csv(paste0("results/CCLE_",network_choice,"/full_results_SL_partner_level_gene_effect.csv"))
}



  
plot_data <- gene_effect_results %>%
  # Combining all badDriver simulations to one mean
  mutate(algorithm=ifelse(str_detect(algorithm,"bad_sim"), "badDriver", algorithm)) %>%
  group_by(algorithm,final_rank) %>%
  summarise(mean_local_z=mean(mean_local_z, na.rm = T), mean_global_z=mean(mean_global_z, na.rm = T)) %>%
  ungroup() %>%
  filter(algorithm %in% alg_of_interest) %>%
  mutate(algorithm = factor(algorithm, levels = names(alg_colours))) %>%
  mutate(point_scale=   ifelse(
    final_rank <= top_colour, top_colour + 1 - final_rank, (1 - (   (final_rank - min(final_rank)) / (max(final_rank) - min(final_rank))   ))
  )) %>%
  arrange(desc(final_rank))
  

ggplot(plot_data, aes(x=mean_local_z, y=mean_global_z, colour=point_scale)) +
  geom_point() +
  scale_colour_gradient(high = "red", low = "white", 
                        breaks = seq(top_colour,0), labels = c(seq(1,top_colour),"100 (cap)")) +
  geom_vline(xintercept = 0, colour = "black", alpha = 0.25) +
  geom_hline(yintercept = 0, colour = "black", alpha = 0.25) +
  guides(colour=guide_legend(title="Rank")) +
  labs(x="Cumulative Average Gene Effect (Local Z-score)", y= "Cumulative Average Gene Effect (Global Z-score)") +
  theme(panel.background = element_rect(fill="lightgrey")) +
  facet_wrap(~algorithm)
  
ggsave(paste0("results/CCLE_",network_choice,"/SL_partner_level_gene_effect.png"),width = 30, height = 30, units = "cm", dpi = 300)


lineage_plot_alg <- c("PhenoDriverR")

plot_data <- gene_effect_results %>%
  filter(algorithm==lineage_plot_alg) %>%
  mutate(point_scale=   ifelse(
    final_rank <= top_colour, top_colour + 1 - final_rank, (1 - (   (final_rank - min(final_rank)) / (max(final_rank) - min(final_rank))   ))
  )) %>%
  arrange(desc(final_rank))


ggplot(plot_data, aes(x=mean_local_z, y=mean_global_z, colour=point_scale)) +
  geom_point() +
  scale_colour_gradient(high = "red", low = "white", 
                        breaks = seq(top_colour,0), labels = c(seq(1,top_colour),"100 (cap)")) +
  geom_vline(xintercept = 0, colour = "black", alpha = 0.25) +
  geom_hline(yintercept = 0, colour = "black", alpha = 0.25) +
  guides(colour=guide_legend(title="Rank")) +
  labs(x="Cumulative Average Gene Effect (Local Z-score)", y= "Cumulative Average Gene Effect (Global Z-score)") +
  theme(panel.background = element_rect(fill="lightgrey")) +
  facet_wrap(lineage~algorithm)

ggsave(paste0("results/CCLE_",network_choice,"/SL_partner_level_gene_effect_by_lineage_",paste(lineage_plot_alg, collapse = "_"),".png"),width = 30, height = 30, units = "cm", dpi = 300)



# Stats Plot

stats_top_n <- 5

stats_plot <- SL_partner_results %>%
  filter(final_rank <= stats_top_n) %>%
  left_join(gene_effect_z_scores,
            by = c("target"="gene_ID","cell_ID","lineage"), 
            relationship = "many-to-one") %>%
  group_by(lineage,cell_ID,algorithm) %>%
  summarise(top10_mean_local_z=mean(local_z_score, na.rm = T), .groups = "drop") %>%
  mutate(algorithm=ifelse(str_detect(algorithm,"bad_sim"), "badDriver", algorithm)) %>%
  filter(algorithm %in% alg_of_interest) %>%
  mutate(algorithm = factor(algorithm, levels = names(alg_colours)))


ggplot(stats_plot, aes(x=algorithm, y=top10_mean_local_z, fill=algorithm)) +
  geom_boxplot() +
  scale_fill_manual(breaks = names(alg_colours), values = alg_colours) +
  guides(fill="none") +
  ggtitle(paste0("Mean Tissue-Specific Gene Effect Z-Score of Top ", stats_top_n," Predicted Essential Genes")) +
  labs(y="Mean Tissue-Specific Gene Effect Z-Score", x="Algorithm" ) +
  stat_compare_means(method = "wilcox.test",ref.group = "badDriver", label = "p.signif", hide.ns = F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave(paste0("results/CCLE_",network_choice,"/SL_partner_level_gene_effect",stats_top_n,".png"),width = 50, height = 20, units = "cm", dpi = 300)







## Quantitative Plot With Mann Whitney Ranked Results



top_colour = 25

gene_effect_z_scores <- fread("validation_data/gene_effect_z_scores.csv")

if(!file.exists(paste0("results/CCLE_",network_choice,"/full_results_mann_whitney_ranked_level_gene_effect.csv"))){
  
  gene_effect_results <- foreach(alg=unique(combined_ranking_rare_predictions$algorithm), .combine = "rbind", .export = c("gene_effect_z_scores", "combined_ranking_rare_predictions"), .packages = c("tidyverse","foreach")) %dopar% {
    
    # Get SL-partner level results for specific algorithm
    alg_results <- combined_ranking_rare_predictions %>% 
      filter(algorithm==alg) %>% 
      filter(rank <= 100) %>%
      # Just keep the necessary info
      dplyr::select(algorithm,lineage,cell_ID,target,rank) %>%
      # Only keep data with > 10 replicates
      group_by(lineage,algorithm,rank) %>%
      filter(n()>=10) %>%
      ungroup()
    
    # Combine this information with the gene effect z-scores
    
    alg_results <- alg_results %>%
      left_join(gene_effect_z_scores,
                by = c("target"="gene_ID","cell_ID","lineage"), 
                relationship = "many-to-one")
    
    # Summarise the information down to an average z-score per n for each lineage
    
    #alg_results <- alg_results %>%
    #  group_by(lineage,algorithm,final_rank) %>%
    #  summarise(mean_local_z = mean(local_z_score, na.rm = T), mean_global_z = mean(global_z_score, na.rm = T)) %>%
    #  ungroup()
    
    #alg_results
    
    
    foreach(rank_n=seq(1,100), .combine = "rbind") %do% {
      # Get cumulative mean z-scores for each rank
      
      message(paste0("Getting cumulative gene effect means for ", alg, " top ", rank_n))
      
      rank_n_results <- alg_results %>%
        filter(rank <= rank_n) %>%
        group_by(algorithm) %>%
        summarise(mean_local_z = mean(local_z_score, na.rm = T), 
                  mean_global_z = mean(global_z_score, na.rm = T), 
                  mean_weighted_z = mean(weighted_average, na.rm = T), 
                  .groups = "drop" ) %>%
        ungroup() %>%
        mutate(final_rank=rank_n)
      
      rank_n_results
      
      
    }
    
    
  }
  write_csv(gene_effect_results, paste0("results/CCLE_",network_choice,"/full_results_mann_whitney_ranked_level_gene_effect.csv"))
}else{
  gene_effect_results <- read_csv(paste0("results/CCLE_",network_choice,"/full_results_mann_whitney_ranked_level_gene_effect.csv"))
}





plot_data <- gene_effect_results %>%
  # Combining all badDriver simulations to one mean
  mutate(algorithm=ifelse(str_detect(algorithm,"bad_sim"), "badDriver", algorithm)) %>%
  group_by(algorithm,final_rank) %>%
  summarise(mean_local_z=mean(mean_local_z, na.rm = T), mean_global_z=mean(mean_global_z, na.rm = T)) %>%
  ungroup() %>%
  filter(algorithm %in% alg_of_interest) %>%
  mutate(algorithm = factor(algorithm, levels = names(alg_colours))) %>%
  mutate(point_scale=   ifelse(
    final_rank <= top_colour, top_colour + 1 - final_rank, (1 - (   (final_rank - min(final_rank)) / (max(final_rank) - min(final_rank))   ))
  )) %>%
  arrange(desc(final_rank))


ggplot(plot_data, aes(x=mean_local_z, y=mean_global_z, colour=point_scale)) +
  geom_point() +
  scale_colour_gradient(high = "red", low = "white", 
                        breaks = seq(top_colour,0), labels = c(seq(1,top_colour),"100 (cap)")) +
  geom_vline(xintercept = 0, colour = "black", alpha = 0.25) +
  geom_hline(yintercept = 0, colour = "black", alpha = 0.25) +
  guides(colour=guide_legend(title="Rank")) +
  labs(x="Cumulative Average Gene Effect (Local Z-score)", y= "Cumulative Average Gene Effect (Global Z-score)") +
  theme(panel.background = element_rect(fill="lightgrey")) +
  facet_wrap(~algorithm)

ggsave(paste0("results/CCLE_",network_choice,"/mann_whitney_ranks_gene_effect.png"),width = 30, height = 30, units = "cm", dpi = 300)



## Quantitative Plot With Raw Essentiality Scores

top_colour = 25

#picklesv3 <- fread("data/PICKLESv3/master_table_Avana_Score_TKOv3_BF_Zscore_Expr_LOF_GOF_18959genes_1162cells_lineage_disease.txt") %>%
#  dplyr::select(cell_ID=stripped_cell_line_name, gene_ID=Gene, Avana_BF)
gene_effect_z_scores <- fread("validation_data/gene_effect_z_scores.csv")



if(!file.exists(paste0("results/CCLE_",network_choice,"/full_results_average_essentiality.csv"))){
  
  essentiality_results <- foreach(alg=unique(combined_ranking_rare_predictions$algorithm), .combine = "rbind", .export = c("combined_ranking_rare_predictions"), .packages = c("tidyverse","foreach")) %dopar% {
    
    message(paste0(
      
      "Getting cumulative means for ", 
      alg,
      " (", 
      which(unique(combined_ranking_rare_predictions$algorithm)==alg), 
      "/",
      length(unique(combined_ranking_rare_predictions$algorithm)),
      ")"
      
    ))
    
    
    # Get SL-partner level results for specific algorithm
    alg_results <- combined_ranking_rare_predictions %>% 
      filter(algorithm==alg) %>% 
      filter(rank <= 100) %>%
      # Just keep the necessary info
      dplyr::select(algorithm,lineage,cell_ID,target,rank) %>%
      # Only keep data with > 10 replicates
      group_by(lineage,algorithm,rank) %>%
      filter(n()>=10) %>%
      ungroup()
    
    # Combine this information with the gene gene essentiality and z-scores
    
    alg_results <- alg_results %>%
      #left_join(picklesv3,
      #          by = c("target"="gene_ID","cell_ID"), 
      #          relationship = "many-to-one") %>%
      left_join(gene_effect_z_scores,
                by = c("target"="gene_ID","cell_ID","lineage"), 
                relationship = "many-to-one")
    
    
    foreach(rank_n=seq(1,100), .combine = "rbind") %do% {
      # Get cumulative mean values for each rank
      
     
      
      rank_n_results <- alg_results %>%
        filter(rank <= rank_n) %>%
        group_by(algorithm) %>%
        summarise(mean_local_z = mean(local_z_score, na.rm = T), 
                  mean_global_z = mean(global_z_score, na.rm = T), 
                  mean_weighted_z = mean(weighted_average, na.rm = T), 
                  mean_gene_effect = mean(gene_effect, na.rm = T),
                  .groups = "drop" ) %>%
        ungroup() %>%
        mutate(final_rank=rank_n)
      
      rank_n_results
      
      
    }
    
    
  }
  write_csv(essentiality_results, paste0("results/CCLE_",network_choice,"/full_results_average_essentiality.csv"))
}else{
  essentiality_results <- read_csv(paste0("results/CCLE_",network_choice,"/full_results_average_essentiality.csv"))
}





plot_data <- essentiality_results %>%
  # Combining all badDriver simulations to one mean
  mutate(algorithm=ifelse(str_detect(algorithm,"bad_sim"), "badDriver", algorithm)) %>%
  group_by(algorithm,final_rank) %>%
  summarise(mean_local_z=mean(mean_local_z, na.rm = T), 
            mean_global_z=mean(mean_global_z, na.rm = T),
            mean_weighted_z=mean(mean_weighted_z, na.rm = T),
            mean_gene_effect=mean(mean_gene_effect, na.rm = T)
            ) %>%
  ungroup() %>%
  filter(algorithm %in% alg_of_interest) %>%
  mutate(algorithm = factor(algorithm, levels = names(alg_colours))) %>%
  mutate(point_scale=   ifelse(
    final_rank <= top_colour, top_colour + 1 - final_rank, (1 - (   (final_rank - min(final_rank)) / (max(final_rank) - min(final_rank))   ))
  )) %>%
  arrange(desc(final_rank))


ggplot(plot_data, aes(x=mean_weighted_z, y=mean_gene_effect, colour=point_scale)) +
  geom_point() +
  scale_colour_gradient(high = "red", low = "white", 
                        breaks = seq(top_colour,0), labels = c(seq(1,top_colour),"100 (cap)")) +
  geom_vline(xintercept = 0, colour = "black", alpha = 0.25) +
  geom_hline(yintercept = 0, colour = "black", alpha = 0.25) +
  guides(colour=guide_legend(title="Rank")) +
  labs(x="Cumulative Average of Weighted Z-Scores", y= "Cumulative Average of Raw Gene Effect") +
  theme(panel.background = element_rect(fill="lightgrey")) +
  facet_wrap(~algorithm)

ggsave(paste0("results/CCLE_",network_choice,"/mann_whitney_ranks_essentiality_vs_z.png"),width = 30, height = 30, units = "cm", dpi = 300)



#Top N weighted z score box plot



alg_of_interest <- c(
  "consensus_topklists_median_CSN_NCUA_CSN_MDS_OncoImpact_sysSVM2",
  "consensus_topklists_geomean_DawnRank_OncoImpact_PersonaDrive",
  "consensus_topklists_geomean_ALL",
  "DawnRank",
  "OncoImpact",
  "PNC",
  "PRODIGY",
  "PersonaDrive",
  "SCS",
  "sysSVM2",
  "PhenoDriverR",
  "CSN_MDS",
  "CSN_NCUA",
  "badDriver"
)


alg_colours <- read_csv("data/algorithm_colours.csv") %>% deframe()




topn <- 5

top_predictions_z_scores <- combined_ranking_rare_predictions %>%
  filter(rank <= topn, cell_ID %in% quant_plot_samples) %>%
  left_join(gene_effect_z_scores,by=c("lineage","cell_ID","target"="gene_ID")) %>%
  mutate(algorithm = ifelse(str_detect(algorithm,"bad"), "badDriver", algorithm)) %>%
  filter(algorithm %in% alg_of_interest) %>%
  mutate(algorithm = factor(algorithm, levels=alg_of_interest))

plot_data <- foreach(i=1:topn, .combine = "rbind") %do% {
  top_predictions_z_scores %>%
    filter(rank <= i) %>%
    dplyr::select(algorithm, lineage, cell_ID, target, weighted_average, local_z_score, global_z_score) %>%
    mutate(top_n = factor(as.character(i), levels=c("1","2","3","4","5")))
}

ggplot(plot_data, aes(x=algorithm, y=weighted_average, colour=top_n)) +
  geom_boxplot() +
  coord_cartesian(ylim=c(-1.5,1.5))


plot_data2 <- plot_data %>%
  group_by(algorithm,top_n) %>%
  summarise(mean=mean(weighted_average, na.rm = T), n=n(), s=sd(weighted_average, na.rm=T)) %>%
  mutate(ci95=qt(0.975,df=n-2)*s/sqrt(n))

plot_colours <- colorRampPalette(c("red","black"))

ggplot(plot_data2, aes(x=algorithm, y=mean,colour=top_n)) +
  geom_point(position = position_dodge(width=0.75)) +
  geom_errorbar(aes(ymin=mean-ci95,ymax=mean+ci95),
                position = position_dodge(width=0.75)) +
  scale_colour_manual(breaks = c("1","2","3","4","5"),
                      values = plot_colours(5)) +
  ylab("Mean Weighted Average Gene Effect Z-Score +/- CI95%") +
  xlab("") +
  theme_bw() +
  geom_vline(xintercept = seq(1.5, length(unique(plot_data2$algorithm))-0.5, 1)) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(colour=guide_legend(title = "Top N Predictions"))


ggsave("plots/paper_fig_quantitative_gene_effect2.png", units = "px", width = 5000, height = 2500, dpi = 300)






###################
# Rank Correlation
###################

if(!exists("SL_partner_results")){
  SL_partner_results <- fread(paste0("results/CCLE_",network_choice,"/aggregated_predicted_essential_genes.csv"))
}

if(!exists("gene_effect_z_scores")){
  gene_effect_z_scores <- fread("validation_data/gene_effect_z_scores.csv")
}


picklesv3 <- fread("data/PICKLESv3/master_table_Avana_Score_TKOv3_BF_Zscore_Expr_LOF_GOF_18959genes_1162cells_lineage_disease.txt") %>%
  dplyr::select(cell_ID=stripped_cell_line_name, gene_ID=Gene, Avana_BF)




gs_ranks <- gene_effect_z_scores %>%
  filter(gene_ID %in% genes) %>%
  arrange(global_z_score) %>%
  group_by(cell_ID) %>%
  mutate(gs_global_rank=row_number()) %>%
  ungroup() %>%
  arrange(local_z_score) %>%
  group_by(cell_ID) %>%
  mutate(gs_local_rank=row_number()) %>%
  ungroup() %>%
  # Add picklesv3 scores (not z-scores, approximately the same as raw gene effect)
  left_join(picklesv3, by = c("cell_ID","gene_ID")) %>%
  ungroup() %>%
  arrange(desc(Avana_BF)) %>%
  group_by(cell_ID) %>%
  mutate(Avana_BF_rank=row_number()) %>%
  ungroup()



rank_cor_all_essential <- SL_partner_results %>%
  dplyr::select(lineage,cell_ID,gene_ID=target,predicted_rank=final_rank,algorithm) %>%
  filter(algorithm %in% alg_of_interest) %>%
  filter(predicted_rank <= 100) %>%
  inner_join(gs_ranks, by = c("lineage","cell_ID","gene_ID"))
  


ggplot(rank_cor_all_essential, aes(x=predicted_rank, y=Avana_BF_rank)) +
  geom_pointdensity() +
  facet_wrap(~algorithm)
  












get_consensus <- function(df, algs,measure="geo.mean",title){
  
  consensus_drivers <- foreach(cell=unique(df$cell_ID), .combine = "rbind") %do% {
    
    indiv_ranks <- foreach(alg=algs) %do% {
      x <- df %>% 
        filter(cell_ID==cell) %>%
        filter(algorithm==alg) %>%
        arrange(rank) %>%
        pull(target)
      
      # NA indicates there were no results from that algorithm for that cell_ID. In this case, replace with empty vector
      if(any(is.na(x))){
        x <- c()
      }
      
      x
      
    }
    
    names(indiv_ranks) <- algs
    
    
    lineage <- df %>% filter(cell_ID==cell) %>% pull(lineage) %>% head(1)
    
    topk <- TopKLists::Borda(indiv_ranks)[[1]]
    
    topk %>%
      mutate(lineage=lineage,cell_ID=cell, rank=row_number()) %>%
      pivot_longer(cols = c(mean,median,geo.mean,l2norm), names_to = "algorithm", values_to = "target") %>%
      filter(algorithm==measure) %>%
      mutate(algorithm = title) %>%
      dplyr::select(lineage,cell_ID,target,rank,algorithm)
    
  }
  
}


