#evaluate_drug_sensitivity.r


# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(foreach, quietly = T))
suppressPackageStartupMessages (library(doParallel, quietly = T))
suppressPackageStartupMessages (library(ggpubr, quietly = T))
suppressPackageStartupMessages (library(coin, quietly = T))




# Handling input arguments
option_list = list(
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network being used", metavar ="Network"),
  make_option(c("-a", "--algorithms"), type="character", default="ALL", 
              help="algorithms to include in comparison separated by spaces, or 'ALL' (Default), or lead with '-' to exclude", metavar ="Algorithms"),
  make_option(c("-t", "--threads"), type="integer", default=4, 
              help="Number of threads to use if running in parallel", metavar ="Threads"),
  make_option(c("-c", "--celltype"), type="character", default="ALL", 
              help="Cell types to include in the analysis separated by semicolons, or 'ALL' (Default)", metavar ="badDriver")
  
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

network_choice <- opt$network
algorithms <- opt$algorithms
threads <- opt$threads
celltypes <- opt$celltype


if(celltypes != "ALL"){
  celltypes <- str_split_1(celltypes, pattern = ";")
}


threads <- 4


if(threads>1){
  #registerDoParallel(cores=threads)
  cl <- makeCluster(threads, outfile = "log/evaluate_drug_sensitivity.log")
  registerDoParallel(cl)
}

#############################
# Sample Info
#############################

sample_info <- read_csv(paste0("validation_data/CCLE_",network_choice,"/sample_info.csv"))
samples <- sample_info$cell_ID %>% sort()




#############################
# Read In Gold Standards
#############################

gold_standard_drug_sensitivity <- read_csv(paste0("validation_data/CCLE_",network_choice,"/gold_standard_drug_sensitivity.csv"))

all_gold_standards <- gold_standard_drug_sensitivity %>%
  filter(global_sensitive) %>%
  dplyr::select(cell_ID,drug_ID) %>%
  unique() %>%
  group_by(cell_ID) %>%
  summarise(sensitive_drugs = list(drug_ID)) %>%
  deframe()

rare_gold_standards <- gold_standard_drug_sensitivity %>%
  filter(rare_sensitive) %>%
  dplyr::select(cell_ID,drug_ID) %>%
  unique() %>%
  group_by(cell_ID) %>%
  summarise(sensitive_drugs = list(drug_ID)) %>%
  deframe()


drug_targets <- read_csv(paste0("validation_data/CCLE_",network_choice,"/drug_targets.csv")) %>%
  dplyr::select(drug_ID,gene_ID) %>%
  unique() %>%
  # Get the number of gene targets for each drug
  group_by(drug_ID) %>%
  mutate(n_targets=length(gene_ID)) %>%
  ungroup()



###################################
# Read In Predicted Essential Genes
###################################


all_essential_genes <- read_csv(paste0("results/CCLE_",network_choice,"/aggregated_predicted_essential_genes.csv"))

rare_essential_genes <- read_csv(paste0("results/CCLE_",network_choice,"/rare_essential_genes_mann_whitney.csv")) %>%
  filter(padj < 0.05) %>%
  group_by(algorithm,lineage,cell_ID) %>%
  arrange(rank) %>%
  mutate(rank=row_number()) %>%
  ungroup() %>%
  arrange(algorithm,lineage,cell_ID,rank)


####################################
# All Predictions
####################################

all_predictions <- all_essential_genes %>%
  dplyr::select(algorithm,lineage,cell_ID,gene_ID=target,gene_rank=final_rank) %>%
  left_join(drug_targets, by = "gene_ID", relationship="many-to-many") %>%
  filter(!is.na(drug_ID)) %>%
  dplyr::select(-gene_ID) %>%
  group_by(algorithm,lineage,cell_ID,drug_ID,n_targets) %>%
  summarise(min_gene_rank=min(gene_rank)) %>%
  ungroup() %>%
  group_by(algorithm,cell_ID) %>%
  # prioritise drugs that target fewer genes
  arrange(min_gene_rank,n_targets) %>%
  filter(!duplicated(drug_ID)) %>%
  mutate(rank=row_number()) %>%
  ungroup() %>%
  arrange(algorithm,lineage,cell_ID,rank)


rare_predictions <- rare_essential_genes %>%
  dplyr::select(algorithm,lineage,cell_ID,gene_ID=target,gene_rank=rank) %>%
  left_join(drug_targets, by = "gene_ID", relationship="many-to-many") %>%
  filter(!is.na(drug_ID)) %>%
  dplyr::select(-gene_ID) %>%
  group_by(algorithm,lineage,cell_ID,drug_ID,n_targets) %>%
  summarise(min_gene_rank=min(gene_rank)) %>%
  ungroup() %>%
  group_by(algorithm,cell_ID) %>%
  # prioritise drugs that target fewer genes
  arrange(min_gene_rank,n_targets) %>%
  filter(!duplicated(drug_ID)) %>%
  mutate(rank=row_number()) %>%
  ungroup() %>%
  arrange(algorithm,lineage,cell_ID,rank)
  
  
  

####################################
# BadDrug Predictions
####################################


suppressWarnings(rm(badDrug_predictions))

for(ct in list.dirs(paste0("bad_driver_simulations/CCLE_",network_choice,"_drugs/"), recursive = F)){
  ct <- gsub(".*/","",ct)
  message("badDrug results for ", ct)
  
  for(run in list.files(paste0("bad_driver_simulations/CCLE_",network_choice,"_drugs/",ct))) {
    alg <- gsub("_[0-9]+.csv","",run)
    run <- str_match(run,"bad_drug_([0-9]+).csv")[2]
    
    indiv_result <- read_csv(paste0("bad_driver_simulations/CCLE_",network_choice,"_drugs/",ct,"/", alg, "_", run, ".csv"), col_types = cols()) %>%
      mutate(algorithm=paste0(alg,"_",run))
    
    if(!exists("badDrug_predictions")){
      badDrug_predictions <- indiv_result
    }else{
      badDrug_predictions <- rbind(badDrug_predictions,indiv_result)
    }
    
  }
  
}

badDrug_predictions <- badDrug_predictions %>% filter(rank <= 100)


all_predictions <- all_predictions %>% rbind(badDrug_predictions %>% mutate(n_targets=NA,min_gene_rank=NA))

rare_predictions <- rare_predictions %>% rbind(badDrug_predictions %>% mutate(n_targets=NA,min_gene_rank=NA))

####################################
# GoodDrug Predictions
####################################

#good_drug_predictions <- read_csv("validation_data/drug_sensitivity.csv") %>%
#  filter(sensitivity_source=="GDSC2") %>%
#  group_by(cell_ID)
  




####################################
# Precily Predictions
####################################

# Precily: Deep neural network based framework for prediction of drug response in both in vivo and in vitro settings
# https://github.com/SmritiChawla/Precily
# The code from fig1d was run separately in a linux environment on the Test_Set.csv file provided
# The cell ID, drug name and predicted drug sensitivities were extracted


## Try out:

### Compare results when only filtered to the same cells and the same drugs
### Compare results when including all drugs from my pipeline, still only the same cells


# Rank precily predictions based on predicted IC50 (false) or predicted IC50 z-scores (true)?
precily_z=T

if(precily_z){
  
  
  get_z_scores <- function(Matrix){
    Sd   <- apply(Matrix, 2, function(x) sd(x, na.rm=T))
    Mean <- apply(Matrix, 2, function(x) mean(x, na.rm=T))
    centered <- suppressWarnings(t(t(Matrix) - Mean))
    z_scores <- suppressWarnings(t(t(centered) / Sd))
    return(z_scores)
  }
  
  precily_fig1d <- read_csv("data/precily_fig1d.csv") %>%
    dplyr::select(cell_ID,drug_ID,predicted_lnIC50) %>%
    pivot_wider(names_from = drug_ID, values_from = predicted_lnIC50) %>%
    column_to_rownames("cell_ID") %>%
    get_z_scores() %>%
    as.data.frame() %>%
    rownames_to_column("cell_ID") %>%
    pivot_longer(cols = -cell_ID, names_to = "drug_ID", values_to = "predicted_z_score") %>%
    group_by(cell_ID) %>%
    arrange(predicted_z_score) %>%
    mutate(rank=row_number()) %>%
    ungroup() %>%
    mutate(algorithm="Precily") %>%
    inner_join(sample_info %>% dplyr::select(cell_ID,lineage), by="cell_ID") %>%
    dplyr::select(-predicted_z_score)
  
}else{
  precily_fig1d <- read_csv("data/precily_fig1d.csv") %>%
    dplyr::select(cell_ID,drug_ID,predicted_lnIC50) %>%
    group_by(cell_ID) %>%
    arrange(predicted_lnIC50) %>%
    mutate(rank=row_number()) %>%
    ungroup() %>%
    mutate(algorithm="Precily") %>%
    inner_join(sample_info %>% dplyr::select(cell_ID,lineage), by="cell_ID") %>%
    dplyr::select(-predicted_lnIC50)
}







all_predictions <- all_predictions %>% rbind(precily_fig1d %>% mutate(n_targets=NA,min_gene_rank=NA))
rare_predictions <- rare_predictions %>% rbind(precily_fig1d %>% mutate(n_targets=NA,min_gene_rank=NA))

precily_samples <- unique(precily_fig1d$cell_ID)

precily_drugs <- unique(precily_fig1d$drug_ID)

length(precily_drugs)
length(precily_drugs[which(precily_drugs%in% drug_targets$drug_ID)])
length(precily_drugs[which(precily_drugs%in% unlist(all_gold_standards))])




if(celltypes[1] != "ALL"){
  all_predictions <- all_predictions %>% filter(lineage %in% celltypes)
}





########################
# All vs All Evaluation
######################## 

tmp_samples <- intersect(all_predictions$cell_ID, names(all_gold_standards))

if(!file.exists(paste0("results/CCLE_",network_choice,"/full_results_all_vs_all_drug_sensitivity.csv"))){
  
  drug_prediction_stats <- foreach(sample=tmp_samples, .combine = "rbind", .packages = c("tidyverse","foreach")) %dopar% {
    
    lineage <- sample_info %>% filter(cell_ID==sample) %>% pull("lineage")
    gs <- all_gold_standards[sample] %>% unlist()
    tmp1 <- all_predictions %>% filter(cell_ID == sample)
    
    foreach(alg=unique(tmp1$algorithm), .combine = "rbind", .packages = c("tidyverse","foreach")) %do% {
      
      tmp2 <- tmp1 %>% filter(algorithm == alg)
      if(nrow(tmp2)==0){break}
      max_drugs <- max(tmp2 %>% pull(rank))
      
      print(paste0("Calculating stats for ", sample, "(",which(tmp_samples == sample),"/",length(tmp_samples),")", " targets using ", alg))
      
      foreach(n=1:max_drugs, .combine = "rbind", .packages = c("tidyverse")) %do% {
        
        predicted <- tmp2 %>% filter(rank <= n) %>% pull(drug_ID)
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
  write_csv(drug_prediction_stats, paste0("results/CCLE_",network_choice,"/full_results_all_vs_all_drug_sensitivity.csv"))
}else{
  message("Drug Prediction stats file already exists. Reading in previous results.")
  message(paste0("To stop this, remove file: ", "results/CCLE_",network_choice,"/full_results_all_vs_all_drug_sensitivity.csv"))
  drug_prediction_stats <- fread(paste0("results/CCLE_",network_choice,"/full_results_all_vs_all_drug_sensitivity.csv"))
}




alg_of_interest <- c(
  #"consensus_topklists_geomean_CSN_NCUA_PersonaDrive_PRODIGY_OncoImpact_sysSVM2_PhenoDriverR", # best for STRINGv11 reference drivers
  #"consensus_topklists_median_SCS_DawnRank_sysSVM2",
  #"consensus_topklists_geomean_CSN_NCUA_OncoImpact_sysSVM2",
  
  #"consensus_topklists_geomean_CSN_NCUA_OncoImpact_sysSVM2", # Best for STRINGv11 SL partners
  #"consensus_topklists_median_SCS_DawnRank_sysSVM2", # Best for STRINGv11 rare SL partners

  "consensus_topklists_median_CSN_NCUA_CSN_MDS_OncoImpact_sysSVM2",
  "consensus_topklists_geomean_DawnRank_OncoImpact_PersonaDrive",
  "consensus_topklists_geomean_ALL",
  
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
  #"CSN_MDS",
  #"CSN_MMS",
  "CSN_NCUA",
  #"LIONESS_DFVS",
  #"LIONESS_MDS",
  #"LIONESS_MMS",
  #"LIONESS_NCUA",
  #"SPCC_DFVS",
  #"SPCC_MDS",
  #"SPCC_MMS",
  "SPCC_NCUA",
  #"SSN_DFVS",
  #"SSN_MDS",
  #"SSN_MMS",
  #"SSN_NCUA",
  
  "badDriver",
  "badDrug",
  "Precily"
)


alg_colours <- read_csv("data/algorithm_colours.csv") %>% deframe()


drug_prediction_stats_summarised <- drug_prediction_stats %>%
  # Combining all badDriver simulations to one mean
  mutate(algorithm=ifelse(str_detect(algorithm,"bad_sim"), "badDriver", algorithm)) %>%
  mutate(algorithm=ifelse(str_detect(algorithm,"bad_drug"), "badDrug", algorithm)) %>%
  filter(algorithm %in% alg_of_interest) %>%
  # Only keeping results for samples with > 10 gold standards
  filter(length_gs >= 10) %>%
  group_by(algorithm,n) %>%
  # Only keep measurements where more than 5 cells are available to calculate mean
  filter(n()>=5) %>%
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
gs_lengths <- drug_prediction_stats %>% dplyr::select(cell_ID,length_gs) %>% unique() %>% pull(length_gs)


if(max_predictions=="median gold standards"){
  ## Calculates 2*median number of sensitive genes for cell lines with >3 sensitive genes
  N_max <- gs_lengths[which(gs_lengths  > 3 )] %>% median()*2
} else {
  N_max <- max_predictions
}




ggplot(mapping=aes(x = n, y = value, color = algorithm)) +
  geom_line(data=drug_prediction_stats_summarised %>% filter(algorithm %in% alg_of_interest) %>% filter(n<=N_max) %>% filter(!str_detect(algorithm, "consensus|bad|Precily")), mapping=aes(alpha = sample_size)) +
  geom_line(data=drug_prediction_stats_summarised %>% filter(algorithm %in% alg_of_interest) %>% filter(n<=N_max) %>% filter(str_detect(algorithm, "consensus|bad|Precily")), linetype = "dashed") +
  scale_color_manual(values = alg_colours) +
  #scale_linetype_manual(
  #  breaks = names(alg_colours),
  #  values = ifelse(str_detect(names(alg_colours),"consensus|badDriver"),"dashed","solid")
  #) +
  ylab("Average Value") +
  xlab("Number of Predicted Sensitive Drugs") +
  guides(colour=guide_legend(title="Algorithm (Colour)"),alpha=guide_legend(title = "Sample Size Remaining (Opacity)")) +
  theme_light() +
  facet_wrap(~measure, nrow = 1, scales = "free")

ggsave(paste0("results/CCLE_",network_choice,"/drug_sensitivity_stats_all_vs_all_max_",max_predictions,"_",paste(celltypes, collapse = "-"),".png"),width = 50, height = 20, units = "cm", dpi = 300)



######### Precily Samples Only Plot #############

precily_plot_all_drugs_data <- drug_prediction_stats %>%
  filter(cell_ID %in% precily_samples) %>%
  # Combining all badDriver simulations to one mean
  mutate(algorithm=ifelse(str_detect(algorithm,"bad_sim"), "badDriver", algorithm)) %>%
  mutate(algorithm=ifelse(str_detect(algorithm,"bad_drug"), "badDrug", algorithm)) %>%
  filter(algorithm %in% alg_of_interest) %>%
  # Only keeping results for samples with > 10 gold standards
  #filter(length_gs >= 10) %>%
  group_by(algorithm,n) %>%
  # Only keep measurements where more than 5 cells are available to calculate mean
  filter(n()>=5) %>%
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


ggplot(mapping=aes(x = n, y = value, color = algorithm)) +
  geom_line(data=precily_plot_all_drugs_data %>% filter(algorithm %in% alg_of_interest) %>% filter(n<=N_max) %>% filter(!str_detect(algorithm, "consensus|bad|Precily")), mapping=aes(alpha = sample_size)) +
  geom_line(data=precily_plot_all_drugs_data %>% filter(algorithm %in% alg_of_interest) %>% filter(n<=N_max) %>% filter(str_detect(algorithm, "consensus|bad|Precily")), linetype = "dashed") +
  scale_color_manual(values = alg_colours) +
  #scale_linetype_manual(
  #  breaks = names(alg_colours),
  #  values = ifelse(str_detect(names(alg_colours),"consensus|badDriver"),"dashed","solid")
  #) +
  ylab("Average Value") +
  xlab("Number of Predicted Sensitive Drugs") +
  guides(colour=guide_legend(title="Algorithm (Colour)"),alpha=guide_legend(title = "Sample Size Remaining (Opacity)")) +
  theme_light() +
  facet_wrap(~measure, nrow = 1, scales = "free")

ggsave(paste0("results/CCLE_",network_choice,"/drug_sensitivity_stats_precily_samples_all_drugs_all_vs_all_max_",max_predictions,"_",paste(celltypes, collapse = "-"),".png"),width = 50, height = 20, units = "cm", dpi = 300)

########################
# All vs Rare Evaluation
######################## 

tmp_samples <- intersect(all_predictions$cell_ID, names(rare_gold_standards))

if(!file.exists(paste0("results/CCLE_",network_choice,"/full_results_all_vs_rare_drug_sensitivity.csv"))){
  
  all_v_rare_drug_prediction_stats <- foreach(sample=tmp_samples, .combine = "rbind", .packages = c("tidyverse","foreach")) %dopar% {
    
    lineage <- sample_info %>% filter(cell_ID==sample) %>% pull("lineage")
    gs <- rare_gold_standards[sample] %>% unlist()
    tmp1 <- all_predictions %>% filter(cell_ID == sample)
    
    foreach(alg=unique(tmp1$algorithm), .combine = "rbind", .packages = c("tidyverse","foreach")) %do% {
      
      tmp2 <- tmp1 %>% filter(algorithm == alg)
      if(nrow(tmp2)==0){break}
      max_drugs <- max(tmp2 %>% pull(rank))
      
      print(paste0("Calculating stats for ", sample, "(",which(tmp_samples == sample),"/",length(tmp_samples),")", " targets using ", alg))
      
      foreach(n=1:max_drugs, .combine = "rbind", .packages = c("tidyverse")) %do% {
        
        predicted <- tmp2 %>% filter(rank <= n) %>% pull(drug_ID)
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
  write_csv(all_v_rare_drug_prediction_stats, paste0("results/CCLE_",network_choice,"/full_results_all_vs_rare_drug_sensitivity.csv"))
}else{
  message("Drug Prediction stats file already exists. Reading in previous results.")
  message(paste0("To stop this, remove file: ", "results/CCLE_",network_choice,"/full_results_all_vs_rare_drug_sensitivity.csv"))
  all_v_rare_drug_prediction_stats <- fread(paste0("results/CCLE_",network_choice,"/full_results_all_vs_rare_drug_sensitivity.csv"))
}






alg_colours <- read_csv("data/algorithm_colours.csv") %>% deframe()


all_v_rare_drug_prediction_stats_summarised <- all_v_rare_drug_prediction_stats %>%
  # Combining all badDriver simulations to one mean
  mutate(algorithm=ifelse(str_detect(algorithm,"bad_sim"), "badDriver", algorithm)) %>%
  mutate(algorithm=ifelse(str_detect(algorithm,"bad_drug"), "badDrug", algorithm)) %>%
  filter(algorithm %in% alg_of_interest) %>%
  # Only keeping results for samples with > 10 gold standards
  #filter(length_gs >= 10) %>%
  group_by(algorithm,n) %>%
  # Only keep measurements where more than 5 cells are available to calculate mean
  filter(n()>=5) %>%
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
gs_lengths <- drug_prediction_stats %>% dplyr::select(cell_ID,length_gs) %>% unique() %>% pull(length_gs)


if(max_predictions=="median gold standards"){
  ## Calculates 2*median number of sensitive genes for cell lines with >3 sensitive genes
  N_max <- gs_lengths[which(gs_lengths  > 3 )] %>% median()*2
} else {
  N_max <- max_predictions
}




ggplot(mapping=aes(x = n, y = value, color = algorithm)) +
  geom_line(data=all_v_rare_drug_prediction_stats_summarised %>% filter(algorithm %in% alg_of_interest) %>% filter(n<=N_max) %>% filter(!str_detect(algorithm, "consensus|bad|Precily")), mapping=aes(alpha = sample_size)) +
  geom_line(data=all_v_rare_drug_prediction_stats_summarised %>% filter(algorithm %in% alg_of_interest) %>% filter(n<=N_max) %>% filter(str_detect(algorithm, "consensus|bad|Precily")), linetype = "dashed") +
  scale_color_manual(values = alg_colours) +
  #scale_linetype_manual(
  #  breaks = names(alg_colours),
  #  values = ifelse(str_detect(names(alg_colours),"consensus|badDriver"),"dashed","solid")
  #) +
  ylab("Average Value") +
  xlab("Number of Predicted Sensitive Drugs") +
  guides(colour=guide_legend(title="Algorithm (Colour)"),alpha=guide_legend(title = "Sample Size Remaining (Opacity)")) +
  theme_light() +
  facet_wrap(~measure, nrow = 1, scales = "free")

ggsave(paste0("results/CCLE_",network_choice,"/drug_sensitivity_stats_all_vs_rare_max_",max_predictions,"_",paste(celltypes, collapse = "-"),".png"),width = 50, height = 20, units = "cm", dpi = 300)



######### Precily Samples Only Plot #############

precily_plot_rare_drugs_data <- all_v_rare_drug_prediction_stats %>%
  filter(cell_ID %in% precily_samples) %>%
  # Combining all badDriver simulations to one mean
  mutate(algorithm=ifelse(str_detect(algorithm,"bad_sim"), "badDriver", algorithm)) %>%
  mutate(algorithm=ifelse(str_detect(algorithm,"bad_drug"), "badDrug", algorithm)) %>%
  filter(algorithm %in% alg_of_interest) %>%
  # Only keeping results for samples with > 10 gold standards
  #filter(length_gs >= 10) %>%
  group_by(algorithm,n) %>%
  # Only keep measurements where more than 5 cells are available to calculate mean
  filter(n()>=5) %>%
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


ggplot(mapping=aes(x = n, y = value, color = algorithm)) +
  geom_line(data=precily_plot_rare_drugs_data %>% filter(algorithm %in% alg_of_interest) %>% filter(n<=N_max) %>% filter(!str_detect(algorithm, "consensus|bad|Precily")), mapping=aes(alpha = sample_size)) +
  geom_line(data=precily_plot_rare_drugs_data %>% filter(algorithm %in% alg_of_interest) %>% filter(n<=N_max) %>% filter(str_detect(algorithm, "consensus|bad|Precily")), linetype = "dashed") +
  scale_color_manual(values = alg_colours) +
  #scale_linetype_manual(
  #  breaks = names(alg_colours),
  #  values = ifelse(str_detect(names(alg_colours),"consensus|badDriver"),"dashed","solid")
  #) +
  ylab("Average Value") +
  xlab("Number of Predicted Sensitive Drugs") +
  guides(colour=guide_legend(title="Algorithm (Colour)"),alpha=guide_legend(title = "Sample Size Remaining (Opacity)")) +
  theme_light() +
  facet_wrap(~measure, nrow = 1, scales = "free")

ggsave(paste0("results/CCLE_",network_choice,"/drug_sensitivity_stats_precily_samples_all_drugs_all_vs_rare_max_",max_predictions,"_",paste(celltypes, collapse = "-"),".png"),width = 50, height = 20, units = "cm", dpi = 300)







######################
#Rare vs Rare
######################



tmp_samples <- intersect(rare_predictions$cell_ID, names(rare_gold_standards))

if(!file.exists(paste0("results/CCLE_",network_choice,"/full_results_rare_vs_rare_drug_sensitivity.csv"))){
  
  rare_v_rare_drug_prediction_stats <- foreach(sample=tmp_samples, .combine = "rbind", .packages = c("tidyverse","foreach")) %dopar% {
    
    lineage <- sample_info %>% filter(cell_ID==sample) %>% pull("lineage")
    gs <- rare_gold_standards[sample] %>% unlist()
    tmp1 <- rare_predictions %>% filter(cell_ID == sample)
    
    foreach(alg=unique(tmp1$algorithm), .combine = "rbind", .packages = c("tidyverse","foreach")) %do% {
      
      tmp2 <- tmp1 %>% filter(algorithm == alg)
      if(nrow(tmp2)==0){break}
      max_drugs <- max(tmp2 %>% pull(rank))
      
      print(paste0("Calculating stats for ", sample, "(",which(tmp_samples == sample),"/",length(tmp_samples),")", " targets using ", alg))
      
      foreach(n=1:max_drugs, .combine = "rbind", .packages = c("tidyverse")) %do% {
        
        predicted <- tmp2 %>% filter(rank <= n) %>% pull(drug_ID)
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
  write_csv(rare_v_rare_drug_prediction_stats, paste0("results/CCLE_",network_choice,"/full_results_rare_vs_rare_drug_sensitivity.csv"))
}else{
  message("Drug Prediction stats file already exists. Reading in previous results.")
  message(paste0("To stop this, remove file: ", "results/CCLE_",network_choice,"/full_results_rare_vs_rare_drug_sensitivity.csv"))
  rare_v_rare_drug_prediction_stats <- fread(paste0("results/CCLE_",network_choice,"/full_results_rare_vs_rare_drug_sensitivity.csv"))
}




alg_colours <- read_csv("data/algorithm_colours.csv") %>% deframe()


rare_v_rare_drug_prediction_stats_summarised <- rare_v_rare_drug_prediction_stats %>%
  # Combining all badDriver simulations to one mean
  mutate(algorithm=ifelse(str_detect(algorithm,"bad_sim"), "badDriver", algorithm)) %>%
  mutate(algorithm=ifelse(str_detect(algorithm,"bad_drug"), "badDrug", algorithm)) %>%
  filter(algorithm %in% alg_of_interest) %>%
  # Only keeping results for samples with > 10 gold standards
  #filter(length_gs >= 10) %>%
  group_by(algorithm,n) %>%
  # Only keep measurements where more than 5 cells are available to calculate mean
  filter(n()>=5) %>%
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
gs_lengths <- rare_v_rare_drug_prediction_stats %>% dplyr::select(cell_ID,length_gs) %>% unique() %>% pull(length_gs)


if(max_predictions=="median gold standards"){
  ## Calculates 2*median number of sensitive genes for cell lines with >3 sensitive genes
  N_max <- gs_lengths[which(gs_lengths  > 3 )] %>% median()*2
} else {
  N_max <- max_predictions
}




ggplot(mapping=aes(x = n, y = value, color = algorithm)) +
  geom_line(data=rare_v_rare_drug_prediction_stats_summarised %>% filter(algorithm %in% alg_of_interest) %>% filter(n<=N_max) %>% filter(!str_detect(algorithm, "consensus|bad|Precily")), mapping=aes(alpha = sample_size)) +
  geom_line(data=rare_v_rare_drug_prediction_stats_summarised %>% filter(algorithm %in% alg_of_interest) %>% filter(n<=N_max) %>% filter(str_detect(algorithm, "consensus|bad|Precily")), linetype = "dashed") +
  scale_color_manual(values = alg_colours) +
  #scale_linetype_manual(
  #  breaks = names(alg_colours),
  #  values = ifelse(str_detect(names(alg_colours),"consensus|badDriver"),"dashed","solid")
  #) +
  ylab("Average Value") +
  xlab("Number of Predicted Sensitive Drugs") +
  guides(colour=guide_legend(title="Algorithm (Colour)"),alpha=guide_legend(title = "Sample Size Remaining (Opacity)")) +
  theme_light() +
  facet_wrap(~measure, nrow = 1, scales = "free")

ggsave(paste0("results/CCLE_",network_choice,"/drug_sensitivity_stats_rare_vs_rare_max_",max_predictions,"_",paste(celltypes, collapse = "-"),".png"),width = 50, height = 20, units = "cm", dpi = 300)



######### Precily Samples Only Plot #############

precily_plot_rare_vs_rare_drugs_data <- rare_v_rare_drug_prediction_stats %>%
  filter(cell_ID %in% precily_samples) %>%
  # Combining all badDriver simulations to one mean
  mutate(algorithm=ifelse(str_detect(algorithm,"bad_sim"), "badDriver", algorithm)) %>%
  mutate(algorithm=ifelse(str_detect(algorithm,"bad_drug"), "badDrug", algorithm)) %>%
  filter(algorithm %in% alg_of_interest) %>%
  # Only keeping results for samples with > 10 gold standards
  #filter(length_gs >= 10) %>%
  group_by(algorithm,n) %>%
  # Only keep measurements where more than 5 cells are available to calculate mean
  filter(n()>=5) %>%
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


ggplot(mapping=aes(x = n, y = value, color = algorithm)) +
  geom_line(data=precily_plot_rare_vs_rare_drugs_data %>% filter(algorithm %in% alg_of_interest) %>% filter(n<=N_max) %>% filter(!str_detect(algorithm, "consensus|bad|Precily")), mapping=aes(alpha = sample_size)) +
  geom_line(data=precily_plot_rare_vs_rare_drugs_data %>% filter(algorithm %in% alg_of_interest) %>% filter(n<=N_max) %>% filter(str_detect(algorithm, "consensus|bad|Precily")), linetype = "dashed") +
  scale_color_manual(values = alg_colours) +
  #scale_linetype_manual(
  #  breaks = names(alg_colours),
  #  values = ifelse(str_detect(names(alg_colours),"consensus|badDriver"),"dashed","solid")
  #) +
  ylab("Average Value") +
  xlab("Number of Predicted Sensitive Drugs") +
  guides(colour=guide_legend(title="Algorithm (Colour)"),alpha=guide_legend(title = "Sample Size Remaining (Opacity)")) +
  theme_light() +
  facet_wrap(~measure, nrow = 1, scales = "free")

ggsave(paste0("results/CCLE_",network_choice,"/drug_sensitivity_stats_precily_samples_all_drugs_rare_vs_rare_max_",max_predictions,"_",paste(celltypes, collapse = "-"),".png"),width = 50, height = 20, units = "cm", dpi = 300)









