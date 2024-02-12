suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(readxl, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(ggpubr, quietly = T))
suppressPackageStartupMessages (library(ggpointdensity, quietly = T))





# Handling input arguments
option_list = list(
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network being used", metavar ="Network")
  
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


network_choice <- opt$network



cell_list <- read_csv("data/cell_list.csv")

#############################
# Get cell line information
#############################

CCLE_sample_info <- read_csv("data/CCLE/Model.csv") %>%
  dplyr::select(DepMap_ID=ModelID, lineage=OncotreeLineage, cell_ID=StrippedCellLineName) %>%
  mutate(lineage=gsub(" |/","_",lineage)) %>%
  filter(cell_ID %in% cell_list$cell_ID) %>%
  unique()



#############################
# Data
#############################

gene_effect <- fread("data/CCLE/CRISPRGeneEffect.csv") %>%
  as.data.frame() %>%
  dplyr::rename(DepMap_ID=V1) %>%
  inner_join(CCLE_sample_info[c("DepMap_ID","cell_ID")], by = "DepMap_ID") %>%
  dplyr::select(-DepMap_ID) %>%
  column_to_rownames("cell_ID")

colnames(gene_effect) <- gsub(" \\(.*","", colnames(gene_effect))
gene_effect <- gene_effect[!duplicated(colnames(gene_effect))]

gene_effect_z_scores <- fread("validation_data/gene_effect_z_scores.csv")


gold_standard <- read_csv(paste0("validation_data/CCLE_",network_choice,"/gold_standards.csv")) %>%
  unique()

rare_gold_standard <- read_csv(paste0("validation_data/CCLE_",network_choice,"/rare_gold_standards.csv")) %>%
  unique()

#############################
# Check Normality
#############################


scaled_gene_effect <- gene_effect %>%
  rownames_to_column("cell_ID") %>%
  pivot_longer(cols = -cell_ID, names_to = "gene_ID", values_to = "gene_effect") %>%
  group_by(gene_ID) %>%
  mutate(scaled_gene_effect= (gene_effect - min(gene_effect))/(max(gene_effect) - min(gene_effect))  )

ggplot(scaled_gene_effect, aes(x=scaled_gene_effect)) +
  geom_histogram(bins = 20) +
  theme_minimal() +
  xlab("Min-Max Scaled Gene Effect") +
  ylab("Count")

ggsave("plots/QC/min_max_scaled_gene_effect_distribution.png", width = 10, height = 10, units = "cm", dpi=300)



ggplot(scaled_gene_effect, aes(x=gene_effect)) +
  geom_density() +
  theme_minimal() +
  xlab("Raw Gene Effect") +
  ylab("Count")



#ggqqplot(scaled_gene_effect$scaled_gene_effect)


#############################
# Check Gold Standards
#############################

dir.create("plots/QC/gold_standard_outlier_examples/local_examples", recursive = T)
dir.create("plots/QC/gold_standard_outlier_examples/global_examples", recursive = T)

if(gs_type == "GESD_gene_effect"){
  


gs_sample <- gold_standard %>%
  filter(gene_ID %in% sample(gold_standard$gene_ID, 12, replace = F)) %>% 
  mutate(is_gold_standard=T)

check_gene_effect_gs <- gene_effect[,unique(gs_sample$gene_ID)] %>% 
  as.data.frame() %>%
  rownames_to_column("cell_ID") %>%
  pivot_longer(-cell_ID, names_to = "gene_ID", values_to = "gene_effect") %>%
  left_join(gs_sample, by = c("cell_ID","gene_ID")) %>%
  mutate(is_gold_standard = ifelse(is.na(is_gold_standard), FALSE, is_gold_standard)) %>%
  left_join(CCLE_sample_info %>% dplyr::select(cell_ID, lineage), by = "cell_ID")

ggplot(check_gene_effect_gs, aes(x=gene_effect)) +
  geom_density() +
  geom_point(aes(y=is_gold_standard, colour = is_gold_standard)) +
  facet_wrap(~gene_ID, scales = "free")

ggsave("plots/QC/eg_gold_standards.png", width = 30, height = 25, units = "cm", dpi=300)



# Select number of plots to make
n = 10

# Local Outliers

# Look only at locally detected outliers
local_eg <- gold_standard %>%
  filter(local_dependent)
# Randomly select 3 gene - lineage combinations
random_selection <- local_eg %>% dplyr::select(lineage,gene_ID) %>% unique() %>% deframe() %>% sample(n)

for(combo in random_selection){
  check_gene_effect_gs <- gene_effect %>% 
    dplyr::select(!!combo) %>%
    rownames_to_column("cell_ID") %>%
    pivot_longer(-cell_ID, names_to = "gene_ID", values_to = "gene_effect") %>%
    left_join(gold_standard %>% dplyr::select(cell_ID, gene_ID, local_dependent), by = c("cell_ID","gene_ID")) %>%
    mutate(local_dependent = ifelse(is.na(local_dependent), FALSE, local_dependent)) %>%
    left_join(CCLE_sample_info %>% dplyr::select(cell_ID, lineage), by = "cell_ID") %>%
    filter(lineage==names(random_selection[which(random_selection==combo)]))
  
  plot <- ggplot(check_gene_effect_gs, aes(x=gene_effect)) +
    geom_density() +
    geom_point(aes(y=local_dependent, colour = local_dependent, size = local_dependent)) +
    scale_colour_manual(breaks = c(TRUE,FALSE), values = c("red","black")) +
    scale_size_manual(breaks = c(TRUE,FALSE), values = c(5,2)) +
    ggtitle(paste0("Example: ", combo, " in ", names(random_selection[which(random_selection==combo)]), " cells")) +
    labs(x = "Raw Gene Effect", y = "Probability Density", caption =  NULL) +
    guides(size="none", colour=guide_legend(title="Local Outlier")) +
    theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),panel.background = element_rect(fill="white"))
    
  print(plot)
  ggsave(paste0("plots/QC/gold_standard_outlier_examples/local_examples/",combo,"_",names(random_selection[which(random_selection==combo)]),".png"), dpi = 100, width = 10, height = 10, units = "cm")
}


# Global Outliers

# Look only at locally detected outliers
random_selection <- gold_standard %>%
  filter(global_dependent) %>%
  pull(gene_ID) %>%
  unique() %>%
  sample(n)

for(gene in random_selection){
  check_gene_effect_gs <- gene_effect %>% 
    dplyr::select(!!gene) %>%
    rownames_to_column("cell_ID") %>%
    pivot_longer(-cell_ID, names_to = "gene_ID", values_to = "gene_effect") %>%
    left_join(gold_standard %>% dplyr::select(cell_ID, gene_ID, global_dependent), by = c("cell_ID","gene_ID")) %>%
    mutate(global_dependent = ifelse(is.na(global_dependent), FALSE, global_dependent))
  
  plot <- ggplot(check_gene_effect_gs, aes(x=gene_effect)) +
    geom_density() +
    geom_point(aes(y=global_dependent, colour = global_dependent, size = global_dependent)) +
    scale_colour_manual(breaks = c(TRUE,FALSE), values = c("red","black")) +
    scale_size_manual(breaks = c(TRUE,FALSE), values = c(5,2)) +
    ggtitle(paste0("Example: ", gene, " in all cells")) +
    labs(x = "Raw Gene Effect", y = "Probability Density", caption =  NULL) +
    guides(size="none", colour=guide_legend(title="Global Outlier")) +
    theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),panel.background = element_rect(fill="white"))
  
  print(plot)
  ggsave(paste0("plots/QC/gold_standard_outlier_examples/global_examples/",gene,".png"), dpi = 100, width = 10, height = 10, units = "cm")
}



test <- gold_standard %>% 
  mutate(cell_gene_combo = paste(cell_ID,gene_ID,sep = "-"))

test <- test %>%
  filter(cell_gene_combo %in% test$cell_gene_combo[which(duplicated(test$cell_gene_combo))])
  
  filter(cell_gene_combo %in% 
           which(gold_standard %>% mutate(cell_gene_combo = paste(cell_ID,gene_ID,sep = "-")) %>% pull(cell_gene_combo) %>% duplicated())
         )
}

if(gs_type == "picklesv3"){
  
  picklesv3 <- fread("data/PICKLESv3/master_table_Avana_Score_TKOv3_BF_Zscore_Expr_LOF_GOF_18959genes_1162cells_lineage_disease.txt")
  
  
  
  ###################
  ##Check Normality##
  ###################
  
  
  scaled_picklesv3 <- picklesv3 %>%
    dplyr::select(gene_ID=Gene,cell_ID=stripped_cell_line_name,Avana_Z,Avana_BF,Avana_Chronos,Score_Z,Score_BF) %>%
    group_by(gene_ID) %>%
    mutate(scaled_Avana_Z= (Avana_Z - min(Avana_Z,na.rm = T))/(max(Avana_Z,na.rm = T) - min(Avana_Z,na.rm = T))  ) %>%
    mutate(scaled_Avana_BF= (Avana_BF - min(Avana_BF,na.rm = T))/(max(Avana_BF,na.rm = T) - min(Avana_BF,na.rm = T))  ) %>%
    mutate(scaled_Avana_Chronos= (Avana_Chronos - min(Avana_Chronos,na.rm = T))/(max(Avana_Chronos,na.rm = T) - min(Avana_Chronos,na.rm = T))  ) %>%
    mutate(scaled_Score_Z= (Score_Z - min(Score_Z,na.rm = T))/(max(Score_Z,na.rm = T) - min(Score_Z,na.rm = T))  ) %>%
    mutate(scaled_Score_BF= (Score_BF - min(Score_BF,na.rm = T))/(max(Score_BF,na.rm = T) - min(Score_BF,na.rm = T))  ) %>%
    pivot_longer(cols = c(scaled_Avana_Z,scaled_Avana_BF,scaled_Avana_Chronos,scaled_Score_Z,scaled_Score_BF), names_to = "measure",values_to = "value")
  
  ggplot(scaled_picklesv3, aes(x=value)) +
    geom_histogram(bins = 20) +
    theme_minimal() +
    xlab("Min-Max Scaled Values") +
    ylab("Count") +
    facet_wrap(~measure)
  
  ggsave("plots/QC/picklesv3_distribution.png", width = 10, height = 10, units = "cm", dpi=300)
  
  
  
  
  
  
  # Select number of plots to make
  n = 10
  
  # Local Outliers
  
  # Look only at locally detected outliers
  local_eg <- rare_gold_standard
  # Randomly select gene - lineage combinations
  random_selection <- local_eg %>% dplyr::select(lineage,gene_ID) %>% unique() %>% deframe() %>% sample(n)
  
  for(combo in random_selection){
    check_gene_effect_gs <- gene_effect %>% 
      # Select the single gene
      dplyr::select(!!combo)
    # Fix colname
    colnames(check_gene_effect_gs) <- combo
    # Join with gold standard info
    check_gene_effect_gs <- check_gene_effect_gs %>%
      rownames_to_column("cell_ID") %>%
      pivot_longer(-cell_ID, names_to = "gene_ID", values_to = "gene_effect") %>%
      left_join(rare_gold_standard %>% dplyr::select(cell_ID, gene_ID) %>% mutate(is_gold_standard = TRUE), by = c("cell_ID","gene_ID")) %>%
      mutate(is_gold_standard = ifelse(is.na(is_gold_standard), FALSE, is_gold_standard)) %>%
      # Keep only cells from same lineage
      left_join(CCLE_sample_info %>% dplyr::select(cell_ID, lineage), by = "cell_ID") %>%
      filter(lineage==names(random_selection[which(random_selection==combo)]))
    
    plot <- ggplot(check_gene_effect_gs, aes(x=gene_effect)) +
      geom_density() +
      geom_point(aes(y=is_gold_standard, colour = is_gold_standard, size = is_gold_standard)) +
      scale_colour_manual(breaks = c(TRUE,FALSE), values = c("red","black")) +
      scale_size_manual(breaks = c(TRUE,FALSE), values = c(5,2)) +
      ggtitle(paste0("Example: ", combo, " in ", names(random_selection[which(random_selection==combo)]), " cells")) +
      labs(x = "Raw Gene Effect", y = "Probability Density", caption =  NULL) +
      guides(size="none", colour=guide_legend(title="Gold Standard")) +
      theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),panel.background = element_rect(fill="white"))
    
    print(plot)
    ggsave(paste0("plots/QC/gold_standard_outlier_examples/local_examples/",combo,"_",names(random_selection[which(random_selection==combo)]),".png"), dpi = 100, width = 10, height = 10, units = "cm")
  }
  
  
  # Global Outliers
  
  random_selection <- gold_standard %>%
    pull(gene_ID) %>%
    unique() %>%
    sample(n)
  
  for(gene in random_selection){
    check_gene_effect_gs <- gene_effect %>% 
      dplyr::select(!!gene) %>%
      rownames_to_column("cell_ID") %>%
      pivot_longer(-cell_ID, names_to = "gene_ID", values_to = "gene_effect") %>%
      left_join(gold_standard %>% dplyr::select(cell_ID, gene_ID) %>% mutate(is_gold_standard = TRUE), by = c("cell_ID","gene_ID")) %>%
      mutate(is_gold_standard = ifelse(is.na(is_gold_standard), FALSE, is_gold_standard))
    
    plot <- ggplot(check_gene_effect_gs, aes(x=gene_effect)) +
      geom_density() +
      geom_point(aes(y=is_gold_standard, colour = is_gold_standard, size = is_gold_standard)) +
      scale_colour_manual(breaks = c(TRUE,FALSE), values = c("red","black")) +
      scale_size_manual(breaks = c(TRUE,FALSE), values = c(5,2)) +
      ggtitle(paste0("Example: ", gene, " in all cells")) +
      labs(x = "Raw Gene Effect", y = "Probability Density", caption =  NULL) +
      guides(size="none", colour=guide_legend(title="Gold Standard")) +
      theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),panel.background = element_rect(fill="white"))
    
    print(plot)
    ggsave(paste0("plots/QC/gold_standard_outlier_examples/global_examples/",gene,".png"), dpi = 100, width = 10, height = 10, units = "cm")
  }
  
}



# Z-score Rank Correlation

qc_z_scores <- gene_effect %>%
  rownames_to_column("cell_ID") %>%
  pivot_longer(-cell_ID, names_to = "gene_ID", values_to = "gene_effect") %>%
  inner_join(gene_effect_z_scores, by = c("cell_ID","gene_ID")) %>%
  na.omit() %>%
  # Raw Ranks
  arrange(gene_effect) %>%
  group_by(cell_ID) %>%
  mutate(raw_rank=row_number()) %>%
  # Global Ranks
  arrange(global_z_score) %>%
  group_by(cell_ID) %>%
  mutate(global_z_rank=row_number())
  

qc_z_scores <- qc_z_scores[sample(1:nrow(qc_z_scores), 10000),]



ggplot(qc_z_scores, aes(x=raw_rank,y=global_z_rank)) +
  geom_pointdensity()


qc_z_scores <- qc_z_scores %>%
  # Add all gold standards
  left_join(gold_standard %>% dplyr::select(cell_ID, gene_ID) %>% mutate(is_gold_standard = TRUE), by = c("cell_ID","gene_ID")) %>%
  mutate(is_gold_standard = ifelse(is.na(is_gold_standard), FALSE, is_gold_standard)) %>%
  # Add rare gold standards
  left_join(rare_gold_standard %>% dplyr::select(cell_ID, gene_ID) %>% mutate(is_rare_gold_standard = TRUE), by = c("cell_ID","gene_ID")) %>%
  mutate(is_rare_gold_standard = ifelse(is.na(is_rare_gold_standard), FALSE, is_rare_gold_standard))


p1 <- ggplot(qc_z_scores, aes(x=gene_effect,colour=is_gold_standard)) +
  geom_density() +
  ggtitle("All Gold Standards vs Raw Gene Effect")

p2 <- ggplot(qc_z_scores, aes(x=global_z_score,colour=is_gold_standard)) +
  geom_density()+
  ggtitle("All Gold Standards vs Gene Effect Global Z-Score (All Cells)")

p3 <- ggplot(qc_z_scores, aes(x=global_z_score,colour=is_rare_gold_standard)) +
  geom_density()+
  ggtitle("Rare Gold Standards vs Gene Effect Global Z-Score (All Cells)")

p4 <- ggplot(qc_z_scores, aes(x=gene_effect,colour=is_rare_gold_standard)) +
  geom_density()+
  ggtitle("Rare Gold Standards vs Raw Gene Effect")

p5 <- ggplot(qc_z_scores, aes(x=local_z_score,colour=is_rare_gold_standard)) +
  geom_density()+
  ggtitle("Rare Gold Standards vs Gene Effect Local Z-Score (Same Tissue)")

ggarrange(p1,p2,p3,p4,p5)
ggsave("plots/QC/Gold_standards_vs_gene_effect.png",height = 20, width = 40, units = "cm",dpi = 300)




picklesv3 <- fread("data/PICKLESv3/master_table_Avana_Score_TKOv3_BF_Zscore_Expr_LOF_GOF_18959genes_1162cells_lineage_disease.txt") %>%
  dplyr::select(cell_ID=stripped_cell_line_name,gene_ID=Gene,Avana_Z,Avana_BF,Avana_Chronos,Score_Z,Score_BF)

picklesv3_vs <- qc_z_scores %>%
  dplyr::select(cell_ID,gene_ID,gene_effect,global_z_score) %>%
  inner_join(picklesv3, by = c("cell_ID","gene_ID")) %>%
  pivot_longer(cols = Avana_Z:Score_BF, names_to = "measure",values_to = "value")

picklesv3_vs <- picklesv3_vs[sample(1:nrow(picklesv3_vs), 10000),]

ggplot(picklesv3_vs, aes(x=gene_effect,y=value)) +
  geom_pointdensity() +
  facet_wrap(~measure, scales = "free")


### Note: CRISPRgeneeffect is different to the picklesv3 Avana Chronos, though it is highly correlated. My guess is the picklesv3 version comes from
### ScreenGeneEffect, that hasn't been batch corrected using Harmonia













ggplot(gene_effect_z_scores %>% filter(gene_ID %in% genes), aes(x=weighted_average)) +
  geom_density()

ggplot(gene_effect_z_scores, aes(x=local_z_score)) +
  geom_density()



