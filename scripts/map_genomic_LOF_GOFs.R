#prepare_LOFGOF_annotations.R


# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))


# Handling input arguments
option_list = list(
  make_option(c("-l", "--LOFGOFfilepath"), type="character", default="data/LOFGOF/logofunc-predictions/LoGoFuncVotingEnsemble_preds_final.csv", 
              help="Path to LoGoFunc Prediction file", metavar ="LOFGOF File Path"),
  make_option(c("-m", "--mutationfilepath"), type="character", default="data/CCLE/OmicsSomaticMutations.csv", 
              help="Path to mutations file", metavar ="Mutations File Path"),
  make_option(c("-o", "--outputfilepath"), type="character", default="data/LOFGOF/annotated_CCLE_mutations.csv", 
              help="Path to output file", metavar ="Output File Path"),
  make_option(c("-t", "--threads"), type="character", default=4, 
              help="Number of threads to use", metavar ="Threads")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

LOFGOFfilepath <- opt$LOFGOFfilepath
mutationfilepath <- opt$mutationfilepath
outputfilepath <- opt$outputfilepath


LOF_GOF_predictions <- data.table::fread(LOFGOFfilepath, sep = "\t",
                                col.names =  c("chr", "pos", "ref", "alt", "id", "prediction", "score_neutral", "score_GOF", "score_LOF"),
                                colClasses =  "cncccfnnn") %>%
  mutate(chr = paste0("chr",chr))


annotated_mutations <- fread(mutationfilepath) %>%
  left_join(LOF_GOF_predictions, 
                       by = c("Chrom" = "chr", 
                              "Pos" = "pos",
                              "Ref" = "ref",
                              "Alt" = "alt"))

write_csv(annotated_mutations,outputfilepath)
