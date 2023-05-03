#!/usr/bin/env Rscript --vanilla

print("running test")
print(paste0("Current directory is", getwd()))

data_dir <- "data/"

read.csv(paste0(data_dir,"cell_list.csv"))