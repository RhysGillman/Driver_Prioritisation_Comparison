# run_ANNOVAR.r

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))

# Handling input arguments
option_list = list(
  make_option(c("-n", "--network"), type="character", default="STRINGv11", 
              help="network being used", metavar ="Network"),
    make_option(c("-a", "--annovarDir"), type="character", default="/home/rhysg/programs/annovar/", 
              help="the directory where annovar is intalled", metavar ="Annovar Directory")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

network_choice <- opt$network
annovar_dir <- opt$annovarDir

# Get sysSVM2 functions

source("scripts/sysSVM2/annotation_functions.R")


# Install ANNOVAR database files
current_dir <- getwd()
setwd(annovar_dir)
if(!file.exists(paste0(annovar_dir,"humandb/hg38_refGene.txt"))){
    message("hg38_refGene database missing")
    system("perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/")
}
if(!file.exists(paste0(annovar_dir,"humandb/hg38_dbnsfp35a.txt"))){
    message("hg38_dbnsfp35a database missing")
    system("perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp35a humandb/")
}
if(!file.exists(paste0(annovar_dir,"humandb/hg38_dbscsnv11.txt"))){
    message("hg38_dbscsnv11 database missing")
    system("perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbscsnv11 humandb/")
}
message("All required ANNOVAR databases available")
setwd(current_dir)

avinput <- "validation_data/CCLE_STRINGv11/ANNOVAR_input/Liver/HEPG2.avinput"

if(!dir.exists("tmp/annovar_tmp")){
    dir.create("tmp/annovar_tmp")
}


ssms_annotated = annotate_ssms(
  avinput=avinput, 
  sample = "HEPG2", 
  annovar_dir, 
  genome_version = "hg38", 
  gene_aliases_entrez = "scripts/sysSVM2/annotation_reference_files/gene_aliases_entrez.tsv", 
  hotspots = "scripts/sysSVM2/annotation_reference_files/tcga_pancancer_hotspots_oncodriveclust.tsv",
  temp_dir = "tmp/annovar_tmp"
)

write_tsv("test_annovar.tsv")



#annovar_cmd = paste0(
#    "perl ", annovar_dir, "/table_annovar.pl ", avinput, 
#    " ", annovar_dir, "/humandb/ -buildver hg38", " -out /home/rhysg/test/", " -remove ", 
#    "-protocol refGene,dbnsfp35a,dbscsnv11 ", 
#    "-operation g,f,f"
#    )

#message("Running ANNOVAR...")
#system(annovar_cmd)
#message("Done")