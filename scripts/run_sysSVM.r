# run_ANNOVAR.r

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))

# Handling input arguments
option_list = list(
    make_option(c("-n", "--network"), type="character", default="STRINGv11", 
                help="network being used", metavar ="Network"),
    make_option(c("-a", "--annovarDir"), type="character", default="/home/rhysg/programs/annovar/", 
                help="the directory where annovar is installed", metavar ="Annovar Directory"),
    make_option(c("-b", "--bedtoolsDir"), type="character", default="/home/rhysg/programs/bedtools/", 
                help="the directory where bedtools is installed", metavar ="Bedtools Directory")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

network_choice <- opt$network
annovar_dir <- opt$annovarDir
bedtools_dir <- opt$bedtoolsDir

# Get sysSVM2 functions

source("scripts/sysSVM2/annotation_functions_modified.R")
source("scripts/sysSVM2/train_predict_functions.R")

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


ssms_annotated <- annotate_ssms(
  avinput=avinput, 
  sample = "HEPG2", 
  annovar_dir, 
  genome_version = "hg38", 
  gene_aliases_entrez = "scripts/sysSVM2/annotation_reference_files/gene_aliases_entrez.tsv", 
  hotspots = "scripts/sysSVM2/annotation_reference_files/tcga_pancancer_hotspots_oncodriveclust.tsv",
  temp_dir = "tmp/annovar_tmp/"
)



segment_cnv <- read_csv("validation_data/CCLE_STRINGv11/cnv_segment.csv", col_names = c("sample","chromosome","start","end","segment_mean"), col_types = "ccnnn") %>%
    filter(sample=="HEPG2")


if(!dir.exists("tmp/bedtools_tmp")){
    dir.create("tmp/bedtools_tmp")
}

cnvs_annotated <- annotate_cnvs(
    cnv_segments=segment_cnv,            # Table with the following columns: sample; chromosome; start; end; and copy_number or segment_mean
    ploidy = NULL,           # Table with the following columns: sample; ploidy. Leave null if unavailable (assumes diploidy)
    ploidy_threshold = 2,    # Threshold for determining gene amplifications: CN >= ploidy_threshold*ploidy
    gene_coords = "scripts/sysSVM2/annotation_reference_files/gene_coords_hg38.tsv",             # gene_coords_hg19.tsv or gene_coords_hg38.tsv from the sysSVM2 GitHub repository
    bedtools_bin_dir = bedtools_dir, # Directory where the bedtools binary executable is located, if not in $PATH
    temp_dir = "tmp/bedtools_tmp"     # Directory for temporary files to be created
)

molecular_data <- make_totalTable(
  ssms_annotated, cnvs_annotated, 
  canonical_drivers = "scripts/sysSVM2/example_data/canonical_drivers.rds"
)

PANCAN_trained_sysSVM <- readRDS("scripts/sysSVM2/trained_models/PANCAN_trained_sysSVM.rds")

predictions <- predict_sysSVM2(
  trained_sysSVM = PANCAN_trained_sysSVM, 
  molecular_data = molecular_data, 
  systemsLevel_data = "scripts/sysSVM2/example_data/systemsLevel_features_allGenes.tsv"
)
write_tsv(predictions, "sysSVM2_test.tsv")