# CRAN Packages
packages <- c(
    "tidyverse",
    "ggrepel",
    "optparse",
    "BiocManager",
    "snowfall",
    "data.table",
    "maxstat",
    "MASS",
    "igraph",
    "ff",
    "BH",
    "httr",
    "Rcpp",
    "visNetwork",
    "devtools",
    "BiocManager",
    "plyr",
    "mixtools",
    "cowplot",
    "doParallel",
    "reshape2",
    "foreach",
    "TopKLists"
    
)

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages], repos = "https://cloud.r-project.org")
}

# BioConductor Packages

packages <- c(
    "DESeq2",
    "graphite",
    "org.Hs.eg.db",
    "topGO",
    "biomaRt"
)

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  BiocManager::install(packages[!installed_packages])
}

# Github Packages

packages <- c(
    "PCSF"="IOR-Bioinformatics/PCSF"
)

installed_packages <- names(packages) %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  devtools::install_github(packages[!installed_packages],
                         dependencies=TRUE, type="source", force=TRUE)
}
