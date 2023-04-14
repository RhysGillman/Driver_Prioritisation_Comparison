# CRAN Packages
packages <- c(
    "tidyverse",
    "ggrepel",
    "optparse",
    "BiocManager"
)

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages], repos = "https://cloud.r-project.org")
}

# BioConductor Packages



packages <- c(
    "DESeq2"
)

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  BiocManager::install(packages[!installed_packages])
}

