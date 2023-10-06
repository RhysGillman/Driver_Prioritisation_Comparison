########################
# General
########################
threads=4

########################
# Dependencies
########################

install_dependencies=0

########################
# Download Data
########################

# Downloads all required data to the /data directory. Set to 1 if the data needs to be downloaded.
# Doesn't include very large files which may be better utilised separately

download_data=0

download_large_data=0


########################
# Gene Interaction Network
########################

# The reference network to be used
# Options --
  # STRINGv11 - STRING v11.0 functional interactions with known directionality
network_choice="STRINGv11"

# Confidence threshold to be used
# Options --
  # Value between 0 and 1. Everything above the threshold will be kept.
network_conf_th=0.4

########################
# Prepare Data
########################

# Performs all re-formatting and filtering required after initial downloading of the data

prepare_data=0

########################
# Cell Types
########################

# CCLE Cell types to analyse
# Options --
  # "ALL" - Loops through all cell types present in the data files
  # Single lineage (eg. "Liver")
  # Multiple lineages (eg. ("Liver" "Brain" "Skin"))

#cell_types=("Liver" "Skin")
cell_types="ALL"
#cell_types="Liver"

########################
# Gold Standards
########################

# Alpha (significance level) to use for identifying global (cell-type centric) outliers in gene dependency and drug sensitivity
# Options --
  # Value between 0 and 1 (Default is 0.1)

local_alpha=0.1


# Alpha (significance level) to use for identifying global (cell-type agnostic) outliers in gene dependency and drug sensitivity
# Options --
  # FALSE - do not use global outliers
  # Value between 0 and 1 (Default is 0.1)

global_alpha=0.1

########################
# Driver Algorithms
########################

# Indicate which algorithms to run
# Options --
  # 0 - Do not run
  # 1 - Run

run_DawnRank=0
run_PRODIGY=0
run_OncoImpact=0
run_PersonaDrive=0
run_SCS=0
run_PNC=0
run_combined_de_novo=0
run_sysSVM2=0
run_PhenoDriverR=1

########################
# badDriver Simulation
########################

# Indicate whether to run
# Options --
  # 0 - Do not run
  # 1 - Run
  
run_badDriver=0

# Use reference set for badDriver simulation
# Options --
  # none - Don't use reference set
  # CGC
  
badDriver_ref_set="none"

# Number of simulations to run

  
badDriver_n_sim=10


########################
# LOF GOF Annotation
########################
# Indicate whether mapping of predicted genomic LOF/GOF annotations is required.

# WARNING!                            This requires ~ 30gb of memory                       WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Options --
  # 0 - Do not run
  # 1 - Run

map_genomic_LOF_GOFs=0


# Indicate whether LOF/GOF annotations need to be run
# Options --
  # 0 - Do not run
  # 1 - Run
  
run_LOF_GOF_annotations=0


########################
# Setup Symbolic Link
########################

# Link Directory
# Some parts of this code struggles to run with complex file paths containing spaces, special characters, etc. If this is the case, 
# the user may optionally create a symbolic link to a new location from which to run the algorithms
# Options --
  # FALSE - do not create link
  # File path to desired link location

use_symbolic_link=true
link_location="/c/Users/jc428796/Driver_Prioritisation_link"
#link_location="/home/rgillman/Driver_Prioritisation_link"
#link_location="/e/Driver_Prioritisation_link"

windows_mode=true






#gurobi_path="C:\gurobi1000\win64\matlab"
gurobi_path="/opt/gurobi1002/linux64/matlab/"  

annovar_path="/home/rhysg/programs/annovar/"
bedtools_path="/home/rhysg/programs/bedtools/"
