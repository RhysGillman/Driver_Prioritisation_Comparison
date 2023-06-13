########################
# Dependencies
########################

install_dependencies=0

########################
# Download Data
########################

# Downloads all required data to the /data directory. Set to 1 if the data needs to be downloaded.

download_data=0

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
# Cell Type
########################

# CCLE Cell type to analyse
# Options --
  # "ALL" - Loops through all cell types present in the data files
  # Single lineage (eg. "Liver", "Breast", "Skin", etc.)

cell_type="Liver"


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
run_SCS=1



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


windows_mode=true



  
