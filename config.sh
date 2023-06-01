# Cell type to analyse
# Options --
  # ALL - Loops through all cell types present in the data files
  
cell_type="ALL"


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



  
