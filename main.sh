#!/bin/bash
# main.sh
# This is the main script for running the personalised driver prioritisation comparison pipeline

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd "$SCRIPT_DIR"


############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "This is the main script for running the personalised driver prioritisation comparison pipeline"
   echo
   echo "Syntax: main.sh [-h|f|q]"
   echo "options:"
   echo "h     Print help"
   echo "c     Config mode (Default), looks for options set in config.sh"
   echo "f     Fresh-run mode , downloads all data and prepares all required input files"
   echo "q     Quick mode, assumes that data is downloaded and set-up scripts have been run"
   echo "t     Test mode"    
   echo
}
############################################################

############################################################
# Options                                                  #
############################################################

# Set Defaults

test_mode=0
config_mode=1
use_symoblic_link=0

# Process Input Options
while getopts ":hcfql" option; do
    case $option in
        h) # display Help
            Help
            exit;;
        c) #Config mode
            config_mode=1;;
        f) # title
            ;;
        q) # title
            ;;
        t) # Test
            test_mode=1;;
        \?) # Invalid option
            echo "Error: Invalid option"
            exit;;
   esac
done


############################################################
# Read Config File                                         #
############################################################

if (($config_mode==1))
then
  source config.sh
fi

############################################################
# Setup Symbolic Link                                      #
############################################################

if (("$use_symbolic_link" == true));
then

  rm "$link_location/Driver_Prioritisation_Comparison_Pipeline"
  echo "Attempting to create link to $link_location"
  mkdir -p "$link_location"
  # The below code is needed to add permissions for symoblic links when working in windows Git Bash
  if(("$windows_mode" == true));
    then
    export MSYS=winsymlinks:nativestrict
  fi

  ln -sf "$SCRIPT_DIR" "$link_location"

  cd "$link_location/Driver_Prioritisation_Comparison_Pipeline"
  echo "Link created to $link_location"
  SCRIPT_DIR="$link_location/Driver_Prioritisation_Comparison_Pipeline"

fi

############################################################
# Setup Directories                                        #
############################################################

mkdir -p results/CCLE_$network_choice


############################################################
# Download Data                                            #
############################################################

if (($download_data==1))
then
    bash scripts/download_data.sh -a
fi

############################################################
# Prepare Data                                             #
############################################################

if (($prepare_data==1))
then
    Rscript --vanilla "scripts/choose_cells.r" -w $SCRIPT_DIR
    Rscript --vanilla "scripts/prepare_all_gold_standards" -l $local_alpha -g $global_alpha
    Rscript --vanilla "scripts/prepare_specific_data.r" -w $SCRIPT_DIR -n $network_choice -c $network_conf_th
fi


############################################################
# Run DawnRank                                             #
############################################################

if (($run_DawnRank==1))
then
    Rscript --vanilla "scripts/run_DawnRank.R" -n $network_choice -c $cell_type
fi

############################################################
# Run PRODIGY                                              #
############################################################

if (($run_PRODIGY==1))
then
    Rscript --vanilla "scripts/run_PRODIGY.R" -n $network_choice -c $cell_type
fi

############################################################
# Run OncoImpact                                           #
############################################################

if (($run_OncoImpact==1))
then
    Rscript --vanilla "scripts/prepare_OncoImpact_data.R" -w "$SCRIPT_DIR" -n $network_choice -c $cell_type
    perl scripts/OncoImpact/oncoIMPACT.pl tmp/tmp_OncoImpact_config.cfg
fi


