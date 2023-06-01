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
config_mode=0

# Process Input Options
while getopts ":hcfq" option; do
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
# Setup Directories                                        #
############################################################

mkdir -p results/CCLE_$network_choice


############################################################
# Run                                                       #
############################################################


if (($test_mode==1))
then
    Rscript --vanilla "scripts/test.r"
fi
