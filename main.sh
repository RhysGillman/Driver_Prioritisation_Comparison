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
# Dependencies                                             #
############################################################

if (($install_dependencies==1))
then
  #python3 -m pip install -r scripts/PersonaDrive/requirements.txt
  Rscript --vanilla "scripts/install_R_dependencies.r"
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


# If $cell_type = ALL then read sample_info and get all cell_types, then make $cell_type = vector of cell types

# For every lineage in cell_type, run all of the algorithms


    memory_usage () {
        local pid=$1
        local maxmem=0
        while [ -d "/proc/${pid}" ]; do
            local mem=`cat /proc/${pid}/status | grep VmRSS | awk '{print $2}'`
            if [[ ${mem} -gt ${maxmem} ]]; then
                local maxmem=${mem}
            fi
            sleep 1
        done
        echo $maxmem
    }

    ############################################################
    # Run DawnRank                                             #
    ############################################################
    
    if (($run_DawnRank==1))
    then
        # Start time
        start=$(date +%s.%N)
        Rscript --vanilla "scripts/run_DawnRank.R" -n $network_choice -c $cell_type > log/DawnRank_$cell_type.log &
        # Get the process ID (PID) of the  script
        pid=$!
        max_mem=$( memory_usage $pid )
        wait $pid
        end=$(date +%s.%N)
        # Measure time difference
        #runtime=$(echo "$end - $start" | bc -l)
        runtime=$(echo "$end $start" | awk '{print $1 - $2}')
        # Save log information to a file
        echo -e "Runtime_sec\tPeak_VmRSS_KiB" > log/DawnRank_${cell_type}_stats.txt
        echo -e "$runtime\t$max_mem" >> log/DawnRank_${cell_type}_stats.txt
    fi
    
    ############################################################
    # Run PRODIGY                                              #
    ############################################################
    
    if (($run_PRODIGY==1))
    then
        # Start time
        start=$(date +%s.%N)
        Rscript --vanilla "scripts/run_PRODIGY.R" -n $network_choice -c $cell_type > log/PRODIGY_$cell_type.log &
        # Get the process ID (PID) of the  script
        pid=$!
        max_mem=$( memory_usage $pid )
        wait $pid
        end=$(date +%s.%N)
        # Measure time difference
        runtime=$(echo "$end $start" | awk '{print $1 - $2}')
        # Save log information to a file
        echo -e "Runtime_sec\tPeak_VmRSS_KiB" > log/PRODIGY_${cell_type}_stats.txt
        echo -e "$runtime\t$max_mem" >> log/PRODIGY_${cell_type}_stats.txt
    fi
    
    ############################################################
    # Run OncoImpact                                           #
    ############################################################
    
    if (($run_OncoImpact==1))
    then
        Rscript --vanilla "scripts/prepare_OncoImpact_data.R" -w "$SCRIPT_DIR" -n $network_choice -c $cell_type
        # Start time
        start=$(date +%s.%N)

        perl scripts/OncoImpact/oncoIMPACT.pl tmp/tmp_OncoImpact_config.cfg > log/OncoImpact_$cell_type.log & 

        # Get the process ID (PID) of the  script
        pid=$!
        max_mem=$( memory_usage $pid )

        wait $pid

        end=$(date +%s.%N)
        # Measure time difference
        runtime=$(echo "$end $start" | awk '{print $1 - $2}')
        # Save log information to a file
        echo -e "Runtime_sec\tPeak_VmRSS_KiB" > log/OncoImpact_${cell_type}_stats.txt
        echo -e "$runtime\t$max_mem" >> log/OncoImpact_${cell_type}_stats.txt


        Rscript --vanilla "scripts/format_OncoImpact_results.R" -n $network_choice -c $cell_type
        rm -rf "results/CCLE_$network_choice/OncoImpact/$cell_type/ANALYSIS"
        rm -rf "results/CCLE_$network_choice/OncoImpact/$cell_type/COMPLETE_SAMPLES"
        rm -rf "results/CCLE_$network_choice/OncoImpact/$cell_type/INCOMPLETE_SAMPLES"
    fi
    
    ############################################################
    # Run PersonaDrive                                         #
    ############################################################
    
    if (($run_PersonaDrive==1))
    then
        Rscript --vanilla "scripts/prepare_PersonaDrive_data.R" -n $network_choice -c $cell_type
        echo "1. Personalized Bipartite Networks (PBNs)..." > log/PersonaDrive_$cell_type.log

        # Start time
        start=$(date +%s.%N)

        python3 scripts/PersonaDrive/constructing_PBNs.py -o "$SCRIPT_DIR/results/CCLE_$network_choice/PersonaDrive/$cell_type" >> log/PersonaDrive_$cell_type.log &
    
        # Get the process ID (PID) of the  script
        pid=$!
        max_mem1=$( memory_usage $pid )
        wait $pid

        echo "2 - Rank Mutated Genes ..." >> log/PersonaDrive_$cell_type.log
        python3 scripts/PersonaDrive/PersonaDrive.py -o "$SCRIPT_DIR/results/CCLE_$network_choice/PersonaDrive/$cell_type" >> log/PersonaDrive_$cell_type.log &
        # Get the process ID (PID) of the  script
        pid=$!
        max_mem2=$( memory_usage $pid )
        wait $pid

        end=$(date +%s.%N)
        # Measure time difference
        runtime=$(echo "$end $start" | awk '{print $1 - $2}')
        # Save log information to a file
        echo -e "Runtime_sec\tPeak_VmRSS" > log/PersonaDrive_${cell_type}_stats.txt
        echo -e "$runtime\t$max_mem" >> log/PersonaDrive_${cell_type}_stats.txt

        Rscript --vanilla "scripts/format_PersonaDrive_results.R" -n $network_choice -c $cell_type
    fi
    
    ############################################################
    # Run SCS                                                  #
    ############################################################
    
    if (($run_SCS==1))
    then
        
        mkdir -p results/CCLE_$network_choice/SCS/$cell_type
        # Prepare input data
        Rscript --vanilla "scripts/prepare_SCS_data.R" -n $network_choice -c $cell_type
        cd scripts
        matlab -batch "create_matlab_network('../validation_data/CCLE_$network_choice/network_directed.csv')"
        cd SCS

        # Start time
        start=$(date +%s.%N)

        matlab -batch "main_SCS('$network_choice', '$cell_type')" > $SCRIPT_DIR/log/SCS_$cell_type.log &
        pid=$!

        max_mem=$( memory_usage $pid )

        wait $pid
        end=$(date +%s.%N)
        # Measure time difference
        runtime=$(echo "$end $start" | awk '{print $1 - $2}')
        # Save log information to a file
        echo -e "Runtime_sec\tPeak_VmRSS_KiB" > $SCRIPT_DIR/log/SCS_${cell_type}_stats.txt
        echo -e "$runtime\t$max_mem" >> $SCRIPT_DIR/log/SCS_${cell_type}_stats.txt

        cd ..
        matlab -batch "get_SCS_results_names('$network_choice', '$cell_type')"
        cd $SCRIPT_DIR
        Rscript --vanilla "scripts/format_SCS_results.R" -n $network_choice -c $cell_type
        rm scripts/SCS/CNV_internate.txt
        rm scripts/SCS/EXPR_internate.txt
        rm scripts/SCS/SNP_internate.txt
    fi
    
    ############################################################
    # Run PNC                                                  #
    ############################################################
    
    if (($run_PNC==1))
    then
        
        mkdir -p results/CCLE_$network_choice/PNC/$cell_type
        # Prepare input data
        Rscript --vanilla "scripts/prepare_PNC_data.R" -n $network_choice -c $cell_type
        cd scripts
        matlab -batch "create_matlab_network('../validation_data/CCLE_$network_choice/network_directed.csv')"
        cd PNC

        # Start time
        start=$(date +%s.%N)

        matlab -batch "main_PNC('$network_choice', '$cell_type', '$gurobi_path')" > $SCRIPT_DIR/log/PNC_$cell_type.log &
        pid=$!

        max_mem=$( memory_usage $pid )

        wait $pid
        end=$(date +%s.%N)
        # Measure time difference
        runtime=$(echo "$end $start" | awk '{print $1 - $2}')
        # Save log information to a file
        echo -e "Runtime_sec\tPeak_VmRSS_KiB" > $SCRIPT_DIR/log/PNC_${cell_type}_stats.txt
        echo -e "$runtime\t$max_mem" >> $SCRIPT_DIR/log/PNC_${cell_type}_stats.txt

        cd $SCRIPT_DIR
        Rscript --vanilla "scripts/format_PNC_results.R" -n $network_choice -c $cell_type
        
    fi
