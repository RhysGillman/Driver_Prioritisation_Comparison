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
  echo -e "\n\n---------------------------"
  echo -e "Reading Config File"
  echo -e "---------------------------\n\n"
  source config.sh
fi

############################################################
# Dependencies                                             #
############################################################

if (($install_dependencies==1))
then
  echo -e "\n\n---------------------------"
  echo -e "Installing Python Packages"
  echo -e "---------------------------\n\n"
  python3 -m pip install -r scripts/PersonaDrive/requirements.txt
  echo -e "\n\n---------------------------"
  echo -e "Installing R Packages"
  echo -e "---------------------------\n\n"
  Rscript --vanilla "scripts/install_R_dependencies.r"
fi


############################################################
# Setup Symbolic Link                                      #
############################################################

if [ "$use_symbolic_link" == true ]
then

  rm "$link_location/Driver_Prioritisation_Comparison_Pipeline"
  echo -e "\n\n---------------------------"
  echo -e "Attempting to create link to $link_location"
  echo -e "---------------------------\n\n"
  mkdir -p "$link_location"
  # The below code is needed to add permissions for symoblic links when working in windows Git Bash
  if [ "$windows_mode" == true ]
    then
    export MSYS=winsymlinks:nativestrict
  fi

  ln -sf "$SCRIPT_DIR" "$link_location"

  cd "$link_location/Driver_Prioritisation_Comparison_Pipeline"
  echo -e "\n\n---------------------------"
  echo -e "Link created to $link_location"
  echo -e "---------------------------\n\n"
  SCRIPT_DIR="$link_location/Driver_Prioritisation_Comparison_Pipeline"

fi

############################################################
# Setup Directories                                        #
############################################################

mkdir -p log
mkdir -p validation_data
mkdir -p plots/QC
mkdir -p tmp
mkdir -p results/CCLE_$network_choice

############################################################
# Download Data                                            #
############################################################

if (($download_data==1))
then
    echo -e "\n\n---------------------------"
    echo -e "Downloading Raw Data"
    echo -e "---------------------------\n\n"
    bash scripts/download_data.sh -a
fi

if (($download_large_data==1))
then
    echo -e "\n\n---------------------------"
    echo -e "Downloading Raw Data"
    echo -e "---------------------------\n\n"
    bash scripts/download_data.sh -L
fi

############################################################
# Prepare Data                                             #
############################################################

if (($prepare_data==1))
then
    echo -e "\n\n---------------------------"
    echo -e "Data Processing"
    echo -e "---------------------------\n\n"
    echo -e "\n\n---------------------------"
    echo -e "Step 1: Choosing cells"
    echo -e "---------------------------\n\n"
    Rscript --vanilla "scripts/choose_cells.r" -w $SCRIPT_DIR > log/prepare_data.log
    echo -e "\n\n---------------------------"
    echo -e "Step 2: Identifying Gold-Standard Sensitive Genes"
    echo -e "---------------------------\n\n"
    Rscript --vanilla "scripts/prepare_all_gold_standards.r" -l $local_alpha -g $global_alpha >> log/prepare_data.log
    echo -e "\n\n---------------------------"
    echo -e "Step 3: Filtering and reformatting data"
    echo -e "---------------------------\n\n"
    Rscript --vanilla "scripts/prepare_specific_data.r" -w $SCRIPT_DIR -n $network_choice -c $network_conf_th >> log/prepare_data.log
fi


# Cell Types to analyse

if [ "$cell_types" == "ALL" ] 
then
 
    dos2unix validation_data/CCLE_${network_choice}/lineages.txt
    readarray -t cell_types < validation_data/CCLE_${network_choice}/lineages.txt

fi



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
    
for cell_type in ${cell_types[@]}; do
    echo -e "\n\n---------------------------"
    echo -e "Commencing Driver Prioritisation for cell type: $cell_type"
    echo -e "---------------------------\n\n"
    

    ############################################################
    # Run DawnRank                                             #
    ############################################################
    
    if (($run_DawnRank==1))
    then
        echo -e "\n\n---------------------------"
        echo -e "Running DawnRank for $cell_type"
        echo -e "---------------------------\n\n"
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
        
        echo -e "\n\n---------------------------"
        echo -e "Finished running DawnRank"
        echo -e "Time Taken: $runtime seconds"
        echo -e "Peak Memory Usage: $max_mem KiB"
        echo -e "---------------------------\n\n"
        
        # Save log information to a file
        echo -e "Runtime_sec\tPeak_VmRSS_KiB" > log/DawnRank_${cell_type}_stats.txt
        echo -e "$runtime\t$max_mem" >> log/DawnRank_${cell_type}_stats.txt
        
        
    fi
    
    ############################################################
    # Run PRODIGY                                              #
    ############################################################
    
    if (($run_PRODIGY==1))
    then
        echo -e "\n\n---------------------------"
        echo -e "Running PRODIGY for $cell_type"
        echo -e "---------------------------\n\n"
        
        # Start time
        start=$(date +%s.%N)
        Rscript --vanilla "scripts/run_PRODIGY.R" -n $network_choice -c $cell_type -t $threads > log/PRODIGY_$cell_type.log &
        # Get the process ID (PID) of the  script
        pid=$!
        max_mem=$( memory_usage $pid )
        wait $pid
        end=$(date +%s.%N)
        # Measure time difference
        runtime=$(echo "$end $start" | awk '{print $1 - $2}')
        
        echo -e "\n\n---------------------------"
        echo -e "Finished running PRODIGY"
        echo -e "Time Taken: $runtime seconds"
        echo -e "Peak Memory Usage: $max_mem KiB"
        echo -e "---------------------------\n\n"
        
        # Save log information to a file
        echo -e "Runtime_sec\tPeak_VmRSS_KiB" > log/PRODIGY_${cell_type}_stats.txt
        echo -e "$runtime\t$max_mem" >> log/PRODIGY_${cell_type}_stats.txt
    fi
    
    ############################################################
    # Run OncoImpact                                           #
    ############################################################
    
    if (($run_OncoImpact==1))
    then
        echo -e "\n\n---------------------------"
        echo -e "Running OncoImpact for $cell_type"
        echo -e "---------------------------\n\n"
        
        # Prepare data files for OncoImpact and store in tmp/
        Rscript --vanilla "scripts/prepare_OncoImpact_data.R" -w "$SCRIPT_DIR" -n $network_choice -c $cell_type -t $threads
        # Start time
        start=$(date +%s.%N)
        ####### Running OncoImpact ########
        perl scripts/OncoImpact/oncoIMPACT.pl tmp/tmp_${cell_type}_OncoImpact_config.cfg &> log/OncoImpact_${cell_type}_${network_choice}.log & 

        # Get the process ID (PID) of the  script
        pid=$!
        # Monitor memory usage while running and wait until the process is done
        max_mem=$( memory_usage $pid )
        wait $pid
        # Get finish time
        end=$(date +%s.%N)
        # Measure time difference
        runtime=$(echo "$end $start" | awk '{print $1 - $2}')
        
        echo -e "\n\n---------------------------"
        echo -e "Finished running OncoImpact"
        echo -e "Time Taken: $runtime seconds"
        echo -e "Peak Memory Usage: $max_mem KiB"
        echo -e "---------------------------\n\n"
        
        # Save log information to a file
        echo -e "Runtime_sec\tPeak_VmRSS_KiB" > log/OncoImpact_${cell_type}_${network_choice}_stats.txt
        echo -e "$runtime\t$max_mem" >> log/OncoImpact_${cell_type}_${network_choice}_stats.txt
  
        # Formatting results
        Rscript --vanilla "scripts/format_OncoImpact_results.R" -n $network_choice -c $cell_type
        # Remove temporary files
        rm -rf "results/CCLE_$network_choice/OncoImpact/$cell_type/ANALYSIS"
        rm -rf "results/CCLE_$network_choice/OncoImpact/$cell_type/COMPLETE_SAMPLES"
        rm -rf "results/CCLE_$network_choice/OncoImpact/$cell_type/INCOMPLETE_SAMPLES"
        rm tmp/tmp_${cell_type}_OncoImpact_cnv.txt
        rm tmp/tmp_${cell_type}_OncoImpact_config.cfg
        rm tmp/tmp_${cell_type}_OncoImpact_EXP.txt
        rm tmp/tmp_${cell_type}_OncoImpact_network.txt
        rm tmp/tmp_${cell_type}_OncoImpact_SNP.txt
    fi
    
    ############################################################
    # Run PersonaDrive                                         #
    ############################################################
    
    if (($run_PersonaDrive==1))
    then
        echo -e "\n\n---------------------------"
        echo -e "Running PersonaDrive for $cell_type"
        echo -e "---------------------------\n\n"
        
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
        # Take greatest maxmem value
        if [ $max_mem1 -gt $max_mem2 ]
        then
            max_mem=$max_mem1
        else
            max_mem=$max_mem2
        fi
        
        echo -e "\n\n---------------------------"
        echo -e "Finished running PersonaDrive"
        echo -e "Time Taken: $runtime seconds"
        echo -e "Peak Memory Usage: $max_mem KiB"
        echo -e "---------------------------\n\n"
        
        # Save log information to a file
        echo -e "Runtime_sec\tPeak_VmRSS" > log/PersonaDrive_${cell_type}_stats.txt
        echo -e "$runtime\t$max_mem" >> log/PersonaDrive_${cell_type}_stats.txt

        Rscript --vanilla "scripts/format_PersonaDrive_results.R" -n $network_choice -c $cell_type
        
        # Remove temporary files
        rm tmp/tmp_PersonaDrive_DEGs.csv
        rm tmp/tmp_PersonaDrive_MUT.csv
        rm tmp/tmp_PersonaDrive_network.tsv
    fi
    
    ############################################################
    # Run SCS                                                  #
    ############################################################
    
    if (($run_SCS==1))
    then
        
        echo -e "\n\n---------------------------"
        echo -e "Running SCS for $cell_type"
        echo -e "---------------------------\n\n"
        
        mkdir -p results/CCLE_$network_choice/SCS/$cell_type
        # Prepare input data
        Rscript --vanilla "scripts/prepare_SCS_data.R" -n $network_choice -c $cell_type
        
        if (($network_choice=="own"))
        then
          cp -f data/own_networks/SCS.mat tmp/tmp_network.mat
          cd scripts
          cd SCS
        else
          cd scripts
          matlab -batch "create_matlab_network('../validation_data/CCLE_$network_choice/network_directed.csv')"
          cd SCS
        fi
        

        # Start time
        start=$(date +%s.%N)

        matlab -batch "main_SCS('$network_choice', '$cell_type')" > $SCRIPT_DIR/log/SCS_$cell_type.log &
        pid=$!

        max_mem=$( memory_usage $pid )

        wait $pid
        end=$(date +%s.%N)
        # Measure time difference
        runtime=$(echo "$end $start" | awk '{print $1 - $2}')
        
        echo -e "\n\n---------------------------"
        echo -e "Finished running SCS"
        echo -e "Time Taken: $runtime seconds"
        echo -e "Peak Memory Usage: $max_mem KiB"
        echo -e "---------------------------\n\n"
        
        # Save log information to a file
        echo -e "Runtime_sec\tPeak_VmRSS_KiB" > $SCRIPT_DIR/log/SCS_${cell_type}_stats.txt
        echo -e "$runtime\t$max_mem" >> $SCRIPT_DIR/log/SCS_${cell_type}_stats.txt

        cd ..
        matlab -batch "get_SCS_results_names('$network_choice', '$cell_type')"
        cd $SCRIPT_DIR
        Rscript --vanilla "scripts/format_SCS_results.R" -n $network_choice -c $cell_type
        
        # Remove temporary files
        rm scripts/SCS/CNV_internate.txt
        rm scripts/SCS/EXPR_internate.txt
        rm scripts/SCS/SNP_internate.txt
        rm tmp/tmp_SCS_cnv.txt
        rm tmp/tmp_SCS_EXP.txt
        rm tmp/tmp_SCS_SNP.txt
        rm tmp/tmp_network.mat
        
    fi
    
    ############################################################
    # Run PNC                                                  #
    ############################################################
    
    if (($run_PNC==1))
    then
        echo -e "\n\n---------------------------"
        echo -e "Running PNC for $cell_type"
        echo -e "---------------------------\n\n"
        
        mkdir -p results/CCLE_$network_choice/PNC/$cell_type
        # Prepare input data
        Rscript --vanilla "scripts/prepare_PNC_data.R" -n $network_choice -c $cell_type


        if (($network_choice=="own"))
        then
          cp -f data/own_networks/PNC.mat tmp/tmp_network.mat
          cd scripts
          cd PNC
        else
          cd scripts
          matlab -batch "create_matlab_network('../validation_data/CCLE_$network_choice/network_directed.csv')"
          cd PNC
        fi

        # Start time
        start=$(date +%s.%N)

        matlab -batch "main_PNC('$network_choice', '$cell_type', '$gurobi_path')" > $SCRIPT_DIR/log/PNC_$cell_type.log &
        pid=$!

        max_mem=$( memory_usage $pid )

        wait $pid
        end=$(date +%s.%N)
        # Measure time difference
        runtime=$(echo "$end $start" | awk '{print $1 - $2}')
        
        echo -e "\n\n---------------------------"
        echo -e "Finished running PNC"
        echo -e "Time Taken: $runtime seconds"
        echo -e "Peak Memory Usage: $max_mem KiB"
        echo -e "---------------------------\n\n"
        
        # Save log information to a file
        echo -e "Runtime_sec\tPeak_VmRSS_KiB" > $SCRIPT_DIR/log/PNC_${cell_type}_stats.txt
        echo -e "$runtime\t$max_mem" >> $SCRIPT_DIR/log/PNC_${cell_type}_stats.txt

        cd $SCRIPT_DIR
        Rscript --vanilla "scripts/format_PNC_results.R" -n $network_choice -c $cell_type
        
        
        # Remove temporary files
        rm tmp/tmp_PNC_pseudonormal_expression.txt
        rm tmp/tmp_PNC_tumour_expression.txt
        rm tmp/tmp_network.mat
    fi
    
    
    ############################################################
    # Run De Novo Methods                                      #
    ############################################################
    
    if (($run_combined_de_novo==1))
    then
        echo -e "\n\n---------------------------"
        echo -e "Running Combined De Novo Methods for $cell_type"
        echo -e "---------------------------\n\n"
        
        mkdir -p results/CCLE_$network_choice/combined_de_novo_methods/$cell_type
        # Prepare input data
        Rscript --vanilla "scripts/prepare_combined_de_novo_methods_data.R" -n $network_choice -c $cell_type


        if (($network_choice=="own"))
        then
          cp -f data/own_networks/de_novo_net2.mat tmp/tmp_network.mat
          cd scripts
          cd combined_de_novo_methods
        else
          cd scripts
          matlab -batch "create_matlab_network('../validation_data/CCLE_$network_choice/network_directed.csv')"
          cd combined_de_novo_methods
        fi

        # Start time
        start=$(date +%s.%N)

        matlab -batch "main_Benchmark_control('$network_choice', '$cell_type', '$gurobi_path')" > $SCRIPT_DIR/log/combined_de_novo_methods_$cell_type.log &
        pid=$!

        max_mem=$( memory_usage $pid )

        wait $pid
        end=$(date +%s.%N)
        # Measure time difference
        runtime=$(echo "$end $start" | awk '{print $1 - $2}')
        
        echo -e "\n\n---------------------------"
        echo -e "Finished running combined de novo methods"
        echo -e "Time Taken: $runtime seconds"
        echo -e "Peak Memory Usage: $max_mem KiB"
        echo -e "---------------------------\n\n"
        
        # Save log information to a file
        echo -e "Runtime_sec\tPeak_VmRSS_KiB" > $SCRIPT_DIR/log/combined_de_novo_methods_${cell_type}_stats.txt
        echo -e "$runtime\t$max_mem" >> $SCRIPT_DIR/log/combined_de_novo_methods_${cell_type}_stats.txt

        cd $SCRIPT_DIR
        Rscript --vanilla "scripts/format_combined_de_novo_methods_results.R" -n $network_choice -c $cell_type
        
        # Remove temporary files
        rm tmp/tmp_combined_de_novo_methods_pseudonormal_expression.txt
        rm tmp/tmp_combined_de_novo_methods_tumour_expression.txt
        rm tmp/tmp_network.mat
        
    fi
    
    
    ############################################################
    # Run sysSVM2                                              #
    ############################################################
    
    if (($run_sysSVM2==1))
    then
        echo -e "\n\n---------------------------"
        echo -e "Running sysSVM2 for $cell_type"
        echo -e "---------------------------\n\n"
        
        mkdir -p validation_data/CCLE_$network_choice/ANNOVAR_input/$cell_type
        mkdir -p results/CCLE_$network_choice/sysSVM2/$cell_type
        
        
        # Prepare input data
        
        Rscript --vanilla "scripts/create_ANNOVAR_input_files.R"  -n $network_choice -c $cell_type
        
        # Start time
        start=$(date +%s.%N)
        
        # Run
        
        Rscript --vanilla "scripts/run_sysSVM.r" -n $network_choice -c $cell_type -a $annovar_path -b $bedtools_path  > $SCRIPT_DIR/log/sysSVM2_${network_choice}_${cell_type}.log &
        
        pid=$!

        max_mem=$( memory_usage $pid )

        wait $pid
        end=$(date +%s.%N)
        # Measure time difference
        runtime=$(echo "$end $start" | awk '{print $1 - $2}')
        
        echo -e "\n\n---------------------------"
        echo -e "Finished running sysSVM2"
        echo -e "Time Taken: $runtime seconds"
        echo -e "Peak Memory Usage: $max_mem KiB"
        echo -e "---------------------------\n\n"
        
        # Save log information to a file
        echo -e "Runtime_sec\tPeak_VmRSS_KiB" > $SCRIPT_DIR/log/sysSVM2_${network_choice}_${cell_type}_stats.txt
        echo -e "$runtime\t$max_mem" >> $SCRIPT_DIR/log/sysSVM2_${network_choice}_${cell_type}_stats.txt
        
        Rscript --vanilla "scripts/format_sysSVM2_results.R" -n $network_choice -c $cell_type
        
    
    
    fi
    
    
    ############################################################
    # Run PhenoDriverR                                         #
    ############################################################
    
    if (($run_PhenoDriverR==1))
    then
        echo -e "\n\n---------------------------"
        echo -e "Running PhenoDriverR for $cell_type"
        echo -e "---------------------------\n\n"
        
        
        # Start time
        start=$(date +%s.%N)
        
        # Run
        
        Rscript --vanilla "scripts/run_PhenoDriverR.R" -n $network_choice -c $cell_type -t $threads  > "$SCRIPT_DIR/log/PhenoDriverR_${network_choice}_${cell_type}.log" &
        
        pid=$!

        max_mem=$( memory_usage $pid )

        wait $pid
        end=$(date +%s.%N)
        # Measure time difference
        runtime=$(echo "$end $start" | awk '{print $1 - $2}')
        
        echo -e "\n\n---------------------------"
        echo -e "Finished running PhenoDriverR"
        echo -e "Time Taken: $runtime seconds"
        echo -e "Peak Memory Usage: $max_mem KiB"
        echo -e "---------------------------\n\n"
        
        # Save log information to a file
        echo -e "Runtime_sec\tPeak_VmRSS_KiB" > "$SCRIPT_DIR/log/PhenoDriverR_${network_choice}_${cell_type}_stats.txt"
        echo -e "$runtime\t$max_mem" >> "$SCRIPT_DIR/log/PhenoDriverR_${network_choice}_${cell_type}_stats.txt"
        
    
    fi
    
    
    
    ############################################################
    # Run badDriver Simulation                                 #
    ############################################################
    
    if (($run_badDriver==1))
    then
        echo -e "\n\n---------------------------"
        echo -e "Running badDriver Simulation for $cell_type"
        echo -e "---------------------------\n\n"
        
        
        # Run
        
        Rscript --vanilla "scripts/badDriver.R" -n $network_choice -c $cell_type -t $threads -r $badDriver_ref_set -s $badDriver_n_sim
        
    
    
    fi
    
    
    done
    

############################################################
# LOF GOF Annotations                                      #
############################################################

if (($map_genomic_LOF_GOFs==1))
then
    echo -e "\n\n---------------------------"
    echo -e "Mapping genomic level LOF GOF annotations to CCLE mutations"
    echo -e "WARNING!   This step requires ~ 30gb of memory"
    echo -e "---------------------------\n\n"
    Rscript --vanilla "scripts/map_genomic_LOF_GOFs.R" -t $threads
fi
