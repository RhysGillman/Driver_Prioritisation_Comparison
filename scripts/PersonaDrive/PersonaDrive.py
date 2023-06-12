#***********************************
# Run the command 'python PersonaDrive.py' in the Terminal.
# Or directly run this file in the IDE.

# Note that, before running this script, you should run the script of constructing_personalized_networks.py to prepare the data for PersonaDrive.
#
#

# This may take several minutes to finish all the samples.
# Exact time mainly depends on your computer and the number of samples in your dataset.
#***********************************
import argparse
from tqdm import tqdm as tqdm
import os
import sys
import pandas as pd
import networkx as nx
from networkx.algorithms import bipartite
import math

#Load pathways
def load_pathways():
    with open("data/kegg_pathways_v1.txt", 'r')as ifile:
        pathways={}
        for line in ifile:
            line = line.rstrip('\n')
            p=line.split('\t\t\t')[0] #pathway ID
            genes=set(line.split('\t\t\t')[1].split(",") [:-1])
            for g in genes:
                if g not in pathways:
                    pathways[g]=set()
                    pathways[g].add(p)
                else:
                    pathways[g].add(p)
    return pathways

def calculate_pps(graphs,outputDir):
    samples=[sample.split('_')[-1].replace('.gml','') for sample in os.listdir(outputDir+'/BP') if '.gml' in sample]

    #load graphs
    D_data={}

    for i in tqdm(os.listdir(outputDir+'/BP')): # for each sample i in S
        if '.gml' in i:
            sample = i.split('_')[-1].replace('.gml','') # retreive samle ID from the graph name

            Bi = nx.read_gml(outputDir+'/BP/'+i) #load Bi graph
            # set of all outliers O_final in the second partition connected to M_i
            D_Bi = [n for n, d in Bi.nodes(data=True) if d['bipartite']==1]
            Di=[]

            # select onmy outliers that corresponds to samlpe_i
            for v_d in D_Bi:
                if v_d.split("_")[0]==sample:
                    Di.append(v_d.split("_")[1])

            D_data[sample]=set(Di)
    pps=pd.DataFrame(0.0,columns=samples,index=samples)
    for i in pps: #for each sample i in S
        D_i=D_data[i]
        for j in pps: #for each sample i in S
            D_j=D_data[j]
            #calculate the similarity score
            if len(D_i)*len(D_j)!=0:
                pps[i][j]=math.pow(len(D_i.intersection(D_j)),2)/(len(D_i)*len(D_j))
    pps.to_csv(outputDir+"/pps_matrix.csv")

# calculate influence scores
def calculate_infl_scores(Bi,Mi,sample_i,similarity,pathways):
    influence_scores={}
    for u in Mi:
        influence_scores[u]=0.0
        for v in Bi.neighbors(u):
            sample_j=v.split("_")[0]
            # check that the mutaated gene u and the outlier gene v exist in the pathways data
            if v.split("_")[1] in pathways and u in pathways:
                # calculate cp score
                ppc_u_v=len(set(pathways[v.split("_")[1]]).intersection(set(pathways[u])))
                # retreive the patient similarity score
                pps_i_j=similarity[sample_i][sample_j]

                influence_scores[u]+= ppc_u_v * pps_i_j

    return influence_scores

#prioritize mutated genes
def PersonaDrive(graphs,pps_matrix,outputDir):
    with open(outputDir+"/PersonaDrive.txt","w") as ofile:
        personalized_drivers={}

        for i in tqdm(os.listdir(outputDir+'/BP')):
            if '.gml' in i:
                #get sample id
                sample = i.split('_')[-1].replace('.gml','')

                #load the graph file
                Bi = nx.read_gml(outputDir+'/BP/'+i)
                #list all mutations
                Mi = [n for n, d in Bi.nodes(data=True) if d['bipartite']==0]
                if len(Mi) == 0:
                    continue
                else:
                    #assign edge weights
                    infl_scores = calculate_infl_scores(Bi, Mi,sample,pps_matrix,pathways)
                    drivers=sorted(infl_scores.items(), key=lambda kv: (-kv[1], kv[0]))
                    ofile.write(sample+"\t")
                    for u in [u_i[0] for u_i in drivers]:
                        ofile.write(u+"\t")
                    ofile.write("\n")


if __name__ == "__main__":

    description="Ranking Mutated Genes in PBNs"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-o', "--outputDir", type=str, required=True, default='', help="Path to output Directory")
    
    args = parser.parse_args()
    outputDir = args.outputDir

    graphs=outputDir+'/BP'

    #~~~~~~~~~~~~~Step 1：~~~~~~~~~~~~~~~~~~
    # load pathways data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pathways = load_pathways()
    print('pathways loaded...')


    #~~~~~~~~~~~~~Step 1：~~~~~~~~~~~~~~~~~~
    # construct pps matrix
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if "pps_matrix.csv" not in os.listdir(outputDir):
        calculate_pps(graphs,outputDir)
    pps_file=outputDir+"/pps_matrix.csv"
    pps_matrix=pd.read_csv(pps_file,index_col=0)
    print('PPS matrix loaded...')

    #~~~~~~~~~~~~~Step 1：~~~~~~~~~~~~~~~~~~
    # run PersonaDrive
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PersonaDrive(graphs,pps_matrix,outputDir)
