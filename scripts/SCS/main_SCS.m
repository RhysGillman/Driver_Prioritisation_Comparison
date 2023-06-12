clc
clear
%   $Id: main_SCS.m Created at 2017-10-22 16:25:22 $
%   by Weifeng Guo, Northwestern Polytechtical University, China
%   Copyright (c) 2014-2018 by Key Laboratory of Information Fusion Technology of Ministry of Education in Northwestern Polytechnical University,
%   and key Laboratory of Systems Biology in Shanghai Institutes for Biological Science; 
%   If any problem,pleasse contact shaonianweifeng@126.com for help.

%**************Part 1��Input the information of samples and the network****
%*******************default network inforamtion
load('network1_information.mat') %N=11648
%load('network2_information.mat') %N=6339

%*****Extract the mutation genes and the differentially expressed genes for
%samples
expression_fileName = 'EXPR_trim.txt';
CNV_fileName = 'CNV_trim.txt';
SNP_fileName = 'SNP_trim.txt';

%%**************Part 2��SCS outputs the patient-specific driver profiles****

[ result_driver_gene_module ] = SCS( edge0,node0,expression_fileName, CNV_fileName,SNP_fileName )

%%**************Part 3��save the result****
save SCS_network1_results_GBM

