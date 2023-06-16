function [] = main_PNC( network, cell_type, gurobi_path)

addpath(gurobi_path)

%   $Id: main_PNC.m Created at 2019-05-29 22:22:20 $
%   by Weifeng Guo, Northwestern Polytechtical University, China
%   Copyright (c) 2014-2019 by Key Laboratory of Information Fusion Technology of Ministry of Education in Northwestern Polytechnical University,
%   and key Laboratory of Systems Biology in Shanghai Institutes for Biological Science; 
%   If any problem,pleasse contact shaonianweifeng@126.com for help.

%Remainder: Please install gurobi before running our code
%Remainder: Please install gurobi before running our code
%Remainder: Please install gurobi before running our code

%**************Part 1:Input the information of samples and network information****
%**************sample information**************
%Example:TCGA-example cancer data (BRCA cancer datasets)
expression_tumor_fileName = '../../tmp/tmp_PNC_tumour_expression.txt';
expression_normal_fileName = '../../tmp/tmp_PNC_pseudonormal_expression.txt';
network_fileName = '../../tmp/tmp_network.mat';


%%**************Part 2:Network control methods output the predicted combinational drugs****
% Rhys note, added code to extract out and in degree for driver ranking as
% well as "tumor" for outputting results correctly

[ PNC_driver_result, out_deg, in_deg, tumor ] = PNC(expression_tumor_fileName,expression_normal_fileName,network_fileName);

%%**************Part 3:save the result****

%save PNC_driver_result 

% Rhys note: Below code added by me

writetable(array2table(PNC_driver_result,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), strcat('../../results/CCLE_',network,'/PNC/',cell_type,'/result.csv'), 'WriteRowNames',true)
writetable(array2table(out_deg,"VariableNames",['gene_ID',tumor.textdata(1,2:end)]), strcat('../../results/CCLE_',network,'/PNC/',cell_type,'/out_deg.csv'), 'WriteRowNames',false)
writetable(array2table(in_deg,"VariableNames",['gene_ID',tumor.textdata(1,2:end)]), strcat('../../results/CCLE_',network,'/PNC/',cell_type,'/in_deg.csv'), 'WriteRowNames',false)

end