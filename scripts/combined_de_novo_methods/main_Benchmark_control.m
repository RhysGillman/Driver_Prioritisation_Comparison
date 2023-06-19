function [] = main_PNC( network, cell_type, gurobi_path)

addpath(gurobi_path)

%   $Id: main_Benchmark_control.m Created at 2020-02-05 21:25:22 $
%   by Weifeng Guo, Zhengzhou University, China
%   Copyright (c) 2019-2023 by School of Electrical Engineering, Zhengzhou University, 
%   and key Laboratory of Systems Biology in Shanghai Institutes for Biological Science; 
%   If any problem,pleasse contact shaonianweifeng@126.com for help.

%Remaind: Please install gurobi before running our code
%Remaind: Please install gurobi before running our code
%Remaind: Please install gurobi before running our code

%**************Part 1:Input the information of samples and network information****
%**************sample information**************

expression_tumor_fileName = '../../tmp/tmp_PNC_tumour_expression.txt';
expression_normal_fileName = '../../tmp/tmp_PNC_pseudonormal_expression.txt';
network_fileName = '../../tmp/tmp_network.mat';


[tumor,~,name_tumor]=importdata(expression_tumor_fileName);
gene_list=tumor.textdata(2:end,1);tumor_data=tumor.data;
[normal,~,name_normal]=importdata(expression_normal_fileName);
Sample_name_normal=normal.textdata(1,2:end);normal_data=normal.data;
data=tumor_data;ref_data=normal_data;

%**************the network construction information****
%if Network_index=1,we use CSN; if Network_index=2,we use SSN
%if Network_index=3,we use SPCC; if Network_index=4,we use LIONESS

Network_method_index=1;


%%**************Part 2:PDC outputs the predicted combinational drugs****
%Note that the input variable "ref_data" only is used by SSN,although it is a input variable in our function;

[ MMS,MDS,NCU,NCD, out_deg, in_deg ] = benchmark_control( data,ref_data,gene_list,Network_method_index,network_fileName );
% Rhys note: NCD = DFVS, NCU = NCUA

%%**************Part 3:save the result****

writetable(array2table(MMS,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), strcat('../../results/CCLE_',network,'/combined_de_novo_methods/',cell_type,'/CSN_MMS_result.csv'),'WriteRowNames',true)
writetable(array2table(MDS,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), strcat('../../results/CCLE_',network,'/combined_de_novo_methods/',cell_type,'/CSN_MDS_result.csv'), 'WriteRowNames',true)
writetable(array2table(NCU,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), strcat('../../results/CCLE_',network,'/combined_de_novo_methods/',cell_type,'/CSN_NCUA_result.csv'), 'WriteRowNames',true)
writetable(array2table(NCD,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), strcat('../../results/CCLE_',network,'/combined_de_novo_methods/',cell_type,'/CSN_DFVS_result.csv'), 'WriteRowNames',true)
writetable(array2table(out_deg,"VariableNames",['gene_ID',tumor.textdata(1,2:end)]), strcat('../../results/CCLE_',network,'/combined_de_novo_methods/',cell_type,'/CSN_out_deg.csv'), 'WriteRowNames',false)
writetable(array2table(in_deg,"VariableNames",['gene_ID',tumor.textdata(1,2:end)]), strcat('../../results/CCLE_',network,'/combined_de_novo_methods/',cell_type,'/CSN_in_deg.csv'), 'WriteRowNames',false)
vars = {'MMS','MDS','NCU','NCD','out_deg','in_deg'};
clear(vars{:})

Network_method_index=2;
[ MMS,MDS,NCU,NCD, out_deg, in_deg ] = benchmark_control( data,ref_data,gene_list,Network_method_index,network_fileName );
writetable(array2table(MMS,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), strcat('../../results/CCLE_',network,'/combined_de_novo_methods/',cell_type,'/SSN_MMS_result.csv'),'WriteRowNames',true)
writetable(array2table(MDS,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), strcat('../../results/CCLE_',network,'/combined_de_novo_methods/',cell_type,'/SSN_MDS_result.csv'), 'WriteRowNames',true)
writetable(array2table(NCU,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), strcat('../../results/CCLE_',network,'/combined_de_novo_methods/',cell_type,'/SSN_NCUA_result.csv'), 'WriteRowNames',true)
writetable(array2table(NCD,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), strcat('../../results/CCLE_',network,'/combined_de_novo_methods/',cell_type,'/SSN_DFVS_result.csv'), 'WriteRowNames',true)
writetable(array2table(out_deg,"VariableNames",['gene_ID',tumor.textdata(1,2:end)]), strcat('../../results/CCLE_',network,'/combined_de_novo_methods/',cell_type,'/SSN_out_deg.csv'), 'WriteRowNames',false)
writetable(array2table(in_deg,"VariableNames",['gene_ID',tumor.textdata(1,2:end)]), strcat('../../results/CCLE_',network,'/combined_de_novo_methods/',cell_type,'/SSN_in_deg.csv'), 'WriteRowNames',false)
vars = {'MMS','MDS','NCU','NCD','out_deg','in_deg'};
clear(vars{:})

Network_method_index=3;
[ MMS,MDS,NCU,NCD, out_deg, in_deg ] = benchmark_control( data,ref_data,gene_list,Network_method_index,network_fileName );
writetable(array2table(MMS,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), strcat('../../results/CCLE_',network,'/combined_de_novo_methods/',cell_type,'/SPCC_MMS_result.csv'),'WriteRowNames',true)
writetable(array2table(MDS,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), strcat('../../results/CCLE_',network,'/combined_de_novo_methods/',cell_type,'/SPCC_MDS_result.csv'), 'WriteRowNames',true)
writetable(array2table(NCU,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), strcat('../../results/CCLE_',network,'/combined_de_novo_methods/',cell_type,'/SPCC_NCUA_result.csv'), 'WriteRowNames',true)
writetable(array2table(NCD,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), strcat('../../results/CCLE_',network,'/combined_de_novo_methods/',cell_type,'/SPCC_DFVS_result.csv'), 'WriteRowNames',true)
writetable(array2table(out_deg,"VariableNames",['gene_ID',tumor.textdata(1,2:end)]), strcat('../../results/CCLE_',network,'/combined_de_novo_methods/',cell_type,'/SPCC_out_deg.csv'), 'WriteRowNames',false)
writetable(array2table(in_deg,"VariableNames",['gene_ID',tumor.textdata(1,2:end)]), strcat('../../results/CCLE_',network,'/combined_de_novo_methods/',cell_type,'/SPCC_in_deg.csv'), 'WriteRowNames',false)
vars = {'MMS','MDS','NCU','NCD','out_deg','in_deg'};
clear(vars{:})

Network_method_index=4;
[ MMS,MDS,NCU,NCD, out_deg, in_deg ] = benchmark_control( data,ref_data,gene_list,Network_method_index,network_fileName );
writetable(array2table(MMS,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), strcat('../../results/CCLE_',network,'/combined_de_novo_methods/',cell_type,'/LIONESS_MMS_result.csv'),'WriteRowNames',true)
writetable(array2table(MDS,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), strcat('../../results/CCLE_',network,'/combined_de_novo_methods/',cell_type,'/LIONESS_MDS_result.csv'), 'WriteRowNames',true)
writetable(array2table(NCU,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), strcat('../../results/CCLE_',network,'/combined_de_novo_methods/',cell_type,'/LIONESS_NCUA_result.csv'), 'WriteRowNames',true)
writetable(array2table(NCD,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), strcat('../../results/CCLE_',network,'/combined_de_novo_methods/',cell_type,'/LIONESS_DFVS_result.csv'), 'WriteRowNames',true)
writetable(array2table(out_deg,"VariableNames",['gene_ID',tumor.textdata(1,2:end)]), strcat('../../results/CCLE_',network,'/combined_de_novo_methods/',cell_type,'/LIONESS_out_deg.csv'), 'WriteRowNames',false)
writetable(array2table(in_deg,"VariableNames",['gene_ID',tumor.textdata(1,2:end)]), strcat('../../results/CCLE_',network,'/combined_de_novo_methods/',cell_type,'/LIONESS_in_deg.csv'), 'WriteRowNames',false)
vars = {'MMS','MDS','NCU','NCD','out_deg','in_deg'};
clear(vars{:})


end