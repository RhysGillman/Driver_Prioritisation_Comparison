clc
clear
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

expression_tumor_fileName = '../../Driver_Gene_Targeting/CCLE_only_validation_data/de_novo_networks_input/tumour_expression.txt';
expression_normal_fileName = '../../Driver_Gene_Targeting/CCLE_only_validation_data/de_novo_networks_input/pseudonormal_expression.txt';


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

[ MMS,MDS,NCU,NCD, out_deg, in_deg ] = benchmark_control( data,ref_data,gene_list,Network_method_index );
% Rhys note: NCD = DFVS, NCU = NCUA

%%**************Part 3:save the result****

writetable(array2table(MMS,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), '../../Driver_Gene_Targeting/CCLE_only_results/de_novo_networks/CSN_MMS_results.csv', 'WriteRowNames',true)
writetable(array2table(MDS,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), '../../Driver_Gene_Targeting/CCLE_only_results/de_novo_networks/CSN_MDS_results.csv', 'WriteRowNames',true)
writetable(array2table(NCU,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), '../../Driver_Gene_Targeting/CCLE_only_results/de_novo_networks/CSN_NCUA_results.csv', 'WriteRowNames',true)
writetable(array2table(NCD,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), '../../Driver_Gene_Targeting/CCLE_only_results/de_novo_networks/CSN_DFVS_results.csv', 'WriteRowNames',true)
writetable(array2table(out_deg,"VariableNames",['gene_ID',tumor.textdata(1,2:end)]), '../../Driver_Gene_Targeting/CCLE_only_results/de_novo_networks/CSN_out_deg.csv', 'WriteRowNames',false)
writetable(array2table(in_deg,"VariableNames",['gene_ID',tumor.textdata(1,2:end)]), '../../Driver_Gene_Targeting/CCLE_only_results/de_novo_networks/CSN_in_deg.csv', 'WriteRowNames',false)
vars = {'MMS','MDS','NCU','NCD','out_deg','in_deg'};
clear(vars{:})

Network_method_index=2;
[ MMS,MDS,NCU,NCD, out_deg, in_deg ] = benchmark_control( data,ref_data,gene_list,Network_method_index );
writetable(array2table(MMS,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), '../../Driver_Gene_Targeting/CCLE_only_results/de_novo_networks/SSN_MMS_results.csv', 'WriteRowNames',true)
writetable(array2table(MDS,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), '../../Driver_Gene_Targeting/CCLE_only_results/de_novo_networks/SSN_MDS_results.csv', 'WriteRowNames',true)
writetable(array2table(NCU,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), '../../Driver_Gene_Targeting/CCLE_only_results/de_novo_networks/SSN_NCUA_results.csv', 'WriteRowNames',true)
writetable(array2table(NCD,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), '../../Driver_Gene_Targeting/CCLE_only_results/de_novo_networks/SSN_DFVS_results.csv', 'WriteRowNames',true)
writetable(array2table(out_deg,"VariableNames",['gene_ID',tumor.textdata(1,2:end)]), '../../Driver_Gene_Targeting/CCLE_only_results/de_novo_networks/SSN_out_deg.csv', 'WriteRowNames',false)
writetable(array2table(in_deg,"VariableNames",['gene_ID',tumor.textdata(1,2:end)]), '../../Driver_Gene_Targeting/CCLE_only_results/de_novo_networks/SSN_in_deg.csv', 'WriteRowNames',false)
vars = {'MMS','MDS','NCU','NCD','out_deg','in_deg'};
clear(vars{:})

Network_method_index=3;
[ MMS,MDS,NCU,NCD, out_deg, in_deg ] = benchmark_control( data,ref_data,gene_list,Network_method_index );
writetable(array2table(MMS,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), '../../Driver_Gene_Targeting/CCLE_only_results/de_novo_networks/SPCC_MMS_results.csv', 'WriteRowNames',true)
writetable(array2table(MDS,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), '../../Driver_Gene_Targeting/CCLE_only_results/de_novo_networks/SPCC_MDS_results.csv', 'WriteRowNames',true)
writetable(array2table(NCU,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), '../../Driver_Gene_Targeting/CCLE_only_results/de_novo_networks/SPCC_NCUA_results.csv', 'WriteRowNames',true)
writetable(array2table(NCD,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), '../../Driver_Gene_Targeting/CCLE_only_results/de_novo_networks/SPCC_DFVS_results.csv', 'WriteRowNames',true)
writetable(array2table(out_deg,"VariableNames",['gene_ID',tumor.textdata(1,2:end)]), '../../Driver_Gene_Targeting/CCLE_only_results/de_novo_networks/SPCC_out_deg.csv', 'WriteRowNames',false)
writetable(array2table(in_deg,"VariableNames",['gene_ID',tumor.textdata(1,2:end)]), '../../Driver_Gene_Targeting/CCLE_only_results/de_novo_networks/SPCC_in_deg.csv', 'WriteRowNames',false)
vars = {'MMS','MDS','NCU','NCD','out_deg','in_deg'};
clear(vars{:})

Network_method_index=4;
[ MMS,MDS,NCU,NCD, out_deg, in_deg ] = benchmark_control( data,ref_data,gene_list,Network_method_index );
writetable(array2table(MMS,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), '../../Driver_Gene_Targeting/CCLE_only_results/de_novo_networks/LIONESS_MMS_results.csv', 'WriteRowNames',true)
writetable(array2table(MDS,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), '../../Driver_Gene_Targeting/CCLE_only_results/de_novo_networks/LIONESS_MDS_results.csv', 'WriteRowNames',true)
writetable(array2table(NCU,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), '../../Driver_Gene_Targeting/CCLE_only_results/de_novo_networks/LIONESS_NCUA_results.csv', 'WriteRowNames',true)
writetable(array2table(NCD,"RowNames",tumor.textdata(2:end,1),"VariableNames",tumor.textdata(1,2:end)), '../../Driver_Gene_Targeting/CCLE_only_results/de_novo_networks/LIONESS_DFVS_results.csv', 'WriteRowNames',true)
writetable(array2table(out_deg,"VariableNames",['gene_ID',tumor.textdata(1,2:end)]), '../../Driver_Gene_Targeting/CCLE_only_results/de_novo_networks/LIONESS_out_deg.csv', 'WriteRowNames',false)
writetable(array2table(in_deg,"VariableNames",['gene_ID',tumor.textdata(1,2:end)]), '../../Driver_Gene_Targeting/CCLE_only_results/de_novo_networks/LIONESS_in_deg.csv', 'WriteRowNames',false)
vars = {'MMS','MDS','NCU','NCD','out_deg','in_deg'};
clear(vars{:})