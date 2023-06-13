function [] = main_SCS( network, cell_type)

MAIN_DIR="../..";
%   $Id: main_SCS.m Created at 2017-10-22 16:25:22 $
%   by Weifeng Guo, Northwestern Polytechtical University, China
%   Copyright (c) 2014-2018 by Key Laboratory of Information Fusion Technology of Ministry of Education in Northwestern Polytechnical University,
%   and key Laboratory of Systems Biology in Shanghai Institutes for Biological Science; 
%   If any problem,pleasse contact shaonianweifeng@126.com for help.

%**************Part 1��Input the information of samples and the network****
%*******************default network inforamtion
load(strcat(MAIN_DIR,'/tmp/tmp_SCS_network.mat')) %N=11648
%load('network2_information.mat') %N=6339

%*****Extract the mutation genes and the differentially expressed genes for
%samples
expression_fileName = strcat(MAIN_DIR,'/tmp/tmp_SCS_EXP.txt');
CNV_fileName = strcat(MAIN_DIR,'/tmp/tmp_SCS_cnv.txt');
SNP_fileName = strcat(MAIN_DIR,'/tmp/tmp_SCS_SNP.txt');

%%**************Part 2��SCS outputs the patient-specific driver profiles****

[ result_driver_gene_module ] = SCS( edge0,node0,expression_fileName, CNV_fileName,SNP_fileName )

%%**************Part 3��save the result****
%save SCS_network1_results_GBM

%This code was created to output the results from SCS into a form readable in R

% Loop through the top-level elements of the nested array
for i = 1:length(result_driver_gene_module)
    % Get the current element
    current_element = result_driver_gene_module{i};
    
    % Check if the current element is non-empty
    if ~isempty(current_element)
        current_element = cell2table(current_element);
        % Create a file name for the current element
        file_name = sprintf(strcat(MAIN_DIR,'/results/CCLE_',network,'/SCS/',cell_type,'/result_sample_%d.csv'), i);
        
        % Write the contents of the current element to a .csv file
        writetable(current_element, file_name);
    end
end

end