function [ new_result_driver_gene_module ] = SCS( edge0,node0,expression_fileName, CNV_fileName,SNP_fileName )
%PSCS outputs the patient-specific driver profiles
%   Input:
%         Network information:edge0 and node0
%         Mutation profiles(  CNV_fileName,SNP_fileName)
%         Expression profiles(expression_fileName)       
%  Output:
%         new_result_driver_gene_module:patient-specific driver profiles

[ CNV_name,CNV_data,SNP_name,SNP_data,EXPR_name,EXPR_data ] = sample_information( node0,expression_fileName, CNV_fileName,SNP_fileName);
%*************************************************************

name_gene=node0;
omiga=1;%Fold change==2
%[ Constrained_B,Targets_C,z ] = data_deal( edge,node,CNV_name,CNV_data,SNP_name,SNP_data,EXPR_name,EXPR_data,omiga);
[ Constrained_B,Targets_C,z ] = data_deal( edge0,node0,CNV_name,CNV_data,SNP_name,SNP_data,EXPR_name,EXPR_data,omiga);
U=[1:length(Constrained_B)]';
tic
for i=1:length(Constrained_B)
%  for i=1:1
    i
    B_fda=Constrained_B{find_U(i,U),1};C=Targets_C{find_U(i,U),1};
    [ z_sub_network ] = sample_mutation_network( z,B_fda );
   sample_network{i,1}=z_sub_network;%save the sample specific network  
    
    if length(z_sub_network)>=3
   
    C=intersect(Targets_C{i,1},unique(z_sub_network));
    %Apply CTCA to identify personal driver mutation
    node_new=unique(z_sub_network);
    [~,x1]=ismember(z_sub_network(:,1),node_new);
    [~,x2]=ismember(z_sub_network(:,2),node_new);
    z_new=[x1 x2];
    [~,C_new]=ismember(C,node_new);
    [~,B_fda_new]=ismember(B_fda,node_new);
    B_fda_new(B_fda_new==0)=[];
    
    if length(C_new)*length(B_fda_new)~=0
        [predict_driver_gene_module ] = single_sample_control( z_new,B_fda_new,C_new,node_new,name_gene);     
        result_driver_gene_module{i,1}=predict_driver_gene_module;
    end
    
    
    end
    
    
end
[ new_result_driver_gene_module ] = filter_drivers( result_driver_gene_module,sample_network,node0 );

toc

end

