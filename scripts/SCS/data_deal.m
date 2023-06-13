function  [ Constrained_B,Targets_C,z ] = data_deal( edge,node,CNV_name,CNV_data,SNP_name,SNP_data,EXPR_name,EXPR_data,omiga)
%%
%Function:deal with the original data
%Input:
%      edge:Interaction between genes
%      gene_ID:The ID of all genes
%      Mutation profiles:CNV_data,SNP_data
%      Expression profiles:EXPR_data
%Output:
%       z:adjacy list
%       Constrained_B:The address number of individual mutations
%       Targets:The address number of differentially expressed genes for each sample
%%
[~,z1]=ismember(edge(:,1),node);
[~,z2]=ismember(edge(:,2),node);
z=[z1 z2];%the element is the address number
z=unique(z,'rows');
name_gene=node;
%find the mutation genes£¨CNV and SNP£©
[gene_row,sample_colunm]=size(EXPR_data);
omiga=1;%Fold change==2
for i=1:sample_colunm
    %i  %the i-th sample
    %obtain the differentially expressed genes as teh target_genes
    expr_i=abs(EXPR_data(:,i));
    [index,~]=find(expr_i>omiga);
    target_genes=EXPR_name(index);
    Targets_genes_sampels{i,1}=target_genes;
       
    %Identify the individual mutation_gene from CNV and SNP
    cnv_i=abs(CNV_data(:,i));
    [index1,~]=find(cnv_i~=0);
    snp_i=abs(SNP_data(:,i));
    [index2,~]=find(snp_i~=0);
    CNV_genes=[];
    SNP_genes=[];
    
    if length(index1)~=0
    CNV_genes=CNV_name(index1);
    end
    
    if length(index2)~=0
    SNP_genes=SNP_name(index2);
    end
    
   
    Mutation_gene=[CNV_genes;SNP_genes];
    Mutation_genes_sampels{i,1}=Mutation_gene;
end

%calculate the number of target_genes and mutation_gene
for i=1:length(Mutation_genes_sampels)
    Number_targets(i,1)=length( Targets_genes_sampels{i,1});
    Number_mutation(i,1)=length( Mutation_genes_sampels{i,1});
end
Num=[ Number_targets Number_mutation];


Number=Number_targets.* Number_mutation;
index_delete_sample=find(Number==0);

%assign the gene name as the address number
for i=1:length(Targets_genes_sampels)
    
    [~,ind]=ismember(Targets_genes_sampels{i,1},name_gene);
    ind(ind==0)=[];
    Targets_C{i,1}=ind;
end


for i=1:length(Mutation_genes_sampels)
    
    [~,ind]=ismember(Mutation_genes_sampels{i,1},name_gene);
    ind(ind==0)=[];
    Constrained_B{i,1}=ind;
end


end

