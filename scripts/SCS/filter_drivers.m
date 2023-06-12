function [ new_result_driver_gene_module ] = filter_drivers( result_driver_gene_module,sample_network,node0 )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%filter the personal dispensible genes
 new_result_driver_gene_module=[];

for i=1:length(result_driver_gene_module)
    
    
    i
    z_sub_network=sample_network{i,1};
    
    if length(result_driver_gene_module{i,1})~=0
    [ indispensible_genes ] = identify_indispensible( z_sub_network,node0 );
 
    if length(indispensible_genes)~=0
    
    k=1;[row,colunm]=size(result_driver_gene_module{i,1});
    
    
   for j=1:row
       
    [a,b]=ismember(result_driver_gene_module{i,1}{j,1},indispensible_genes);
    if b~=0
      
        new_result_driver_gene_module{i,1}{k,1}=result_driver_gene_module{i,1}{j,1};
        new_result_driver_gene_module{i,1}{k,2}=result_driver_gene_module{i,1}{j,2};
        new_result_driver_gene_module{i,1}{k,3}=result_driver_gene_module{i,1}{j,3};
        k=k+1;
        
    end
   end
    
    end
    
    end
    
end    

 
 

end

