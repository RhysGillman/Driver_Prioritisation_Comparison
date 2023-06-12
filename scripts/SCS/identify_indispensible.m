function [ indispensible_genes ] = identify_indispensible( z_sub_network,node0 )
%function:identify indispensible genes
%   Input:
%        z:network structure
%        node:gene name     

node_new=unique(z_sub_network);
[~,x1]=ismember(z_sub_network(:,1),node_new);
[~,x2]=ismember(z_sub_network(:,2),node_new);
z_new=[x1 x2];
[manner0,Nd0]=full_control(z_new);%full control of network
N=max(max(z_new));
all_driver_number=[];

parfor i=1:N
    
    z0=z_new;
    index0=[find(z0(:,1)==i);find(z0(:,2)==i)];
    if length(index0)~=0
    index=unique(index0);
    z0(index,:)=[];
    if length(z0)~=0
    [manner,Nd]=full_control(z0);%full control of network
    end
    
    if length(z0)==0
      manner=[];
      Nd=0;
    end
    
    all_manner{i,1}=manner;
    all_Nd{i,1}=Nd;
    all_driver_number(i,1)=length(Nd)-length(Nd0);
    end
    
    if length(index0)==0
        
       all_manner{i,1}=manner0;
       all_Nd{i,1}=Nd0;
       all_driver_number(i,1)=0;
        
    end
   
end


indispensible_genes=node0(node_new(find(all_driver_number>1)));


end

