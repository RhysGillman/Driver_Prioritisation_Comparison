function [ driver_control_module] = re_check_subspace( TK,B_fda)
%Function:Construct control tree for identifying the target control subspace
%Input:   TK:The maximum matching in the linking and dynamic graph
%         B_fda:mutations
%Output:driver_control_module,the first colunm is the ID number of mutation£¬
%The second colunm is the control path from the driver to the targets£¬the
%third colunm is the controlled targets
%%
%index_control£ºthe element in the row is the iterated maching£¬the colunm is the mutation
if length(TK)==1
    driver_control_module=[];
    Control_B_fda=intersect(B_fda,TK{1,1}(:,1));
    for i=1:length(Control_B_fda)
     driver_control_module{i,1}=Control_B_fda(i,1);
      driver_control_module{i,2}=TK{1,1};
      driver_control_module{i,3}=TK{1,1}(:,2);
    end
    
    
end


if length(TK)>1

index_control=zeros(length(TK),length(B_fda));
for k=1:length(B_fda)
    
    candidate_mutation=B_fda(k,1);
    for j=1:length(TK)
        if length(intersect(candidate_mutation,TK{j,1}(:,1)))~=0
            index_control(j,k)=1;
        end
    end
end   
find_control_mutation=sum(index_control);
index_control_mutation=find(find_control_mutation~=0)';
%Find the control path from the mutation
%%
%For TK in every step of the maximum matching, reform a new bipartite graph
driver_control_module=[];
Target_node=[];
for i=1:length(index_control_mutation)
    
    ind=index_control_mutation(i,1);
    [row_ind_TK,~]=find(index_control(:,ind)==1);
    Module_source=[];Controlled=[];
    for j=1:length(row_ind_TK)
         
      source=B_fda(ind,1);
      source_address=row_ind_TK(j,1);
      %Generated the control trees z from the source nodes (in the position of source_address) to , and find controlled targeted node
      [z_layer_source,Target_of_source ] = construct_control_tree( TK,source,source_address );
      %Then find the path from the source to the destination node
      [path_source] = find_control_path( z_layer_source ); 
      for iu=1:length(path_source)
          xv=path_source{iu,1};
      if length(xv)>1
         Module=[xv(1:end-1) xv(2:end)];
      end
      if length(xv)==1
         Module=[xv xv];
      end
      Module_source=[Module_source;Module];
      end
      
      Controlled=[Controlled;Target_of_source];  
    end
    %forming the driver control module
    driver_control_module{i,1}=B_fda(index_control_mutation(i,1),1);
    driver_control_module{i,2}=unique(Module_source,'rows');
    driver_control_module{i,3}=unique(Controlled);
end



end

end




