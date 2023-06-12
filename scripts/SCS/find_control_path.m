function [ path_source ] = find_control_path( z_layer_source )
%Function:Adjoining two elements in N step and  N + 1 step
%Input:
%     z_layer_source:each row represents that the left nodes control which right nodes in each iterated bipartite 
%Output:
%     path_source:the first colunm is the control path£¬the second colunm is
%     the controlled targets
%simple example£º
% z_layer_source=[];
% z_layer_source{1,1}={2,3,3;2,4,[3;4];};
% z_layer_source{2,1}={3,4,4;3,5,[4;5];4,5,5;};
% z_layer_source{3,1}={5,6,6;4,5,5;4,6,[5;6];5,7,[6;7];};
%%
%Firstly process the cell, the first and second column represents the end and the first element in the control path of first N step
%the third column represents the first N step control path
re_z_layer=[];
[row_re_z_layer,~]=size(z_layer_source);   
for i=1:row_re_z_layer
    
    [a,b]=size(z_layer_source{i,1});
    X=[];
    for j=1:a
    X{j,1}=z_layer_source{i,1}{j,1};
    X{j,2}=z_layer_source{i,1}{j,2};
    X{j,3}=[z_layer_source{i,1}{j,1};cell2mat(z_layer_source{i,1}(j,3))];
    end
    re_z_layer{i,1}=X;
end

%Process from the first row element , move forward on the basis of the wheel, until finish all of the elements
inter_path=re_z_layer{1,1};
[row_inter_path,~]=size(inter_path);
   path=[];
  for i=1:row_inter_path
      
      path{i,1}=inter_path{i,3};
      path{i,2}=inter_path{i,3}(end,1);
  end

  
[row_re_z_layer,~]=size(re_z_layer);   
for ii=1:row_re_z_layer-1
    
    
    PP=cell2mat(path(:,2));
    ZZ=cell2mat(re_z_layer{ii+1,1}(:,1));
    
     [~,ind]=ismember(ZZ,PP);
    inter_path=[];
    for jj=1:length(ind)
        xx=cell2mat(path(ind(jj,1),1));
        yy=cell2mat(re_z_layer{ii+1,1}(jj,3));
        inter_path{jj,1}=[xx(1:end-1,1);yy];
    end
    
   %Eliminate duplicate elements in inter_path  

    mark=[];
    for kk=1:length(inter_path)
        for kj=kk+1:length(inter_path)
          if isequal(inter_path{kk,1},inter_path{kj,1})
            mark(kk,1)=1;
          end
        end
    end
   inter_path(find(mark==1))=[];
   

   path=[];
  for i=1:length(inter_path)
      path{i,1}=inter_path{i,1};
      path{i,2}=inter_path{i,1}(end,1);
  end
    
  
end


path_source=path(:,1);
            
  
end
            
    


   
 


