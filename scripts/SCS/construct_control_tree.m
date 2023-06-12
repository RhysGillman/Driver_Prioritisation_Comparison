function [ z_layer_source,Target_of_source ] = construct_control_tree( TK,source,source_address )
 %Function:construct control tree
 %Input:   TK:Maximum matching in the linking and dynamic graph
 %         source:candidate mutation gene
 %         source_address:The address number of muattion gene in the linking and dynamic graph 
 %Output:
 %         z_layer_source:control tree
 %         Target_of_source:controlled targets by source
% %  %simple example:
%  source=2;source_address=2;TK={[2,3;3,4;4,5];[1,2;2,3;3,4];[1,3];};
%%
%*************construct control tree at a given gene*********************************
        q1=1;q2=1;%initial parameter
        num_iter=0;
        layer=[];
        while q1 
       
        %Obtain the first neighborhood of the source gene£¬the first colunm is the source gene£¬the second colunm is the controlled target genes
        %the third colunm is the iternate genes
        %save it in the cell
        all_layer_connect=[]; layer_connect=[];
        for ii=1:length(source)
            
            iter_source=source(ii,1);
            source_ind=find(TK{source_address-num_iter,1}(:,1)==iter_source);
            target=TK{source_address-num_iter,1}(source_ind,2);
            search=target;intermidate=target;
        while q2
            
         if length(search)==0
              q2=0;
              break;
         end    
       
         
         ind_intermidate=find(TK{source_address-num_iter,1}(:,1)==search);
   
         if length(ind_intermidate)==0
             q2=0;
         end    
         
         if length(ind_intermidate)~=0
            target=TK{source_address-num_iter,1}(ind_intermidate,2);
            intermidate=[intermidate;target];
            search=target;
         end
         
         if length(unique(intermidate))<length(intermidate)
             q2=0;
         end
          
        end
        q2=1;
        
        %intermidate is the first layer in the control tree
        
        for i=1:length(intermidate)
            layer_connect{i,1}=iter_source;
            layer_connect{i,2}=intermidate(i,1);
            layer_connect{i,3}=intermidate(1:i,1);
        end
        
        all_layer_connect=[all_layer_connect;layer_connect];
        layer_connect=[];
        
        end
        
        %save the controlled targeted genes
        num_iter=num_iter+1;
        layer{num_iter,1}=all_layer_connect;
        
        
        source=intermidate;
      
        if (source_address-num_iter)==0
            q1=0;
            break;
        end
        
        end
        
         %delete the null in the layer
         mark=[];
         for io=1:length(layer)
             if length(layer{io,1})==0
             mark(io,1)=1;
             
             end
         end
         layer(find(mark==1))=[];
         
          z_layer_source=layer;%control tree
         Target_of_source=unique(cell2mat(layer{end,1}(:,2)));
      


end

