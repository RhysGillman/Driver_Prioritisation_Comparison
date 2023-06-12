function  [ TK,predict_driver,omiga,driver_module ] = New_Markov_chain_driver( z,B_fda,C,NM,TK0,interesting_ind0)
%%
%function:
%        Based on the previous Markov chain,obtain the new Markov chain
%        Note：choose v uniformly at random from Mt and choose w
%                 uniformly from the matched links to be replaced ,M(t+1)=M(t)-{v}+{w}          
%****************************************************************************************
%Input:
%      z----the connection list
%      B----ALL nodes in the network
%      B_fda--constrained control nodes(individual mutations)
%      C------targeted genes
%      NM=1，the default parameter in the Markov chain sampling
%      TK0----previous Markov chain
%****************************************************************************************
%Output：
%      final_predict_driver1------new driver set
%      omiga1----the weight of Markov chain
%      TK1-------the new Markov chain
%****************************************************************************************
%%
if interesting_ind0==1
    TK0=TK0;
end

if interesting_ind0==0
    New_TK0=[];
    for i=1:length(TK0)
        cand=TK0{i,1};
        New_TK0{length(TK0)-i+1,1}=[cand(:,2) cand(:,1)];
    end
    TK0=New_TK0;
end

Out_ind=1;TK=[];
while Out_ind

[row_TK0,colunm_TK0]=size(TK0);
P1 = randperm(row_TK0)';
P1(find(P1==1))=[];

index=1;i=1;
while index
    
  if length(P1)==0
       break;
    index=0;
    Out_ind=0;
  end
      
if i>length(P1)
    break;
    index=0;
    Out_ind=0;
end

    k=P1(i,1); 
   [TK,eva_index] = check_re_Greedy( z,C,NM,TK0,k);
   if length(TK)==0
       TK=TK0;
   end
  

if eva_index==0
    i=i+1;
end

    

if eva_index==1
    index=0;
end
    
if length(P1)==0
   TK=TK0;
   index=0;
end

if i>length(P1)
   TK=TK0;
   index=0;
end

end


if interesting_ind0==1
    
    TK=TK;
    
end

if interesting_ind0==0
    New_TK=[];
    for i=1:length(TK)
        cand=TK{i,1};
        New_TK{length(TK)-i+1,1}=[cand(:,2) cand(:,1)];
    end
    TK=New_TK;
end


%%
%****************************************************************************************
%Based on TK，identify the control path from driver mutation
%driver_control_module
if length(TK)~=0
[ driver_control_module] = re_check_subspace( TK,B_fda);
%****************************************************************************************
%According to the driver_control_module,identify teh minimum cover set as predict_driver and the contrallable target nodes and control module

[predict_driver,omiga,Out_ind]=check_IIP(driver_control_module);

[ driver_module ] = form_driver_profile( driver_control_module,predict_driver );%生成target control module
    
end

if length(TK)==0
TK=TK0;
[ driver_control_module] = re_check_subspace( TK,B_fda);
%****************************************************************************************
%According to the driver_control_module,identify teh minimum cover set as predict_driver and the contrallable target nodes and control module

[predict_driver,omiga,Out_ind]=check_IIP(driver_control_module);

[ driver_module ] = form_driver_profile( driver_control_module,predict_driver );%生成target control module
 

end






end


