function [ TK,predict_driver,omiga,driver_module,interesting_ind ] = Markov_chain_driver( z,B_fda,C,NM )
%obtain the Markov chain to obtain the driver set
%Input:
%      z----the connection list
%      B----ALL nodes in the network
%      B_fda--constrained control nodes(individual mutations)
%      C------targeted genes
%      NM=1，the default parameter in the Markov chain sampling
%Output：
%      TK-------linking and dynamic graph的最大匹配
%      path----driver nodes的控制路径
%      predict_driver--预测的driver set
%      omiga----Markov chain的权重
%
interesting_ind=1;
N=max(max(z));
[TK ] = check_Greedy( z,C,NM );
predict_driver=[];
omiga=[];driver_module=[];

if length(TK)~=0
%****************************************************************************************
%****************************************************************************************
%Based on TK，identify the control path from driver mutation
%driver_control_module
[ driver_control_module] = re_check_subspace( TK,B_fda);
if length(driver_control_module)==0
    interesting_ind=0;
    z_inter(:,1)=z(:,2);z_inter(:,2)=z(:,1);
    [TK_inter] = check_Greedy( z_inter,B_fda,NM );
    [row_TK_inter,~]=size(TK_inter);
    %Deal with TK_inter to form TK
    for iy=1:row_TK_inter
        candidate=TK_inter{iy,1};
        New_TK{row_TK_inter-iy+1,1}=[candidate(:,2) candidate(:,1)];
    end
    TK=New_TK;
end

[ driver_control_module] = re_check_subspace( TK,B_fda);
%****************************************************************************************
%According to the driver_control_module,identify teh minimum cover set as predict_driver and the contrallable target nodes and control module
[predict_driver,omiga]=check_IIP(driver_control_module);
[ driver_module ] = form_driver_profile( driver_control_module,predict_driver );%生成target control module

end

end

