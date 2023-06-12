function [predict_driver,omiga,Out_ind]=check_IIP(driver_control_module)
%%Function:calculate the minimum cover set
% Input£º
%      driver_control_module:target controllable subspace of each gene
% Output£º
%      predict_driver:minimum cover set
%      omiga:the weight of the minimum cover set
%      Out_ind:index to evaluate whether the target controllable subspace is empty 
%x = bintprog(f,A,b,Aeq,beq) solves the preceding problem with the additional equality constraint
% min f'*x
% A*x<=b,
% Aeq¡¤x = beq,x binary
%%
Out_ind=1;
Target_subspace=[];predict_driver=[];omiga=[];
[row_driver_control_module,~]=size(driver_control_module);
for i=1:row_driver_control_module
    inter=[ones(length(driver_control_module{i,3}),1)*driver_control_module{i,1} driver_control_module{i,3}];
    Target_subspace=[Target_subspace;inter];
end

if length(Target_subspace)~=0
CC=unique(Target_subspace(:,2));
Out_ind=0;
%renumber the network node to reduce the complexity
add_subspace=Target_subspace;
re_form0=[add_subspace(:,1) add_subspace(:,2)];
re_form=unique(re_form0);
[row_add_subspace,colunm_add_subspace]=size(add_subspace);
re_edge=[];
for i=1:row_add_subspace
    [~,n1]=ismember(add_subspace(i,1),re_form);
    [~,n2]=ismember(add_subspace(i,2),re_form);
    re_edge0=[n1 n2];
    re_edge=[re_edge;re_edge0];
end
%*************************** min f'*x**************************************
N1=max(max(re_edge));
a=zeros(N1(1,1));a1=a(:,1);
f = a1;%
function_node=unique(re_edge(:,1));
f(function_node)=1;
intcon=1:N1;
%**************************A*x<=b***************************************
check_B=zeros(N1);
b0=zeros(N1,1);
[row_re_edge,colunm_re_edge]=size(re_edge);
for i=1:row_re_edge
        check_B(re_edge(i,2),re_edge(i,1))=1;
        b0(re_edge(i,2),1)=1;
end
self_control_node=intersect(re_edge(:,2),re_edge(:,1));
for i=1:length(self_control_node)
check_B(self_control_node(i,1),self_control_node(i,1))=1;
end

A=-check_B;
b=-b0;




%***********************solve the problem**************************************

%**********************MATLAB2012******************************
% OptionsBint=optimset('MaxRLPIter',100000,'NodeSearchStrategy','bn',...
%     'MaxTime',inf,'TolXInteger', 1.e-10);
% x= bintprog(f,A,b,[],[],[],OptionsBint);

%**********************MATLAB2014******************************
lb=zeros(N1,1);
ub=ones(N1,1);
x = intlinprog(f,intcon,A,b,[],[],lb,ub);
%***********************************************************************
x=full(x);
index=find(x(:,1)==1);

predict_driver= re_form(index);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
omiga=length(predict_driver)/length(CC);

end





end

