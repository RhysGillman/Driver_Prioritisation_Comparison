function [ z_sub_network ] = sample_mutation_network( z,B_fda )
%Function:
%         forming the mutation network for each sample
%Input: 
%         z:original networks
%         mutation：样本的变异基因编号
%Output：
%         z_sub_network：样本的变异特异性网络 
mutation=B_fda;
N=max(max(z));
%邻接矩阵
A=zeros(N);
for i=1:N
    A(z(i,2),z(i,1))=1;
end
%*************构建样本特异性mutation 网络**************************
[ q_new ] = RWR( A,mutation );%重启随机游走算法
%随机置乱网络的结构，重新计算mutation gene在网络中游走的概率
%并统计显著性，选取有显著差异的节点，构成sample-specific的mutation network
number=1;ind=1;
num_max=50;
while ind
dir_srand=dir_generate_srand(A);
[ q_rand ] = RWR( dir_srand,mutation );%重启随机游走算法
Q_rand(:,number)=q_rand;
number=number+1;
if number>num_max
    ind=0;
end
end
%统计100次产生节点相似性的平均值和方差
%并计算p_value
[row,colunm]=size(Q_rand);
for i=1:row
    p(i,1)=(max(Q_rand(i,:))<q_new(i,1));
end
p=double(p);
q_new0=q_new;
q_new0(q_new0>0.0001)=1;
p_value=q_new0.*p;
down_stream=find(p_value>0);
subnetwork_nodes=unique([mutation;down_stream]');
%提取出子网络的结构
[x1,y1]=ismember(z(:,1),subnetwork_nodes);
[x2,y2]=ismember(z(:,2),subnetwork_nodes);
ind=double(x1).*double(x2);
z_sub_network=z(find(ind~=0),:);


end

