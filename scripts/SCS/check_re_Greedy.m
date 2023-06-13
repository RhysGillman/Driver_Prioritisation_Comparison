function [TK,eva_index] = check_re_Greedy( z,C,NM,TK0,k)
%%
%Function:randomly find new maximum matching in iterated bipartite graph
%Input:
%     TKO:maximum matching in linking and dynamic graph
%     z:the adjacency list 
%     C:targeted nodes
%     NM:The maximum iterated number
%     k:the address of generated new maximum matching
%Output:
%     TK:new maximum matching in linking and dynamic graph
%     eva_index:evaluate whether we find new maximum matching
% A simple example;
% z=[1 2;1 3;3 4;3 5;3 6;5 6;6 7];R0=[2 4 6 7];
%%
TK(1:k-1,1)=TK0(1:k-1,1);
now_matching=TK0{k,1};
Pre_R0=TK0{k-1,1};
% randomly choose the address of the R0 to produce a new match
R0=Pre_R0(:,1);
% The first parameter: m
 m=now_matching;

 N=max(max(z));
 m(:,2)= m(:,2)+N;
 
 m=m(:,1:2);
% The second parameter:Z
Z=z;
 num=1;E=[];
for i=1:length(z)
    [x1,y1]=ismember(z(i,2),R0);
    if x1~=0
        E(num,:)=z(i,:);num=num+1;
    end    
end
Z=[];
Z(:,1)=E(:,1);Z(:,2)=E(:,2)+N;
[check_PP,~]=ismember(Z(:,2)-N,R0);


%The third parameter:M
 M=zeros(N);
 PP=now_matching;
[row_P,colunm_P]=size(PP);
for i=1:row_P
    M(PP(i,1),PP(i,2))=1;
end
 
 % the main function of Markov chain
 
[ final_PP,eva_index] = check_Markov_sampling( m,Z,M);


if eva_index==1

    
%**********************replace the late updated bipartite graph***************************   
    
    
 R0=final_PP(:,1);
 match=final_PP;
TK{k,1}= match;

N=max(max(z));
M=size(z);
T=M(1,1);
zz=zeros(N);
yy=[];
r=1;
for i=1:T
zz(z(i,2),z(i,1))=1;
end
g=zz;
KK=[];
MKK1=[];
Q=1;
CC=C;

while Q
    
    
%R0 is the right targeted nodes of the updated bipartite graph in the N step
%unmatched nodes in the N step
m3=length(R0);
kK=1;E=[];
for i=1:length(z)
    [x1,y1]=ismember(z(i,2),R0);
    if x1~=0
        E(kK,:)=z(i,:);kK=kK+1;
    end    
end


%iteration condition
if sum(sum(E))==0
    Q=0;
end

if sum(sum(E))~=0
% *************Use Hungary algorithm,to calculate the maximum matching ***********************
Z=[];
N=length(g);
Z(:,1)=E(:,1);Z(:,2)=E(:,2)+N;

%*********************************************************************** 
%recalculate the match, and the purpose is to reduce the complexity
[ PP ] = transform_matching( Z,N );

%************************************************************************

M=zeros(N);
[oo1,oo2]=size(PP);
for i=1:oo1
    M(PP(i,1),PP(i,2))=1;
end


% [M,result] = Greedy_markov( F,CC,Z,E,g,N,NM );
   I=1;K=[];
    [a,b]=size(M);
for i=1:a
    [x,y]=find(M(i,:)~=0);
    if length(x)~=0
        K(I,1)=i; K(I,2)=y;
        I=I+1;
    end
end
Q=1;
m1=length(MKK1);
R0=K(:,1);%R0is the right target nides in the N step

CC=E(:,1);%left matched nodes

%unmatched nodes in the N+1 step
m4=length(R0);

TK{k+1,1}=K;

KK=[KK;K];
MKK=unique(KK,'rows');
MKK1=MKK;
k=k+1;
%avoid circle iteration
m2=length(MKK1);


if (m1==m2) && (m3==m4) 
    Q=0;
end

end


end

%delete the null element in TK
 Mark=[];
for i=1:length(TK)
       
               if length(TK{i,1})==0
                   Mark(i,1)=1;
               end
               if length(TK{i,1})~=0
                   Mark(i,1)=0;
               end
 end
    [sv1,sg1]=find(Mark==1);
    
    
   TK(sv1)=[];
  
 %**************************************************************************  
   
   
end



if eva_index==0
    TK=[];
    
end

end



