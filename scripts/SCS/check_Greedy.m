function [ TK] = check_Greedy( z,C,NM )
%%
%Function:find the maximum matching in the linking and dyanmic graph
%Input:   z:the structure of network
%         C:target nodes
%         NM:the maximum iterted number
%Output:
%         TK:the maximum matching in the linking and dyanmic graph
%%
TK=cell(NM,1);
R0=C;
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
    
    
%R0is the right target nides in the N step
m3=length(R0);
k=1;E=[];
for i=1:length(z)
    [x1,y1]=ismember(z(i,2),R0);
    if x1~=0
        E(k,:)=z(i,:);k=k+1;
    end    
end


%the stop condition
if sum(sum(E))==0
    Q=0;
end

if sum(sum(E))~=0
% *************find the maximum mathcing in each iterted bipartite graph***********************
Z=[];
N=length(g);
Z(:,1)=E(:,1);Z(:,2)=E(:,2)+N;

[ PP ] = transform_matching( Z,N );

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
R0=K(:,1);%The left nodes on the N step K(:,1)£¬as the right mathced nodes in N+1 step

CC=E(:,1);%left matched nodes
%unmatched nodes in the N+1 step
m4=length(R0);

TK{r,1}=K;

KK=[KK;K];
MKK=unique(KK,'rows');
MKK1=MKK;
r=r+1
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


   
  
   
end

