function [ PP ] = transform_matching( Z,N )
%Function£ºWill be renumbered Z , recalculate the amximum matching, the purpose is to reduce the amount of calculation
%Input£º   initial maximum matching Z
%Output£º  matrix denoting the maximum mathcing PP
ZZ=Z;
edge=[ZZ(:,1);ZZ(:,2)];
re_edge=unique(edge);
[row_ZZ,colunm_ZZ]=size(ZZ);
re_ZZ=[];
for i=1:row_ZZ
   [~,pp1]=ismember(ZZ(i,1),re_edge);
   re_ZZ(i,1)=pp1;
   [~,pp2]=ismember(ZZ(i,2),re_edge);
   re_ZZ(i,2)=pp2;
end
% [check_PP,~]=ismember(re_edge(re_ZZ(:,2))-N,R0)

%renumber The network nodes ,to reduce the storage and computation
BB=[];
[c,d]=size(re_ZZ);
for i=1:c
    BB(re_ZZ(i,1),re_ZZ(i,2))=1;
    BB(re_ZZ(i,2),re_ZZ(i,1))=1;
end
BB=sparse(BB);
% the maximum matching related with the targeted nodes
F = matching(BB);
%Apply the Markov sampling to identify different maximum matching
%M is the initial matching£¬result is the output of markov sampling denoting
%different maximum mathcing 
f=[];
N1=length(BB);
f(:,1)=[1:N1]';
f(:,2)=F(1:N1,1);
PP0=[];kk=1;
for i=1:length(f)
    if f(i,2)>f(i,1)
        PP0(kk,1)=f(i,1);PP0(kk,2)=f(i,2);
        kk=kk+1;
    end
end
PP=[re_edge(PP0(:,1)) re_edge(PP0(:,2))-N];



end

