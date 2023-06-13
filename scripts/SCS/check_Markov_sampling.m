function [ final_PP,eva_index] = check_Markov_sampling( m,Z,M)
%%Function: The process of Markov sampling
%Input??
%m:initial maximum matching connections
%Z:roginal connections between nodes
%M??original markov chain
%Output??
%M??different maximum matching in the markov sampling
%eva_index:the index denotes whether we suceessfully find the new maximum matching
%%
final_PP=[];
[row_m,colunm_m]=size(m);
 M_new=M;
M1=cell(row_m,1);M1{1,1}=M; 
eva_index=0;
Z1=[];
[row_Z,colunm_Z]=size(Z);
for t=1:row_Z
  
    Z1(t,1)=Z(t,1);Z1(t,2)=Z(t,2);
    [er1,et1]=ismember(Z(t,1),m(:,1));
    if et1~=0
    [er2,et2]=ismember(Z(t,2),m(et1,2));
    er=et1*et2;
    end
    if et1==0
        er=et1;
    end
    if er~=0
        Z1(t,3)=-1;
    end
       if er==0
        Z1(t,3)=1;
       end
end
        
N=length(M);

[row_m,~]=size(m);
p2 = randperm(row_m)';
[row_p2,colunm_p2]=size(p2);

for iu=1:row_p2
   
    
  ip=p2(iu,1);%choose the address of matched edges
  mm=m;
 ZZ=Z;

%search for each new matched edge with length(m)??
[x1,y1]=find(Z(:,1)==m(ip,1));[x2,y2]=find(Z(:,2)==m(ip,2));
delete_left_node=m(ip,1);
delete_right_node=m(ip,2);

xx=intersect(x1,x2);

ZZ(xx,:)=[];
% [check_PP,~]=ismember(ZZ(:,2)-N,R0)

%adept the new strategy
if length(ZZ)==0
    break;
end

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

%renumber the network node to reduce the complexity
BB=[];
[c,d]=size(re_ZZ);
for i=1:c
    BB(re_ZZ(i,1),re_ZZ(i,2))=1;
    BB(re_ZZ(i,2),re_ZZ(i,1))=1;
end
BB=sparse(BB);
% find matching in iterated biaprtite graph
F = matching(BB);
%identify different matching 
f=[];
N1=length(BB);
f(:,1)=[1:N1]';
f(:,2)=F(1:N1,1);
PP=[];kk=1;
for i=1:length(f)
    if f(i,2)>f(i,1)
        PP(kk,1)=f(i,1);PP(kk,2)=f(i,2);
        kk=kk+1;
    end
end


final_PP=[re_edge(PP(:,1)) re_edge(PP(:,2))-N];
[row_PP,colunm_PP]=size(final_PP);
for i=1:row_PP
    %i
[check_PP,index]=ismember(final_PP(i,2),mm(:,2));
if index~=0
oi=final_PP(i,1)==mm(index,1);
if oi==0
     eva_index=1;
    break;
end
end

if index==0
     eva_index=1;
    break;
end

end


if i<row_PP-1
    eva_index=1;
    break;
end

end

%***************************************************************************        

end





