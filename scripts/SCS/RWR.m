function [ q_new ] = RWR( A,mutation )
%Function:random walk to obtain the similarity between other nodes and mutation
%genes
%Input:
%     A£ºthe adjacency matrix of the network
%     mutation£ºthe num address of the mutation gene
%Column normalization of the matrix
leg=sum(A);
W=bsxfun(@rdivide,A,leg);%Normalized adjacency matrix
W(isnan(W))=0;

r=0.6;index=1;
initial=zeros(length(A),1);
initial(mutation,1)=1/length(mutation);
q0=initial;q_old=q0;
num=1;
while index
    
    
    q_new=(1-r)*W*q_old+r*q0;
    if norm(q_old-q_new,1)<10^(-18)
        index=0;
    end 
    q=q_new+q_old;
    q_old=q_new;
    num=num+1;
    if num>100
        index=0;
    end
    
end
 

end

