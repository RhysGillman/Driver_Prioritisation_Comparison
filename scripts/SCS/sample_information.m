function [ CNV_name,CNV_data,SNP_name,SNP_data,EXPR_name,EXPR_data ] = sample_information( node0,expression_fileName, CNV_fileName,SNP_fileName)
%Function:main function of our SCS
%Input:
%      expression_fileName:expression data 
%      CNV_fileName:copy number of genes 
%      SNP_fileName:point mutations data
%Output:
%     Constrained_B:the mutation number
%     Targets_C:differentially expressed genes number
%CNV data
name_gene=node0;
data=importdata(CNV_fileName);
tumor_gene=data.textdata(:,1);
fid=fopen(CNV_fileName,'r');
[~,memeber]=ismember(name_gene,tumor_gene);
fidw=fopen('CNV_internate.txt','w');
n=length(data.textdata);
for i=1:n
a{i}=fgetl(fid);%read the row line
end 
fclose(fid);
for i=1:n 
     if i==1
        fprintf(fidw,'%s\n',a{i});
    end
    if ismember(i,memeber)~=0
       fprintf(fidw,'%s\n',a{i});
    end
end
fclose(fidw);
%reload the data
[original,~,~]=importdata('CNV_internate.txt');
CNV_data=original.data;
CNV_name=original.textdata(2:end,1);
%write the line related with the gene in PPI network EXPR_internal.txt
%EXPR
data=importdata(expression_fileName);
tumor_gene=data.textdata(:,1);
[~,memeber]=ismember(name_gene,tumor_gene);
fid=fopen(expression_fileName,'r');
fidw=fopen('EXPR_internate.txt','w');
n=length(data.textdata);
for i=1:n
a{i}=fgetl(fid);
end 
fclose(fid);
for i=1:n 
     if i==1
        fprintf(fidw,'%s\n',a{i});
    end
    if ismember(i,memeber)~=0
       fprintf(fidw,'%s\n',a{i});
    end
end
fclose(fidw);
[original,~,~]=importdata('EXPR_internate.txt');
EXPR_data=original.data;
EXPR_name=original.textdata(2:end,1);


%SNP
data=importdata(SNP_fileName);
tumor_gene=data.textdata(:,1);
[~,memeber]=ismember(name_gene,tumor_gene);
fid=fopen(SNP_fileName,'r');
fidw=fopen('SNP_internate.txt','w');
n=length(data.textdata);
for i=1:n
a{i}=fgetl(fid);
end 
fclose(fid);
for i=1:n 
     if i==1
        fprintf(fidw,'%s\n',a{i});
    end
    if ismember(i,memeber)~=0
       fprintf(fidw,'%s\n',a{i});
    end
end
fclose(fidw);
[original,~,~]=importdata('SNP_internate.txt');
SNP_data=original.data;
SNP_name=original.textdata(2:end,1);


end

