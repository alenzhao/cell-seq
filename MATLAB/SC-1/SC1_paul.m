% 06/18/2015
% script for analyzing 28 cells from SC-1 dataset, paul's pre-processing
% and alignment

load SC1_paul.mat
datam = params.allsamples;
datagenes = params.gene_names';

%calculate correlation coefficients between paired ends
n = size(datam,2);

corrcoeff_m = []; 
avg_reads_m = [];

for i=1:n    
    curr_avg_reads = params.aligned_read_num(i);
    avg_reads_m = [avg_reads_m, curr_avg_reads];
end

%remove zero entries in datam and gene names
[datam zero_indices] = remove_zero_read_genes(datam);
datagenes(zero_indices) = [];

%normalize data to read counts
dataprobs_m = counts_to_pdf(datam);

%clip values for each gene to top 5% and bottom 5%
tailperc = 5;
topclip = prctile(dataprobs_m', 100-tailperc);
sum(topclip==0)
botclip = prctile(dataprobs_m', tailperc);
sum(botclip==0)

dataprobs_mc=zeros(size(dataprobs_m));

for i=1:size(dataprobs_m,1)
    curr_row = dataprobs_m(i,:);
    curr_row(curr_row>topclip(i))=topclip(i);
    curr_row(curr_row<botclip(i))=botclip(i);
    dataprobs_mc(i,:) = curr_row;
end

%sort data by read count
%first sum all data in samples
s1_ind=[2:13];
s2_ind=[14:25];
s3_ind=[26:29];

s1=sum(dataprobs_mc(:,s1_ind),2);
s2=sum(dataprobs_mc(:,s2_ind),2);
s3=sum(dataprobs_mc(:,s3_ind),2);

[s1_sorted, inds_1] = sortrows(s1);
[s2_sorted, inds_2] = sortrows(s2);
[s3_sorted, inds_3] = sortrows(s3);

s1genes_sorted = flipud(datagenes(inds_1));
s2genes_sorted = flipud(datagenes(inds_2));
s3genes_sorted = flipud(datagenes(inds_3));

%cluster analysis

y= pdist(dataprobs_m', 'euclidean');
tree = linkage(y);
dendrogram(tree);

%use nmf to cluster data
[w,h] = nnmf(dataprobs_mc(:,2:29), 10);
figure;
imagesc(h);
title('nnmf SC1 paul');
xlabel('samples');

%from the graph, it looks like w1 and w2 are the most variable between 6h
%and 3day samples. let's find the top 20 genes in those vectors and see
%what they are: 

%find top 20 genes in w1: 
[w1_sorted, w1_inds] = sortrows(w(:,1));
w1_genes_sorted=flipud(datagenes(w1_inds));
[w2_sorted, w2_inds] = sortrows(w(:,2));
w2_genes_sorted=flipud(datagenes(w2_inds));
[w7_sorted, w7_inds] = sortrows(w(:,7));
w7_genes_sorted=flipud(datagenes(w7_inds));
[w9_sorted, w9_inds] = sortrows(w(:,9));
w9_genes_sorted=flipud(datagenes(w9_inds));
[w3_sorted, w3_inds] = sortrows(w(:,3));
w3_genes_sorted=flipud(datagenes(w3_inds));
[w10_sorted, w10_inds] = sortrows(w(:,10));
w10_genes_sorted=flipud(datagenes(w10_inds));
