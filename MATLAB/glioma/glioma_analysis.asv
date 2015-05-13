
%load data file
load glioma.mat

%first thing to do is to normalize all columns by total count to get
%probabilities

counts = params.allsamples;
sumtots = sum(counts, 1);
num_genes = size(counts,1);
probs = counts./repmat(sumtots, gene_num, 1);

%filter out genes with zero reads 

