function [ counts, zero_inds ] = remove_zero_read_genes( counts )
%REMOVE_ZERO_READ_GENES 
%   counts - raw read counts over all experments (genes x experiments)

genes_l1_norm = sum(counts,2);
zero_inds = genes_l1_norm==0;
counts(zero_inds, :) = [];

end

