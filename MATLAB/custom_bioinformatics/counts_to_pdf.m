function probs  = counts_to_pdf( counts )
%counts_to_pdf Summary of this function goes here
%   Detailed explanation goes here

[gene_num, ~] = size(counts);

probs = counts./repmat(sum(counts), gene_num, 1);

end


