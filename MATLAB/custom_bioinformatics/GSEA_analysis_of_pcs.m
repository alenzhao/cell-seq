function [terms_pos, terms_neg] = GSEA_analysis_of_pcs(pcs, gene_symbols, gene_lists)
%GSEA_ANALYSIS_OF_PCS compute GSEA enrichment on pcs
% [terms_pos, terms_neg] = GSEA_analysis_of_pcs(pcs, gene_symbols, gene_lists)
%   terms_pos - -log10 of p-values for positive genes in the pc
%                 - matrix has dimensions terms x pcs
%   terms_neg - -log10 of p-values for negative genes in the pc
%                 - matrix has dimensions terms x pcs

% set significance threshold for being included in query list
zscore_thresh = 2;

pos_sig_genes = zscore(pcs)>zscore_thresh;
neg_sig_genes = zscore(pcs)<-zscore_thresh;

if (size(gene_symbols,1) == 1) && (size(gene_symbols,2) > 1)
   gene_symbols = gene_symbols';
end


for pc_ind = 1:size(pcs,2)
   % postive part of the pc
   query_genes = gene_symbols(pos_sig_genes(:,pc_ind));
	p_vals_pos(:,pc_ind) = pathway_enrichments(gene_lists, query_genes);

   % negative part of the pc
   query_genes = gene_symbols(neg_sig_genes(:,pc_ind));
	p_vals_neg(:,pc_ind) = pathway_enrichments(gene_lists, query_genes);     
end

terms_pos = -log10(p_vals_pos);
terms_neg = -log10(p_vals_neg);


end

