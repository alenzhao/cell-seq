function [ enrichments, gene_hits ] = pathway_enrichments( gene_lists, goi_lists)
%PATHWAY_ENRICHMENT tests lists of genes for enrichment against pathways
%   There is a significant speedup if matlabpool is open, allowing the code
%   to use parfor loops
%
% [ enrichments, gene_hits ] = pathway_enrichments( gene_lists, goi_list)
%
% INPUTS:
%   'gene_lists'    - cell of gene lists of genes involved in each pathway 
%                               (KEGG pathways are in kegg_pathways_yeast.mat)
%   'goi_list'    - cell, lists (columns) of yeast ORF IDs to test for enrichment
%                               make sure that you use the right naming
%                               system
%                               yeast - yorf id
%                               mouse - symbolic name
% OUTPUTS:
%   'enrich_scores' - enrichment statistic against each pathway
%   'gene_hits'     - genes in query that are involved in each pathway



% you need an associated .mat file, eg:
% load ~/MATLAB/functions/custom_bioinformatics/kegg_pathways_yeast.mat;


% build hash table containing genes and unique indeces
path_genes = unique(vertcat(gene_lists{:}));
path_gene_map = containers.Map(path_genes, 1:numel(path_genes));

% number of genes with assocated path terms
num_genes_all_paths = numel(path_genes);
% number of genes in each pathway
num_genes_path = cellfun(@numel, gene_lists);
% number of different lists to query
num_queries = size(goi_lists,2);
% number of pathways in gene sets
num_pathways = numel(gene_lists);
% inds of pathways with >0 genes
nonzero_path_inds = cellfun(@numel, gene_lists)>0;

% initialize cell to collect genes that were found in pathways
gene_hits = cell(1, num_queries);
enrichments = zeros(num_pathways, num_queries);

% generate hash table for each pathway
path_hashes = cell(1,numel(nonzero_path_inds));
for path_ind = find(nonzero_path_inds)
    genes_in_path = gene_lists{path_ind};
    path_hashes{path_ind} = containers.Map(genes_in_path, 1:numel(genes_in_path));
end


% loop over all querying lists
for list_ind = 1:num_queries
    
    enrich_scores = zeros(1,num_pathways);
    gene_hits{list_ind} = cell(1,num_pathways);
    
    % current query list
    goi_list = goi_lists(:,list_ind);
    
    % find genes in query that are associated with any pathway
    key_inds = isKey(path_gene_map, goi_list);
    
    % update current query to only include genes in path
    goi_list = goi_list(key_inds);
    % number of genes (which are in any pathway) in the query
    num_genes_query = numel(goi_list);
    
    % zero the number of hits
    hit_count = zeros(1,num_pathways);
        
    in_path_inds = cell(1,numel(nonzero_path_inds));

    parfor path_ind = 1:numel(nonzero_path_inds)
        %if nonzero_path_inds(path_ind) == 1
            % count hits
            in_path_inds{path_ind} = isKey(path_hashes{path_ind}, goi_list);
            hit_count(path_ind) = sum(in_path_inds{path_ind});
            
            if hit_count(path_ind)>0
                % hypergeometric test for enrichment
                enrich_scores(path_ind) = 1 - hygecdf( ...
                    hit_count(path_ind)-1, ...
                    num_genes_all_paths, ...
                    num_genes_path(path_ind), ...
                    num_genes_query);
            end
        %end
    end
    
    % replace entries with no hits with a p-value of "1"
    enrich_scores(hit_count==0) = 1;
    enrichments(:,list_ind) = enrich_scores;
    
    % get genes assocated with each pathway
    for path_ind = nonzero_path_inds
        if nonzero_path_inds(path_ind) == 1
            gene_hits{list_ind}{path_ind} = goi_list(in_path_inds{path_ind});
        end
    end
    
    
end


%%%%%%%%%%%%%  OLD CODE
%         % hypergeometric test for enrichment    OLD (slightly slower)
%         enrich_scores(path_ind) = sum(hygepdf( ...
%             hit_count(path_ind):min(num_genes_path(path_ind), ...
%                                     num_genes_query), ...
%             num_genes_all_paths, ...
%             num_genes_path(path_ind), ...
%             num_genes_query));
end