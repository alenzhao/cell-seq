function [ gene_names, gene_desc, gene_lengths ] = get_mm9_info(mm9_annot, mm9_list )
%GET_MM9_INFO Returns information about a gene (also converts names to systematic names)
%   Input:
%       mm9_list - List of gene ids
%   Output:
%       gene_names - common names
%       gene_desc - gene descriptions
%       gene_lengths - gene length


GeneSystematicName = cellfun(@char, mm9_annot.GeneSystematicName, 'UniformOutput', false);
GeneStandardName = cellfun(@char, mm9_annot.GeneStandardName, 'UniformOutput', false);
%GeneBriefDescription = mm9_annot.GeneBriefDescription;
GeneLength = mm9_annot.GeneLength;

% switch to single row cell (from single column cell)
if size(mm9_list,1)>1 && size(mm9_list,2)==1
    mm9_list = mm9_list';
end


list_len = numel(mm9_list);        % number of genes in query list

% create hash mapping systematic names to list index
id2ind = containers.Map(GeneSystematicName, 1:numel(GeneSystematicName));
name2ind = containers.Map(GeneStandardName, 1:numel(GeneStandardName));

% obtain index for each systematic name
known_gene_name_inds = isKey(id2ind, mm9_list);

query_inds = zeros(1,list_len);
for i = find(known_gene_name_inds)
    query_inds(i) = id2ind(mm9_list{i});
end

% pull out standard name for each index
gene_names(known_gene_name_inds) = GeneStandardName(query_inds(query_inds>0));

% if the standard name is missing, use inputed name
gene_names(~known_gene_name_inds) = mm9_list(~known_gene_name_inds);

for i = find(isKey(name2ind, mm9_list))
    query_inds(i) = name2ind(mm9_list{i});
end


% pull out short description for each index
%gene_desc = cell(1,list_len);
%gene_desc(query_inds>0) = GeneBriefDescription(query_inds(query_inds>0));

% pull out gene length for each gene
gene_lengths = zeros(1,list_len);
gene_lengths(query_inds>0) = GeneLength(query_inds(query_inds>0));


% use systematic names for unnamed genes
no_name_inds = strcmp('#N/A', gene_names);
gene_names(no_name_inds) = mm9_list(no_name_inds);
