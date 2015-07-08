%% input data from hundred stamps, mouse: 

%% Initialize variables.
filename = '/home/iamcam/Documents/cell-seq/MATLAB/inDrop/GSM1599494_ES_d0_main.csv';

%% Format string for each line of text:
g_filename = '/home/iamcam/Documents/cell-seq/MATLAB/inDrop/GSM1599494_ES_d0_genes.csv';
formatSpec = '%s';
fileID = fopen(g_filename,'r');
startRow=0;
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
m_datam = csvread(filename);
m_genes = dataArray{1};

%% Clear out unused data
clear dataArray fileID formatSpec delimter startRow

%% Remove zero read genes: 

[m_datam zero_indices] = remove_zero_read_genes(m_datam);
m_genes(zero_indices) = [];

%% Normalize data to read counts
m_dataprob = counts_to_pdf(m_datam);

% create summed data
m_dataprob_sum = sum(m_dataprob,2);
[s1data_sorted, inds] = sortrows(m_dataprob_sum);
s1data_sorted = flipud(s1data_sorted);
s1genes_sorted = flipud(m_genes(inds));

%% nnmf

[w, h] = nnmf(m_dataprob, 10 );
figure;
imagesc(h);
title('nnmf mouse inDrop');