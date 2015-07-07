%% input data from hundred stamps, mouse: 

%% Initialize variables.
filename = '/home/iamcam/Documents/cell-seq/MATLAB/inDrop/GSM1599494_ES_d0_main.csv';
g_filename = '/home/iamcam/Documents/cell-seq/MATLAB/inDrop/GSM1599494_ES_d0_genes.csv';


%% Format string for each line of text:
 
formatSpec = '%s'

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
m_datam = [dataArray{2:end-1}];
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
title('nnmf mouse drop-test 100 STAMPs');