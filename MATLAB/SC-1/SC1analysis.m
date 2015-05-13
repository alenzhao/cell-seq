% Sisi Chen
% 05/11/2015
% script for analyzing 28 cells from SC-1 dataset

load SC1.mat
datam = params.allsamples;
datagenes = params.gene_names';

%calculate correlation coefficients between paired ends
n = size(datam,2)/2;
if mod(n,1)
    error('not a paired end data set')
end

corrcoeff_m = []; 
avg_reads_m = [];

for i=1:n    
    correlation = corrcoef(datam(:,(2*i-1):2*i));
    corrcoeff_m = [corrcoeff_m, correlation(1,2)];
    curr_avg_reads = mean(params.aligned_read_num((2*i-1):2*i));
    avg_reads_m = [avg_reads_m, curr_avg_reads];
end

% check if there is a correlation between number of reads and corrcoeff_m
y = corrcoeff_m;
x = avg_reads_m;

scatter(x,y);hold on;

[p,S] = polyfit(x,y, 1);
yfit = p(1)*x + p(2);
yresid = y- yfit;
SSresid = sum(yresid.^2);
SStotal =(length(y)-1)*var(y);
rsq = 1-SSresid/SStotal;
rsq_adj = 1 - SSresid/SStotal * (length(y)-1)/(length(y)-length(p));

plot(x,yfit);
ylabel('correlation coefficient');
xlabel('average num of reads');
text(250000,0.75, ['y=' , num2str(p(1)) , ' * x + ', num2str(p(2)), ', R^2=', num2str(rsq)])

%remove zero entries in datam and gene names
[datam zero_indices] = remove_zero_read_genes(datam);
datagenes(zero_indices) = [];

%normalize data to read counts
dataprobs_m = counts_to_pdf(datam);

%split data into 5' reads and 3' reads
odds = (1:2:n*2);
evens = (2:2:n*2);
P5_data_m = remove_zero_read_genes(dataprobs_m(:,odds));
P7_data_m = remove_zero_read_genes(dataprobs_m(:,evens));

%sort data by read count
s1data_m = dataprobs_m(:,1);
[s1data_sorted, inds] = sortrows(s1data_m);
s1data_sorted = flipud(s1data_sorted);
s1genes_sorted = flipud(datagenes(inds));

s1data_P7_m = dataprobs_m(:,2);
[s1data_P7_sorted, inds] = sortrows(s1data_P7_m);
s1data_P7_sorted = flipud(s1data_P7_sorted);
s1genes_P7_sorted = flipud(datagenes(inds));

s13data_m = dataprobs_m(:,25);
[s13data_sorted, inds] = sortrows(s13data_m);
s13data_sorted = flipud(s13data_sorted);
s13genes_sorted = flipud(datagenes(inds));

s13data_P7_m = dataprobs_m(:,26);
[s13data_P7_sorted, inds] = sortrows(s13data_P7_m);
s13data_P7_sorted = flipud(s13data_P7_sorted);
s13genes_P7_sorted = flipud(datagenes(inds));

%use nmf to cluster data

[w5 h5] = nnmf(P5_data_m, 10, );
figure;
imagesc(h5);
title('nnmf p5');

[w7 h7] = nnmf(P7_data_m, 10);
figure;
imagesc(h7);
title('nnmf p7');
