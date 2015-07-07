% script for plotting the number of reads to the number of barcodes

filename = 'matrix.txt';
fid = fopen(filename);
data = textscan(fid,'%s%f');
barcodes = data{1};
counts = [data{2}];
numbarcodes = size(counts,1);

% generate bins for the histogram
Y = hist(counts,[1:100:10000]);
hist(counts,[1:100:10000]);
 
% rank the barcodes in order
[counts_s, indices] = sort(counts);
barcodes_s = {barcodes{indices}};

% plot the number of reads
plot(counts_s);

% plot the cumsum of the reads: 
counts_sum = cumsum(counts_s);
rev_counts_sum = cumsum(flipud(counts_s));

h = figure;
plot(rev_counts_sum(1:1500));
xlabel('sorted barcodes');
ylabel('read counts');
hgsave('read_cumsum.fig')

% % fit the curve to an exponential
% 
% f1 = fit([1:1500]', rev_counts_sum(1:1500), 'exp2');
% figure;plot(rev_counts_sum(1:1500));hold on;plot(f1);
% 
% % try to find the symbolic inflection point, but there isn't really an
% % inflection point in the curve. 
% 
% syms x
% f_sym = f.a*exp(f.b*x) + f.c*exp(f.d*x);
% diff2 = diff(diff(f_sym));
% pt=solve(diff2);
% 

