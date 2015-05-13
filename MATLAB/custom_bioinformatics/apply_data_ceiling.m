function probs = apply_data_ceiling(probs, percentile)
%APPLY_DATA_CEILING windserize data
%   probs - probability of reads to genes (gene x experiments)
%   percentile - percentile of data to be considered as ceiling


% upper_thresh = sort(probs(:), 'descend');
% upper_thresh = upper_thresh(floor(numel(upper_thresh)*percentile));
% probs(probs>upper_thresh) = upper_thresh;
% probs = counts_to_pdf(probs);

for i=1:size(probs,2)
   upper_thresh = sort(probs(:,i), 'descend');
   upper_thresh = upper_thresh(floor(numel(upper_thresh)*percentile));
   probs(probs(:,i)>upper_thresh,i) = upper_thresh;
end
probs = counts_to_pdf(probs);


end

