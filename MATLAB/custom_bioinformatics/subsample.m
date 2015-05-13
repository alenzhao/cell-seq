function [ sub_counts ] = subsample( probs, depth )
%SUBSAMPLE sample reads from a probability distribution at a desired depth
% [ sub_counts ] = subsample( probs, depth )
% function [ sub_counts ] = subsample( probs, depth )
%   probs - probability matrix to sample from
%   depth - desired read count after sampling


sub_counts = zeros(size(probs));

parfor exp=1:size(probs,2)
%    samples = randsample(1:size(probs,1), depth, 'true', probs(:,exp));
%    sub_counts(:,exp) = histc(samples, 1:size(probs,1));
   sub_counts(:,exp) = mnrnd(depth, probs(:,exp));


end


end

