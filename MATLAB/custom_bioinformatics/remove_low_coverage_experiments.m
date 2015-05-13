function [counts, low_depth_inds] = remove_low_coverage_experiments(counts, required_read_depth)
%REMOVE_LOW_COVERAGE_EXPERIMENTS 
%   counts - raw read counts over all experments (genes x experiments)
%   required_read_depth - minimum allowable read depth

low_depth_inds = sum(counts)<required_read_depth;
counts(:,low_depth_inds) = [];

end

