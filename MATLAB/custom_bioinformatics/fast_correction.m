function [corrected_pcs, correction, numerator, denominator] = fast_correction(pcs, eig_vals, noisy_pcs, noisy_eig_vals, correct_num)
%FAST_CORRECTION computes theoretical bound using speedup from low-rank noise
% [corrected_pcs, correction] = fast_correction(pcs, eig_vals, noisy_pcs, noisy_eig_vals, correct_num)
%   pcs - eigen values of deep data
%           -must be complete/full basis
%   eig_vals - eigenvales of deep data
%                 -1 for each vector in complete basis
%                 -should be scaled up for read depth (squared)
%                 -must be colun vector
%   noisy_pcs - PCs of noisy covariance matrix
%                 -incomplete set OK
%   correct_num - number of PCs to compute correction for



%%%% dimension of matrices
dim = numel(eig_vals);
true_dim = numel(noisy_eig_vals);

%%%% Numerator
L_hat = diag(noisy_eig_vals);
bases_mult = pcs'*noisy_pcs;
numerator = bases_mult*L_hat*bases_mult(1:correct_num, 1:true_dim)';
numerator(logical(eye(dim,correct_num))) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% Denominator
% compute pairwise differences of eigenvalues
denominator = repmat(eig_vals(1:correct_num)', dim, 1) - repmat(eig_vals, 1, correct_num);
% bring contribution of diagonal to 0
denominator(logical(eye(size(denominator)))) = Inf;
denominator = 1./denominator;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% Putting them together
eig_multipliers = numerator(:,1:correct_num).*denominator;
correction = pcs*eig_multipliers(:,1:correct_num);
corrected_pcs = normc(pcs(:, 1:correct_num) + correction);




end



