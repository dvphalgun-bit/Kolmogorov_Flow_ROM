function [V, wbar, sv] = pod_basis(W, rmax, stride)
% POD_BASIS  Mean-subtracted POD basis of a snapshot matrix.
%   [V, wbar, sv] = pod_basis(W, rmax, stride)
%   W      N^2 x m snapshot matrix (single or double)
%   rmax   number of modes returned
%   stride use every stride-th snapshot for the SVD (default 1); with
%          dt_save = 0.1, stride = 5 gives dt = 0.5 spacing, which spans the
%          same subspace while keeping the SVD cheap
%   V      N^2 x rmax orthonormal POD modes, wbar = time-mean field,
%   sv     all singular values of the (subsampled) fluctuation matrix
if nargin < 3, stride = 1; end
m = size(W, 2);
wbar = zeros(size(W,1), 1);
for j0 = 1:500:m                                   % mean over ALL snapshots, chunked
    j1 = min(m, j0+499);
    wbar = wbar + sum(double(W(:, j0:j1)), 2);
end
wbar = wbar / m;
X = double(W(:, 1:stride:end)) - wbar;
[U, S, ~] = svd(X, 'econ');
sv = diag(S);
V = U(:, 1:min(rmax, size(U,2)));
end
