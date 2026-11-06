function ops = kf_ops(N, n)
% KF_OPS  Spectral operators of the Kolmogorov-flow full-order model.
%
% Exactly the same discretisation as kolmogorov_solver.m (Fourier, 2/3-rule
% dealiasing), but written for fields stored as columns of an N^2 x m matrix
% (column j = reshape(w, [], 1), rows = y, columns = x, as in the snapshots).
%
%   ops = kf_ops(N, n)
%   ops.f            forcing vector  -n cos(n y)                 (N^2 x 1)
%   ops.lap(W)       Laplacian of each column
%   ops.psi(W)       streamfunction: lap(psi) = -w
%   ops.adv(W1, W2)  dealiased -u(W1).grad(W2)   (W1 may be one column)
%   ops.rhs(W, Re)   full right-hand side  f + lap(W)/Re + adv(W, W)
%   ops.E(W)         kinetic energy   0.5 <|u|^2>   of each column
%   ops.D(W, Re)     dissipation rate <w^2>/Re      of each column
%
% Used by the POD-Galerkin ROM (intrusive) and by the diagnostics.
% Operator Inference does NOT use this file.

L = 2*pi;  x = (0:N-1)' * L/N;
[~, Y] = meshgrid(x, x);
k = [0:N/2-1, -N/2:-1]';
[KX, KY] = meshgrid(k, k);
K2 = KX.^2 + KY.^2;
K2i = 1 ./ K2;  K2i(1,1) = 0;
mask = (abs(KX) < N/3) & (abs(KY) < N/3);

ops.N = N;  ops.n = n;  ops.x = x;
ops.f = reshape(-n * cos(n*Y), [], 1);

ops.lap = @(W) col(real(ifft2p(-K2 .* fft2p(pg(W, N)))));
ops.psi = @(W) col(real(ifft2p( K2i .* fft2p(pg(W, N)))));
ops.adv = @(W1, W2) adv(W1, W2, N, KX, KY, K2i, mask);
ops.rhs = @(W, Re) ops.f + ops.lap(W)/Re + ops.adv(W, W);
ops.E   = @(W) 0.5 * mean(ops.psi(W) .* W, 1);     % <|u|^2> = <psi w> on a torus
ops.D   = @(W, Re) mean(W.^2, 1) / Re;
end

% ---------------------------------------------------------------- helpers
function A = adv(W1, W2, N, KX, KY, K2i, mask)
    w1h = fft2p(pg(W1, N));   w2h = fft2p(pg(W2, N));
    ph  = w1h .* K2i;
    u   = real(ifft2p( 1i*KY .* ph));
    v   = real(ifft2p(-1i*KX .* ph));
    wx  = real(ifft2p( 1i*KX .* w2h));
    wy  = real(ifft2p( 1i*KY .* w2h));
    A   = col(real(ifft2p(-fft2p(u.*wx + v.*wy) .* mask)));   % broadcasting over pages
end
function P = pg(W, N),  P = reshape(W, N, N, []);  end          % columns -> pages
function W = col(P),    W = reshape(P, [], size(P, 3));  end    % pages -> columns
function F = fft2p(P),  F = fft(fft(P, [], 1), [], 2);  end
function P = ifft2p(F), P = ifft(ifft(F, [], 1), [], 2); end
