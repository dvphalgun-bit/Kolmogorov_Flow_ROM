function [V, wbar, ev, info] = pod_basis_inv(W, rmax, stride, N)
% POD_BASIS_INV  Translation-invariant POD basis (x-shift symmetry).
%
% Kolmogorov flow is equivariant under shifts in x, so its patterns drift in x
% (Re = 20: a relative periodic orbit travelling ~0.12 rad per time unit).
% A standard POD basis only contains the x-positions seen in training and its
% projection error grows in the test window. The invariant POD is the POD of
% the training data together with ALL its x-shifted copies; it is computed
% exactly and cheaply wavenumber by wavenumber:
%   modes = phi_{kx,j}(y) e^{i kx x}  ->  real pairs  Re(.), Im(.)
% where phi_{kx,j} are eigenvectors of the y-correlation matrix of the kx-th
% Fourier coefficient. The mean is the x-average of the time mean.
%
%   [V, wbar, ev, info] = pod_basis_inv(W, rmax, stride, N)
%   V      N^2 x r orthonormal modes (r >= rmax, so that cos/sin pairs are complete)
%   wbar   mean field (depends on y only)
%   ev     modal energies (descending), one per mode
%   info.kx    wavenumber of each mode;  info.pc/ps/pk  cos/sin index pairs and kx
%   info.cuts  valid truncations r (pairs never split)

if nargin < 3 || isempty(stride), stride = 1; end
m = size(W, 2);
% x-averaged time mean
wm = zeros(N^2, 1);
for j0 = 1:500:m
    j1 = min(m, j0+499);  wm = wm + sum(double(W(:, j0:j1)), 2);
end
wm = reshape(wm / m, N, N);                      % [y, x]
wbar = reshape(repmat(mean(wm, 2), 1, N), [], 1);

X  = double(W(:, 1:stride:end)) - wbar;
ms = size(X, 2);
F  = fft(reshape(X, N, N, ms), [], 2) / N;       % Fourier in x:  (y, kx, snapshot)
clear X
nper = 30;                                       % candidates per wavenumber
E = [];  KX = [];  PHI = zeros(N, 0);
for kx = 0:N/2-1
    Fk = reshape(F(:, kx+1, :), N, ms);
    C  = (Fk * Fk') / ms;   C = (C + C') / 2;
    [P, L] = eig(C);
    [lam, o] = sort(real(diag(L)), 'descend');
    k = 1:min(nper, N);
    mult = 1 + (kx > 0);                         % +kx and -kx carry equal energy
    E   = [E; lam(k) * mult];                    %#ok<AGROW>
    KX  = [KX; kx * ones(numel(k), 1)];          %#ok<AGROW>
    PHI = [PHI, P(:, o(k))];                     %#ok<AGROW>
end
clear F
[E, o] = sort(E, 'descend');  KX = KX(o);  PHI = PHI(:, o);

x = (0:N-1) * 2*pi / N;
V = zeros(N^2, 0);  kxm = [];  ev = [];  pc = [];  ps = [];  pk = [];  cuts = [];
g = 0;
while size(V, 2) < rmax
    g = g + 1;  phi = PHI(:, g);  kx = KX(g);
    if kx == 0
        f = real(phi) * ones(1, N);
        V = [V, f(:) / norm(f(:))];  kxm = [kxm; 0];  ev = [ev; E(g)];        %#ok<AGROW>
    else
        Z = phi * exp(1i * kx * x);
        c = real(Z);  s = imag(Z);
        V = [V, c(:)/norm(c(:)), s(:)/norm(s(:))];                          %#ok<AGROW>
        kxm = [kxm; kx; kx];  ev = [ev; E(g)/2; E(g)/2];                     %#ok<AGROW>
        pc = [pc; size(V,2)-1];  ps = [ps; size(V,2)];  pk = [pk; kx];       %#ok<AGROW>
    end
    cuts = [cuts, size(V, 2)];                                               %#ok<AGROW>
end
info = struct('kx', kxm, 'pc', pc, 'ps', ps, 'pk', pk, 'cuts', cuts, ...
              'energy_total', sum(E) * 1, 'n_snap', ms);
end
