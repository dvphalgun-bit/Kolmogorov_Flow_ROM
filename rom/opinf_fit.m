function [rom, info] = opinf_fit(Ahat, Adot, lam1, lam2, Reinv)
% OPINF_FIT  Non-intrusive Operator Inference of the quadratic ROM.
%
% Finds c, A, H in  da/dt = c + A a + H q(a)  from projected snapshots alone:
%   min sum_k || c + A a_k + H q(a_k) - adot_k ||^2
%       + g1 (||c||^2 + ||A||_F^2) + g2 ||H||_F^2
% solved as one stacked least-squares problem (QR, no normal equations).
% lam1, lam2 are dimensionless: g1 = lam1 * s_lin, g2 = lam2 * s_quad, where
% s_lin, s_quad are the mean squared column norms of the linear and quadratic
% blocks of the data matrix, so the same lambda grid works for any Re / r.
%
%   [rom, info] = opinf_fit(Ahat, Adot, lam1, lam2)          single Re
%   [rom, info] = opinf_fit(Ahat, Adot, lam1, lam2, Reinv)   Re-affine:
%       c = c0 + c1/Re, A = A0 + A1/Re (Reinv = 1/Re of each sample, 1 x K)
%
%   Ahat   r x K projected states,  Adot  r x K time derivatives
%   info.cond   condition number of the column-normalised data matrix
%   info.ratio  number of samples / number of unknowns per row
%   info.resid  relative training residual of the fit

[r, K] = size(Ahat);
[I, J] = quad_pairs(r);
Q = Ahat(I, :) .* Ahat(J, :);
param = nargin >= 5 && ~isempty(Reinv);
if param
    Reinv = Reinv(:)';
    D = [ones(1,K); Ahat; Q; Reinv; Ahat .* Reinv]';          % K x p
else
    D = [ones(1,K); Ahat; Q]';
end
nl = 1 + r;   nq = size(Q, 1);
cn2 = sum(D.^2, 1);
s_lin  = mean(cn2(2:nl));   s_quad = mean(cn2(nl+1:nl+nq));
reg = [lam1*s_lin*ones(1, nl), lam2*s_quad*ones(1, nq)];
if param, reg = [reg, lam1*mean(cn2(nl+nq+1:end))*ones(1, nl)]; end
p = size(D, 2);
O = [D; diag(sqrt(reg))] \ [Adot'; zeros(p, r)];               % p x r
Ot = O';
rom.r = r;   rom.I = I;   rom.J = J;   rom.type = 'OpInf';
rom.c0 = Ot(:, 1);   rom.A0 = Ot(:, 2:r+1);   rom.H = Ot(:, r+2:r+1+nq);
if param
    rom.c1 = Ot(:, r+2+nq);   rom.A1 = Ot(:, r+3+nq:end);
else
    rom.c1 = zeros(r, 1);     rom.A1 = zeros(r);
    rom = rom_set_Re(rom, Inf);   rom.Re = [];
end
cn = sqrt(sum(D.^2, 1));  cn(cn == 0) = 1;
s = svd(D ./ cn);
info.cond  = s(1) / s(end);
info.ratio = K / p;
info.p     = p;
info.resid = norm(D*O - Adot', 'fro') / norm(Adot, 'fro');
info.lam   = [lam1 lam2];
end
