function rom = galerkin_rom(V, wbar, ops, Re)
% GALERKIN_ROM  Intrusive POD-Galerkin ROM of the Kolmogorov flow.
%
% Substituting w = wbar + V a into  dw/dt = f + lap(w)/Re + adv(w,w)  and
% projecting onto V gives the quadratic ROM
%     da/dt = c + A a + H q(a),   q(a) = [a_i a_j]_{i<=j}
% with Re-affine parts
%     c = c0 + c1/Re,  A = A0 + A1/Re,  H independent of Re.
% Needs the full-order operators (ops = kf_ops(N, n)): intrusive.
%
%   rom = galerkin_rom(V, wbar, ops, Re)   (Re optional; sets rom.c, rom.A)

r = size(V, 2);
[I, J] = quad_pairs(r);
rom.r  = r;   rom.I = I;   rom.J = J;   rom.type = 'POD-Galerkin';
rom.c0 = V' * (ops.f + ops.adv(wbar, wbar));                 % forcing + mean advection
rom.c1 = V' * ops.lap(wbar);                                 % mean viscous term
rom.A0 = V' * (ops.adv(wbar, V) + ops.adv(V, wbar));          % advection by / of the mean
rom.A1 = V' * ops.lap(V);                                    % viscous term of the modes
H = zeros(r, numel(I));
B = 150;                                                     % pairs per batch (memory)
for b0 = 1:B:numel(I)
    b = b0:min(numel(I), b0+B-1);
    T = ops.adv(V(:, I(b)), V(:, J(b)));
    off = I(b) ~= J(b);                                      % i<j: both orderings
    T(:, off) = T(:, off) + ops.adv(V(:, J(b(off))), V(:, I(b(off))));
    H(:, b) = V' * T;
end
rom.H = H;
if nargin >= 4 && ~isempty(Re), rom = rom_set_Re(rom, Re); end
end
