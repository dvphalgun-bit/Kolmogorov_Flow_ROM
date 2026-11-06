function [E, D] = rom_ED(A, pre, Re)
% ROM_ED  Kinetic energy and dissipation of ROM states (columns of A, r rows).
% Uses the first r entries of the precomputed quantities (POD prefix property).
r = size(A, 1);
D = (pre.ww + 2*pre.dv(1:r)'*A + sum(A.^2, 1)) / (pre.N2 * Re);
E = 0.5 * (pre.wpw + 2*pre.pw(1:r)'*A + sum(A .* (pre.P(1:r,1:r)*A), 1)) / pre.N2;
end
