function pre = rom_ED_pre(V, wbar, ops)
% ROM_ED_PRE  Precompute quantities to get E and D directly from ROM coefficients.
% With w = wbar + V a (V orthonormal, grid mean <.> = sum/N^2):
%   D(a) = <w^2>/Re   = (wbar'wbar + 2 (V'wbar)'a + a'a) / (N^2 Re)
%   E(a) = <psi w>/2  = (wbar'P wbar + 2 (V'P wbar)'a + a'(V'P V)a) / (2 N^2)
% where P = streamfunction operator (symmetric).  Use rom_ED(A, pre, r, Re).
pre.N2  = numel(wbar);
pre.ww  = wbar' * wbar;
pre.dv  = V' * wbar;
PV      = ops.psi(V);
pre.P   = V' * PV;
pre.pw  = PV' * wbar;
pre.wpw = wbar' * ops.psi(wbar);
end
