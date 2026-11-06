function [err, proj] = rom_error(Ahat, Aproj, nw2)
% ROM_ERROR  Relative state error of a ROM trajectory, without rebuilding fields.
%   err(t)  = ||w'(t) - V a(t)|| / ||w'(t)||,   w' = w_FOM - wbar
%   proj(t) = ||w'(t) - V V' w'(t)|| / ||w'(t)||  (projection error, lower bound)
% Uses ||w' - V a||^2 = ||w'||^2 - 2 a.(V'w') + ||a||^2  (V orthonormal).
%   Ahat   r x m ROM coefficients (NaN after a blow-up)
%   Aproj  r x m projected FOM coefficients V'w'
%   nw2    1 x m squared norms ||w'||^2
e2 = nw2 - 2*sum(Ahat .* Aproj, 1) + sum(Ahat.^2, 1);
e2(e2 < 0) = 0;                                   % round-off; NaN stays NaN
err = sqrt(e2 ./ nw2);
p2 = nw2 - sum(Aproj.^2, 1);  p2(p2 < 0) = 0;
proj = sqrt(p2 ./ nw2);
end
