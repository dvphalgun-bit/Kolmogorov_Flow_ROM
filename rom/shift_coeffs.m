function A2 = shift_coeffs(A, info, theta)
% SHIFT_COEFFS  Coefficients of the x-shifted field w(x - theta) in the
% translation-invariant POD basis (pod_basis_inv).
% Each (cos, sin) pair with wavenumber kx is rotated by the angle kx*theta;
% kx = 0 modes are unchanged. Because the Kolmogorov equations are shift-
% equivariant, a shifted trajectory is again an exact trajectory, so shifted
% copies of the training data are valid extra training data (augmentation).
%   A      r x K coefficients (first r modes of the basis)
%   theta  shift in x (radians)
r = size(A, 1);
A2 = A;
use = info.pc <= r & info.ps <= r;
for q = find(use)'
    i = info.pc(q);  j = info.ps(q);  ph = info.pk(q) * theta;
    A2(i, :) = cos(ph) * A(i, :) - sin(ph) * A(j, :);
    A2(j, :) = sin(ph) * A(i, :) + cos(ph) * A(j, :);
end
end
