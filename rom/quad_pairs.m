function [I, J] = quad_pairs(r)
% QUAD_PAIRS  Index pairs (i <= j) of the r(r+1)/2 unique quadratic terms a_i a_j.
% The SAME ordering is used by the Galerkin operator H, the Operator Inference
% data matrix and the ROM right-hand side, so both ROMs have identical form.
[I, J] = find(triu(ones(r)));
end
