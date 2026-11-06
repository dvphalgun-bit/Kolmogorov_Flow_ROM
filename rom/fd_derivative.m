function [Ad, keep] = fd_derivative(A, dt)
% FD_DERIVATIVE  Fourth-order central difference in time of the columns of A.
%   [Ad, keep] = fd_derivative(A, dt): Ad(:,k) approximates dA/dt at column
%   keep(k); the first and last two columns are dropped.
m = size(A, 2);
keep = 3:m-2;
Ad = (-A(:, keep+2) + 8*A(:, keep+1) - 8*A(:, keep-1) + A(:, keep-2)) / (12*dt);
end
