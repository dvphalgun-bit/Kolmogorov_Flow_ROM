function [Ahat, status, tstop, cpu] = rom_integrate(rom, a0, teval, amax, tol)
% ROM_INTEGRATE  Integrate the quadratic ROM with ode45 and detect blow-up.
%   [Ahat, status, tstop, cpu] = rom_integrate(rom, a0, teval, amax, tol)
%   teval   output times (first entry = initial time)
%   amax    blow-up threshold on ||a|| (e.g. 10x the largest training norm)
%   tol     ode45 RelTol (default 1e-8; AbsTol = tol/100)
%   Ahat    r x numel(teval); columns after a blow-up are NaN
%   status  'ok' or 'blowup';  tstop = time of blow-up (Inf if none)
if nargin < 5 || isempty(tol), tol = 1e-8; end
f   = @(t, a) rom_rhs(rom, a);
ev  = @(t, a) blowup_event(t, a, amax);
opt = odeset('RelTol', tol, 'AbsTol', tol/100, 'Events', ev);
t0 = tic;
[tt, aa, te] = ode45(f, teval, a0, opt);
cpu = toc(t0);
Ahat = nan(numel(a0), numel(teval));
nk = numel(tt);
if numel(teval) == 2 && nk > 2                 % ode45 returns all steps for 2-point tspan
    Ahat(:, 1) = aa(1, :)';   Ahat(:, 2) = aa(end, :)';   nk = 2 * (tt(end) >= teval(end));
else
    Ahat(:, 1:nk) = aa';
end
if isempty(te)
    status = 'ok';  tstop = Inf;
else
    status = 'blowup';  tstop = te(1);
    Ahat(:, teval > tstop) = NaN;
end
end

function [value, isterminal, direction] = blowup_event(~, a, amax)
% event: stop the integration when ||a|| reaches amax
value = norm(a) - amax;   isterminal = 1;   direction = 0;
end
