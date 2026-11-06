function [e_mean, ok, curves, ec_mean, ccurves] = rom_short_error(rom, Aproj, nw2, t, starts, horizon, amax, arms)
% ROM_SHORT_ERROR  Mean short-horizon prediction error from several start times.
% For chaotic flow a single long trajectory error saturates after a few
% Lyapunov times, so accuracy is measured as the average error over
% [t0, t0 + horizon] for several start times t0 (indices 'starts' into t).
%   e_mean   mean FIELD error (vs the FOM field; includes the projection error)
%   ec_mean  mean COEFFICIENT error ||a_rom - a_proj|| / arms (dynamics error only)
%   ok       true if no run blew up (blown-up entries count as error 2)
if nargin < 8 || isempty(arms), arms = sqrt(mean(sum(Aproj.^2, 1))); end
dt = t(2) - t(1);
nstep = round(horizon / dt);
curves = nan(numel(starts), nstep + 1);  ccurves = curves;
ok = true;
for s = 1:numel(starts)
    k = starts(s) + (0:nstep);
    k = k(k <= numel(t));
    [Ah, st] = rom_integrate(rom, Aproj(:, k(1)), t(k), amax, 1e-6);
    if ~strcmp(st, 'ok'), ok = false; end
    curves(s, 1:numel(k))  = rom_error(Ah, Aproj(:, k), nw2(k));
    ccurves(s, 1:numel(k)) = sqrt(sum((Ah - Aproj(:, k)).^2, 1)) / arms;
end
c1 = curves;  c1(isnan(c1)) = 2;   e_mean  = mean(c1(:));
c2 = ccurves; c2(isnan(c2)) = 2;   ec_mean = mean(c2(:));
end
