function s = burst_stats(t, Dn, thr)
% BURST_STATS  Statistics of the normalised dissipation Dn = D/D_lam over time t.
%   s.mean, s.std, s.max      moments of Dn
%   s.rate                    bursts (upward crossings of thr) per 100 time units
%   s.frac                    fraction of time with Dn > thr
%   s.fdom, s.pdom            dominant frequency of Dn and its share of power
if nargin < 3, thr = 0.15; end
ok = ~isnan(Dn);  t = t(ok);  Dn = Dn(ok);
s.mean = mean(Dn);  s.std = std(Dn);  s.max = max(Dn);
up = find(Dn(1:end-1) <= thr & Dn(2:end) > thr);
T = t(end) - t(1);
s.rate = 100 * numel(up) / max(T, eps);
s.frac = mean(Dn > thr);
d = Dn - mean(Dn);
if numel(d) > 16 && std(d) > 1e-8 * abs(mean(Dn))
    P = abs(fft(d)).^2;  P = P(2:floor(end/2));
    fr = (1:numel(P)) / (numel(d) * (t(2) - t(1)));
    [pm, im] = max(P);  s.fdom = fr(im);  s.pdom = pm / sum(P);
else
    s.fdom = NaN;  s.pdom = NaN;
end
end
