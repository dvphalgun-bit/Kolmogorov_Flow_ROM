function P = run_rom_param(cfg)
% RUN_ROM_PARAM  Step 8: prediction at Reynolds numbers not seen in training.
%
% One global POD basis is built from the training windows of all training Re
% (equal number of snapshots per Re, global mean = average of the per-Re means).
%   POD-Galerkin: the Re dependence is exact, c = c0 + c1/Re, A = A0 + A1/Re.
%   OpInf:        the same Re-affine form is learned from all training Re at
%                 once (opinf_fit with 1/Re features); lambda is chosen on the
%                 training Re only (bounded at every training Re, smallest mean
%                 short-horizon error).
% Both ROMs are then run at the unseen Re and compared with the FOM.
%
% cfg: N, n, train_files, train_Re, t_train (rows), pod_counts, deriv_rows,
%      n_aug, test_files, test_Re, r_list, lam1s, lam2s, horizon, Hv, nstarts,
%      tv_thr, tol, outfile
% Basis: translation-invariant POD (pod_basis_inv); OpInf data augmented with
% n_aug x-shifted copies (shift_coeffs), as in run_rom_case.

ops = kf_ops(cfg.N, cfg.n);
nt = numel(cfg.train_Re);
say(cfg, '\n===== parametric ROM: train Re = %s, test Re = %s =====\n', ...
    mat2str(cfg.train_Re), mat2str(cfg.test_Re));

% ---------------------------------------------------------------- global POD
% translation-invariant POD of an equal number of training snapshots per Re
% (so the mean is the equal-weight average of the per-Re means)
X = [];  rec = cell(1, nt);
for q = 1:nt
    S = load(cfg.train_files{q}, 'W', 't_snap', 'diagn');
    tr = S.t_snap >= cfg.t_train(q, 1) & S.t_snap <= cfg.t_train(q, 2);
    idx = find(tr);  st = max(1, floor(numel(idx) / cfg.pod_counts));
    idx = idx(1:st:end);  idx = idx(1:min(end, cfg.pod_counts));
    X = [X, S.W(:, idx)];                                        %#ok<AGROW>
    rec{q} = struct('W', S.W, 't', S.t_snap, 'train', tr, 'diagn', S.diagn);
    clear S
end
rmax = max(cfg.r_list);
[V, wbar, ev, binfo] = pod_basis_inv(X, rmax, 1, cfg.N);  clear X
sv = sqrt(ev);
r_list = arrayfun(@(r) min(binfo.cuts(binfo.cuts >= r)), cfg.r_list);
pre = rom_ED_pre(V, wbar, ops);
rng(0);  thetas = 2*pi*rand(1, cfg.n_aug);

% project every record onto the global basis
for q = 1:nt
    [rec{q}.A, rec{q}.nw2] = project_chunked(rec{q}.W, V, wbar);
    rec{q} = rmfield(rec{q}, 'W');
end
tests = cell(1, numel(cfg.test_Re));
for q = 1:numel(cfg.test_Re)
    S = load(cfg.test_files{q}, 'W', 't_snap', 'diagn');
    [A, nw2] = project_chunked(S.W, V, wbar);
    tests{q} = struct('A', A, 'nw2', nw2, 't', S.t_snap, 'diagn', S.diagn);
    clear S
end

P = struct('cfg', cfg, 'sv', sv, 'V', V, 'wbar', wbar, 'binfo', binfo, 'r_list', r_list);
for ir = 1:numel(r_list)
    r = r_list(ir);
    % ------------------------------------------------ training data for OpInf
    Atr = [];  Ad = [];  Rinv = [];  amax = 0;  e2tr = zeros(1, nt);  armq = e2tr;  vtr = e2tr;
    for q = 1:nt
        dt = rec{q}.t(2) - rec{q}.t(1);
        Aq = rec{q}.A(1:r, rec{q}.train);
        [Adq, keep] = fd_derivative(Aq, dt);  Aq = Aq(:, keep);
        % x-shifted copies (exact trajectories by symmetry), then equal rows per Re
        Xa = Aq;  Xd = Adq;
        for th = thetas
            Xa = [Xa, shift_coeffs(Aq, binfo, th)];  Xd = [Xd, shift_coeffs(Adq, binfo, th)];  %#ok<AGROW>
        end
        rs = max(1, floor(size(Xa, 2) / cfg.deriv_rows));
        Atr = [Atr, Xa(:, 1:rs:end)];  Ad = [Ad, Xd(:, 1:rs:end)];   %#ok<AGROW>
        Rinv = [Rinv, ones(1, numel(1:rs:size(Xa,2))) / cfg.train_Re(q)];  %#ok<AGROW>
        amax = max(amax, 10 * max(sqrt(sum(Aq.^2, 1))));
        e2tr(q) = mean(sum(Aq.^2, 1));  armq(q) = sqrt(e2tr(q));
        vtr(q)  = mean(sum((Aq - mean(Aq, 2)).^2, 1));
    end
    % ------------------------------------------------ Galerkin (exact Re dependence)
    t0 = tic;  G = galerkin_rom(V(:, 1:r), wbar, ops);  tG = toc(t0);
    % ------------------------------------------------ OpInf: select lambda on training Re
    t0 = tic;  best = [];  ebest = Inf;  infoO = struct();
    for i = 1:numel(cfg.lam1s)
        for j = 1:numel(cfg.lam2s)
            [O, inf_ij] = opinf_fit(Atr, Ad, cfg.lam1s(i), cfg.lam2s(j), Rinv);
            % same score as opinf_select (coefficient error + energy mismatch),
            % averaged over the training Re; any blow-up rejects the pair
            e = 0;
            for q = 1:nt
                Oq = rom_set_Re(O, cfg.train_Re(q));
                tq = rec{q}.t;  it = find(rec{q}.train);
                tfq = tq(it(1)) + (0:1:tq(end)-tq(it(1)));
                [Af, st] = rom_integrate(Oq, rec{q}.A(1:r, it(1)), tfq, amax, 1e-6);
                if ~strcmp(st, 'ok'), e = Inf; break; end
                Af = Af(:, tfq > tfq(1) & tfq <= tq(it(end)));   % training window only
                e_stat = abs(log(mean(sum(Af.^2, 1)) / e2tr(q))) + ...
                         abs(log(mean(sum((Af - mean(Af, 2)).^2, 1)) / vtr(q)));
                starts = round(linspace(it(1), it(end) - round(cfg.horizon/(tq(2)-tq(1))), cfg.nstarts));
                [~, ok, ~, ec] = rom_short_error(Oq, rec{q}.A(1:r, :), rec{q}.nw2, tq, starts, cfg.horizon, amax, armq(q));
                if ~ok, e = Inf; break; end
                e = e + (ec + e_stat) / nt;
            end
            if e < ebest, ebest = e; best = O; infoO = inf_ij; infoO.lam = [cfg.lam1s(i) cfg.lam2s(j)]; end
        end
    end
    tO = toc(t0);
    say(cfg, 'r = %2d: Galerkin built in %.1f s; OpInf lam = %s, cond %.2e, rows/unknowns %.1f, fit resid %.3f (%.0f s)\n', ...
        r, tG, mat2str(infoO.lam), infoO.cond, infoO.ratio, infoO.resid, tO);
    % ------------------------------------------------ evaluate at every Re
    allRe = [cfg.train_Re, cfg.test_Re];
    for q = 1:numel(allRe)
        if q <= nt, d = rec{q}; else, d = tests{q - nt}; end
        Re = allRe(q);  Dlam = Re / (2*cfg.n^2);
        t = d.t;  dt = t(2) - t(1);
        Dref = interp1(d.diagn.t_diag, d.diagn.D, t) / Dlam;
        A = d.A(1:r, :);  arms = sqrt(mean(sum(A.^2, 1)));
        Hn = round(cfg.Hv / dt);
        starts = round(linspace(1, numel(t) - Hn, cfg.nstarts));
        pc = nan(numel(starts), Hn + 1);
        for s = 1:numel(starts)
            k = starts(s) + (0:Hn);
            pc(s, :) = sqrt(sum((A(:, k(1)) - A(:, k)).^2, 1)) / arms;
        end
        P.pers{ir, q} = struct('coef_short', mean(rowmean(pc(:, 1:round(cfg.horizon/dt)+1))), ...
                               't_valid', vtime(pc, dt, cfg.tv_thr));
        for m = 1:2
            if m == 1, rom = rom_set_Re(G, Re); else, rom = rom_set_Re(best, Re); end
            curves = nan(numel(starts), Hn + 1);  cc = curves;
            for s = 1:numel(starts)
                k = starts(s) + (0:Hn);
                Ah = rom_integrate(rom, A(:, k(1)), t(k), amax, cfg.tol);
                curves(s, :) = rom_error(Ah, A(:, k), d.nw2(k));
                cc(s, :) = sqrt(sum((Ah - A(:, k)).^2, 1)) / arms;
            end
            hk = 1:round(cfg.horizon/dt) + 1;
            [Ah, st] = rom_integrate(rom, A(:, 1), t, amax, 1e-6);
            [~, Dr] = rom_ED(Ah, pre, Re);
            res = struct('r', r, 'Re', Re, 'seen', q <= nt, 'method', rom.type, ...
                'err_short', mean(rowmean(curves(:, hk))), 'coef_short', mean(rowmean(cc(:, hk))), ...
                't_valid', vtime(cc, dt, cfg.tv_thr), 'free_status', st, 'D_free', Dr / Dlam, ...
                'stats', burst_stats(t, Dr / Dlam), 'fom_stats', burst_stats(t, Dref), ...
                'D_fom', Dref, 't', t, 'proj_err', mean(proj_err(A, d.nw2)));
            say(cfg, '   Re = %2g (%-6s) %-12s field err %.3f coef err %.3f t_v %5.1f free %-6s <D> %.3f (FOM %.3f) std %.3f (FOM %.3f)\n', ...
                Re, ternary(q <= nt, 'train', 'UNSEEN'), rom.type, res.err_short, res.coef_short, ...
                res.t_valid, st, res.stats.mean, res.fom_stats.mean, res.stats.std, res.fom_stats.std);
            P.res{ir, q, m} = res;
        end
        say(cfg, '   Re = %2g            persistence  coef err %.3f t_v %5.1f\n', Re, P.pers{ir, q}.coef_short, P.pers{ir, q}.t_valid);
    end
    P.roms{ir} = struct('G', G, 'O', best, 'infoO', infoO, 'tG', tG, 'tO', tO);
end
P.allRe = [cfg.train_Re, cfg.test_Re];
save(cfg.outfile, 'P', '-v7.3');
say(cfg, 'saved %s\n', cfg.outfile);
end

% ---------------------------------------------------------------- helpers
function [A, nw2] = project_chunked(W, V, wbar)
m = size(W, 2);  A = zeros(size(V, 2), m);  nw2 = zeros(1, m);
for j0 = 1:500:m
    j = j0:min(m, j0+499);
    X = double(W(:, j)) - wbar;
    A(:, j) = V' * X;   nw2(j) = sum(X.^2, 1);
end
end
function e = proj_err(A, nw2)
p2 = nw2 - sum(A.^2, 1);  p2(p2 < 0) = 0;  e = sqrt(p2 ./ nw2);
end
function m = rowmean(X)
X(isnan(X)) = 2;  m = mean(X, 2);
end
function tv = vtime(curves, dt, thr)
tv = zeros(size(curves, 1), 1);
for s = 1:size(curves, 1)
    k = find(~(curves(s, :) <= thr), 1);  if isempty(k), k = size(curves, 2); end
    tv(s) = (k - 1) * dt;
end
tv = mean(tv);
end
function o = ternary(c, a, b)
if c, o = a; else, o = b; end
end

function say(cfg, varargin)
% print to the console and append to cfg.logfile (live progress log)
s = sprintf(varargin{:});
fprintf('%s', s);
if isfield(cfg, 'logfile') && ~isempty(cfg.logfile)
    fid = fopen(cfg.logfile, 'a');  fprintf(fid, '%s', s);  fclose(fid);
end
end
