function R = run_rom_case(cfg)
% RUN_ROM_CASE  Steps 3-6 for one Reynolds number: POD basis, POD-Galerkin and
% Operator Inference ROMs for every r in cfg.r_list, and their evaluation
% against the same full-order reference.
%
% cfg fields
%   name, Re, N, n           case label and FOM parameters
%   train_file               fine-sampled FOM data (W, t_snap, diagn)
%   t_train, t_test          [t0 t1] training window / future-time test window
%   ic_file                  FOM run from a different initial condition ('' = none)
%   r_list                   reduced dimensions
%   ratio_min                OpInf only where samples/unknowns >= ratio_min
%   stride                   snapshot stride for the SVD
%   lam1s, lam2s             OpInf regularisation grid
%   horizon                  short-horizon window for errors / lambda selection
%   Hv                       horizon for the valid prediction time
%   nstarts_sel, nstarts_test number of start times (selection / test)
%   tv_thr                   error threshold defining the valid prediction time
%   tol                      ode45 RelTol for evaluation runs
%   basis                    'invariant' (x-shift symmetric POD, default) or 'standard'
%   n_aug, rows_max          OpInf training data augmented with n_aug random
%                            x-shifted copies (exact trajectories by symmetry),
%                            subsampled to at most rows_max rows
%   outfile                  results .mat file

ops  = kf_ops(cfg.N, cfg.n);
Dlam = cfg.Re / (2*cfg.n^2);
say(cfg, '\n===== %s (Re = %g) =====\n', cfg.name, cfg.Re);

% ---------------------------------------------------------------- data
S  = load(cfg.train_file, 'W', 't_snap', 'diagn');
t  = S.t_snap;   dt = t(2) - t(1);
train = t >= cfg.t_train(1) & t <= cfg.t_train(2);
test  = t >  cfg.t_test(1)  & t <= cfg.t_test(2);
D_fom = interp1(S.diagn.t_diag, S.diagn.D, t) / Dlam;
E_fom = interp1(S.diagn.t_diag, S.diagn.E, t);

% ---------------------------------------------------------------- Step 3: POD
if ~isfield(cfg, 'basis'), cfg.basis = 'invariant'; end
if ~isfield(cfg, 'n_aug'), cfg.n_aug = 0; end
t0 = tic;
rmax = max(cfg.r_list);
if strcmp(cfg.basis, 'invariant')
    [V, wbar, ev, binfo] = pod_basis_inv(S.W(:, train), rmax, cfg.stride, cfg.N);
    sv = sqrt(ev);
    r_list = arrayfun(@(r) min(binfo.cuts(binfo.cuts >= r)), cfg.r_list);  % keep pairs whole
else
    [V, wbar, sv] = pod_basis(S.W(:, train), rmax, cfg.stride);
    binfo = [];  r_list = cfg.r_list;
end
t_pod = toc(t0);
[Aall, nw2] = project_chunked(S.W, V, wbar);
clear S
pre = rom_ED_pre(V, wbar, ops);
say(cfg, 'POD (%s): %d training snapshots (stride %d), %.1f s; r = %s\n', ...
    cfg.basis, nnz(train), cfg.stride, t_pod, mat2str(r_list));

% derivatives of the training coefficients (for OpInf)
[Ad_all, keep] = fd_derivative(Aall(:, train), dt);
Atr_all = Aall(:, train);  Atr_all = Atr_all(:, keep);
K = size(Atr_all, 2);

% new-initial-condition reference
hasIC = isfield(cfg, 'ic_file') && ~isempty(cfg.ic_file);
if hasIC
    Sic = load(cfg.ic_file, 'W', 't_snap', 'diagn');
    tic_ = Sic.t_snap;
    [Aic, nw2ic] = project_chunked(Sic.W, V, wbar);
    D_ic = interp1(Sic.diagn.t_diag, Sic.diagn.D, tic_) / Dlam;
    clear Sic
end

it  = find(test);
Hn  = round(cfg.Hv / dt);
starts_test = round(linspace(it(1), it(end) - Hn, cfg.nstarts_test));
tfree = t(train);  tfree = [tfree(1), t(end)];          % free run over the whole record
kfree = find(t >= tfree(1));

R = struct('cfg', cfg, 'sv', sv, 't', t, 'D_fom', D_fom, 'E_fom', E_fom, ...
           't_pod', t_pod, 'V', V, 'wbar', wbar, 'train', train, 'test', test);
R.binfo = binfo;  R.r_list = r_list;
R.proj_err_test = nan(1, numel(r_list));  R.proj_err_train = R.proj_err_test;
rng(0);  thetas = 2*pi*rand(1, cfg.n_aug);              % fixed random x-shifts
R.fom_stats = burst_stats(t(kfree), D_fom(kfree));
if hasIC, R.t_ic = tic_; R.D_ic = D_ic; R.ic_stats = burst_stats(tic_, D_ic); end

methods = {'Galerkin', 'OpInf'};
for ir = 1:numel(r_list)
    r = r_list(ir);
    A  = Aall(1:r, :);
    [~, proj] = rom_error(zeros(r, numel(t)), A, nw2);
    R.proj_err_test(ir) = mean(proj(test));  R.proj_err_train(ir) = mean(proj(train));
    amax = 10 * max(sqrt(sum(Atr_all(1:r, :).^2, 1)));
    arms = sqrt(mean(sum(Atr_all(1:r, :).^2, 1)));
    p = 1 + r + r*(r+1)/2;
    % persistence baseline: a(t) = a(t0)  ("the flow does not change")
    pc = nan(numel(starts_test), Hn + 1);  pf = pc;
    for s = 1:numel(starts_test)
        k = starts_test(s) + (0:Hn);
        Ah = repmat(A(:, k(1)), 1, numel(k));
        pc(s, :) = sqrt(sum((Ah - A(:, k)).^2, 1)) / arms;
        pf(s, :) = rom_error(Ah, A(:, k), nw2(k));
    end
    hk = 1:round(cfg.horizon/dt) + 1;
    R.pers{ir} = struct('coef_curves', pc, 'coef_short', mean(mean(pc(:, hk), 2)), ...
        'err_short', mean(mean(pf(:, hk), 2)), 't_valid', valid_time(pc, dt, cfg.tv_thr));
    for im = 1:2
        res = struct('r', r, 'method', methods{im}, 'built', false);
        t0 = tic;
        if im == 1
            rom = galerkin_rom(V(:, 1:r), wbar, ops, cfg.Re);
            res.built = true;  res.info = struct();
        elseif K / p >= cfg.ratio_min
            % fit data: original + x-shifted copies (symmetry augmentation)
            [Xa, Xd] = augment(Atr_all(1:r, :), Ad_all(1:r, :), binfo, thetas, cfg);
            Ar = Atr_all(1:r, :);
            opts = struct('horizon', cfg.horizon, 'nstarts', cfg.nstarts_sel, ...
                          'amax', amax, 'Tfree', tfree(2) - tfree(1), 'arms', arms, ...
                          'e2_train', mean(sum(Ar.^2, 1)), ...
                          'var_train', mean(sum((Ar - mean(Ar, 2)).^2, 1)));
            [rom, info, tab] = opinf_select(Xa, Xd, A, nw2, t, train, cfg.lam1s, cfg.lam2s, opts);
            res.built = ~isempty(rom);  res.info = info;  res.lam_table = tab;
            if res.built, res.info.rows = size(Xa, 2); res.info.ratio_orig = K / p; end
        else
            rom = [];  res.info = struct('ratio', K/p);
        end
        res.t_offline = toc(t0);
        if ~res.built
            say(cfg, '  r = %2d  %-8s  not built (K/p = %.1f)\n', r, methods{im}, K/p);
            R.res{ir, im} = res;  continue;
        end
        res.rom = rom;
        % ---- (a) future-time test: short-horizon error and valid prediction time
        % field error  = vs the full FOM field (includes the basis/projection error)
        % coef. error  = vs the projected FOM coefficients, normalised by the RMS
        %                coefficient norm: isolates the error of the ROM dynamics
        curves = nan(numel(starts_test), Hn + 1);  ccurves = curves;
        for s = 1:numel(starts_test)
            k = starts_test(s) + (0:Hn);
            Ah = rom_integrate(rom, A(:, k(1)), t(k), amax, cfg.tol);
            curves(s, :)  = rom_error(Ah, A(:, k), nw2(k));
            ccurves(s, :) = sqrt(sum((Ah - A(:, k)).^2, 1)) / arms;
        end
        hk = 1:round(cfg.horizon/dt) + 1;
        res.err_curves  = curves;   res.coef_curves = ccurves;
        res.err_short   = mean(nanmean_rows(curves(:, hk)));
        res.coef_short  = mean(nanmean_rows(ccurves(:, hk)));
        res.t_valid     = valid_time(ccurves, dt, cfg.tv_thr);
        % ---- (b) long free run from the first training state over the whole record
        [Ah, st, tstop, cpu] = rom_integrate(rom, A(:, kfree(1)), t(kfree), amax, 1e-6);
        [E_r, D_r] = rom_ED(Ah, pre, cfg.Re);
        res.free_status = st;   res.free_tstop = tstop;
        res.cpu_per_time = cpu / (min(tstop, t(end)) - t(kfree(1)));
        res.D_free = D_r / Dlam;   res.E_free = E_r;
        res.free_stats = burst_stats(t(kfree), res.D_free);
        res.err_free = rom_error(Ah, A(:, kfree), nw2(kfree));
        % ---- (c) new initial condition
        if hasIC
            [Ah, st2] = rom_integrate(rom, Aic(1:r, 1), tic_, amax, 1e-6);
            [~, D_r] = rom_ED(Ah, pre, cfg.Re);
            res.ic_status = st2;   res.D_ic = D_r / Dlam;
            res.ic_stats  = burst_stats(tic_, res.D_ic);
            res.err_ic    = rom_error(Ah, Aic(1:r, :), nw2ic);
        end
        say(cfg, '  r = %2d  %-8s  field err %.3f  coef err %.3f  t_v %5.1f  free run %-6s  <D> %.3f (FOM %.3f)  offline %.1fs  online %.3fs/t\n', ...
            r, methods{im}, res.err_short, res.coef_short, res.t_valid, st, res.free_stats.mean, ...
            R.fom_stats.mean, res.t_offline, res.cpu_per_time);
        if im == 2
            say(cfg, '                  OpInf lam = [%g %g], cond %.2e, rows %d (K/p %.1f before augmentation), fit resid %.3f, score %.3f = %.3f + %.3f\n', ...
                res.info.lam, res.info.cond, res.info.rows, res.info.ratio_orig, res.info.resid, ...
                res.info.score, res.info.e_short, res.info.e_stat);
        end
        R.res{ir, im} = res;
    end
    say(cfg, '  r = %2d  persistence baseline: field err %.3f  coef err %.3f  t_v %5.1f | projection error (test) %.3f\n', ...
        r, R.pers{ir}.err_short, R.pers{ir}.coef_short, R.pers{ir}.t_valid, R.proj_err_test(ir));
end
R.Aall_test = Aall(:, test);  R.nw2 = nw2;
save(cfg.outfile, 'R', '-v7.3');
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

function [Xa, Xd] = augment(X, Xdot, binfo, thetas, cfg)
% original training data plus x-shifted copies, subsampled to <= cfg.rows_max rows
Xa = X;  Xd = Xdot;
if isempty(binfo) || isempty(thetas), return; end
for q = 1:numel(thetas)
    Xa = [Xa, shift_coeffs(X, binfo, thetas(q))];      %#ok<AGROW>
    Xd = [Xd, shift_coeffs(Xdot, binfo, thetas(q))];   %#ok<AGROW>
end
st = max(1, ceil(size(Xa, 2) / cfg.rows_max));
Xa = Xa(:, 1:st:end);  Xd = Xd(:, 1:st:end);
end

function m = nanmean_rows(X)
% row means; NaN entries (after a blow-up) count as error 2, i.e. worse than an
% uncorrelated prediction (sqrt(2)), so a blow-up is penalised, not ignored
X(isnan(X)) = 2;
m = mean(X, 2);
end

function tv = valid_time(curves, dt, thr)
% mean over start times of the first time the error exceeds thr
tv = zeros(size(curves, 1), 1);
for s = 1:size(curves, 1)
    k = find(~(curves(s, :) <= thr), 1);         % NaN (blow-up) also counts
    if isempty(k), k = size(curves, 2); end
    tv(s) = (k - 1) * dt;
end
tv = mean(tv);
end

function say(cfg, varargin)
% print to the console and append to cfg.logfile (live progress log)
s = sprintf(varargin{:});
fprintf('%s', s);
if isfield(cfg, 'logfile') && ~isempty(cfg.logfile)
    fid = fopen(cfg.logfile, 'a');  fprintf(fid, '%s', s);  fclose(fid);
end
end
