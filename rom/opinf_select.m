function [best, info, table] = opinf_select(Atr, Ad, Aproj, nw2, t, train, lam1s, lam2s, opts)
% OPINF_SELECT  Choose the Operator Inference regularisation by grid search,
% using ONLY training data.
%
% Score of a (lam1, lam2) pair = e_short + e_stat, where
%   e_short = mean short-horizon COEFFICIENT error from opts.nstarts start
%             times inside the training window (dynamics accuracy), and
%   e_stat  = |log(<|a|^2>_free / <|a|^2>_train)| + |log(var_free / var_train)|:
%             mismatch of the mean energy and of the temporal variance of a free
%             run started from the first training state, both evaluated over the
%             SAME time window as the training data (so transients are compared
%             like with like). The free run is continued to length opts.Tfree
%             only to check that it stays bounded.
% A model whose free run blows up is rejected. The e_stat term rejects energy
% growth AND two degenerate over-regularised models that e_short alone would
% favour in chaotic flow: (i) dynamics decaying to the mean (error ~1, less than
% the ~sqrt(2) of realistic but decorrelated dynamics) and (ii) a "frozen"
% model with da/dt ~ 0, which keeps a typical energy but has no variance and
% scores well whenever the flow changes slowly over the short horizon.
%
%   Atr, Ad        fit data (r x K states / derivatives, may include shifted copies)
%   Aproj, nw2, t  projected FOM coefficients, ||w'||^2, times (full record)
%   train          logical index of the training window in t
%   opts: horizon, nstarts, amax, Tfree, arms, e2_train, var_train;
%         optional Reinv / Re (parametric fit evaluated at one Re)
if ~isfield(opts, 'Reinv'), opts.Reinv = []; end
it = find(train);  dt = t(2) - t(1);
starts = round(linspace(it(1), it(end) - round(opts.horizon/dt), opts.nstarts));
tf = t(it(1)) + (0:1:opts.Tfree);                    % free run, output every time unit
ttr = t(it(end));                                    % end of the training window
table = struct('score', nan(numel(lam1s), numel(lam2s)));
[table.e_short, table.e_stat, table.resid, table.cond] = deal(table.score);
table.lam1s = lam1s;  table.lam2s = lam2s;
best = [];  info = [];  sbest = Inf;
for i = 1:numel(lam1s)
    for j = 1:numel(lam2s)
        [rom, inf_ij] = opinf_fit(Atr, Ad, lam1s(i), lam2s(j), opts.Reinv);
        if isfield(opts, 'Re') && ~isempty(opts.Re), rom = rom_set_Re(rom, opts.Re); end
        table.resid(i, j) = inf_ij.resid;  table.cond(i, j) = inf_ij.cond;
        [Af, st] = rom_integrate(rom, Aproj(:, it(1)), tf, opts.amax, 1e-6);
        if ~strcmp(st, 'ok'), table.score(i, j) = Inf; continue; end
        Aw = Af(:, tf > tf(1) & tf <= ttr);              % training window only
        e_stat = abs(log(mean(sum(Aw.^2, 1)) / opts.e2_train)) + ...
                 abs(log(tvar(Aw) / opts.var_train));
        [~, ok, ~, e_short] = rom_short_error(rom, Aproj, nw2, t, starts, opts.horizon, opts.amax, opts.arms);
        table.e_short(i, j) = e_short;  table.e_stat(i, j) = e_stat;
        if ~ok, table.score(i, j) = Inf; continue; end
        sc = e_short + e_stat;
        table.score(i, j) = sc;
        if sc < sbest
            sbest = sc;  best = rom;  info = inf_ij;
            info.lam = [lam1s(i) lam2s(j)];  info.score = sc;
            info.e_short = e_short;  info.e_stat = e_stat;
        end
    end
end
if isempty(best), warning('opinf_select: no bounded model on the grid'); end
end

function v = tvar(A)
% total temporal variance  <|a - <a>_t|^2>_t
v = mean(sum((A - mean(A, 2)).^2, 1));
end
