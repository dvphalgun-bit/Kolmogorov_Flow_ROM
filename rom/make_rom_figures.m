%% MAKE_ROM_FIGURES  Figures and summary tables of the ROM study.
% Reads ../rom_results/rom_Re20.mat, rom_Re30.mat, rom_Re40.mat, diag_Re40.mat,
% rom_param.mat (and baseline_standard_POD/ for the symmetry comparison) and
% writes PNG figures to ../rom_results/figures/ and a text summary.

clear; close all;
here = fileparts(mfilename('fullpath'));  root = fileparts(here);  addpath(here);
res_dir = fullfile(root, 'rom_results');  fig_dir = fullfile(res_dir, 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
names = {'Re20', 'Re30', 'Re40'};  Res = [20 30 40];
for c = 1:3
    L = load(fullfile(res_dir, ['rom_' names{c} '.mat']), 'R');  RR{c} = L.R;  %#ok<SAGROW>
end
cG = [0.85 0.33 0.10];  cO = [0.00 0.45 0.74];  cF = [0.2 0.2 0.2];  cP = [0.5 0.5 0.5];
sumf = fopen(fullfile(res_dir, 'summary_tables.txt'), 'w');
P = @(varargin) fprintf(sumf, varargin{:});

%% ---------------------------------------------------------------- tables
for c = 1:3
    R = RR{c};
    P('\n==== %s  (train t in [%g,%g], test t in (%g,%g]) ====\n', names{c}, R.cfg.t_train, R.cfg.t_test);
    P(' r  | proj err train/test | coef err  Gal / OpInf / persist. | t_v Gal / OpInf / persist. | field err Gal / OpInf | free run <D>/D_lam Gal / OpInf (FOM %.3f) | std Gal / OpInf (FOM %.3f) | OpInf lam, cond, fit resid\n', ...
        R.fom_stats.mean, R.fom_stats.std);
    for ir = 1:numel(R.r_list)
        g = R.res{ir, 1};  o = R.res{ir, 2};
        if o.built
            os = sprintf('%.3f', o.coef_short);  ot = sprintf('%5.1f', o.t_valid);  of = sprintf('%.3f', o.err_short);
            od = sprintf('%.3f%s', o.free_stats.mean, flag(o));  osd = sprintf('%.3f', o.free_stats.std);
            oi = sprintf('[%g %g], %.1e, %.3f', o.info.lam, o.info.cond, o.info.resid);
        else
            [os, ot, of, od, osd, oi] = deal('  -  ');
        end
        pe = R.pers{ir};
        P(' %2d |   %.3f / %.3f     |  %.3f / %s / %.3f  |  %5.1f / %s / %5.1f  |   %.3f / %s     |   %.3f%s / %s   |  %.3f / %s  | %s\n', R.r_list(ir), ...
            R.proj_err_train(ir), R.proj_err_test(ir), g.coef_short, os, pe.coef_short, g.t_valid, ot, pe.t_valid, g.err_short, of, ...
            g.free_stats.mean, flag(g), od, g.free_stats.std, osd, oi);
    end
end

%% ---------------------------------------------------------------- Fig 1: symmetry
f = figure('Position', [50 50 1100 380], 'Color', 'w');
for c = [1 3]
    subplot(1, 2, 1 + (c == 3));
    R = RR{c};  hold on;
    bfile = fullfile(res_dir, 'baseline_standard_POD', ['rom_' names{c} '.mat']);
    if exist(bfile, 'file')
        B = load(bfile, 'R');  rb = cellfun(@(x) x.r, B.R.res(:, 1));
        plot(rb, B.R.proj_err_test, 's--', 'Color', cG, 'LineWidth', 1.3, 'DisplayName', 'standard POD, test window');
    end
    plot(R.r_list, R.proj_err_train, 'o:', 'Color', cO, 'LineWidth', 1.3, 'DisplayName', 'invariant POD, training window');
    plot(R.r_list, R.proj_err_test, 'o-', 'Color', cO, 'LineWidth', 1.5, 'DisplayName', 'invariant POD, test window');
    grid on;  xlabel('number of modes r');  ylabel('relative projection error');
    title(sprintf('Re = %d: POD projection error', Res(c)));  legend('Location', 'northeast');  ylim([0 1]);
end
save_fig(f, fig_dir, 'fig_r1_symmetry_projection');

%% ---------------------------------------------------------------- Fig 2: accuracy vs r
f = figure('Position', [50 50 1300 700], 'Color', 'w');
for c = 1:3
    R = RR{c};  [cg, co, tg, to, fg, fo, tp, fp] = deal(nan(1, numel(R.r_list)));
    for ir = 1:numel(R.r_list)
        g = R.res{ir, 1};  o = R.res{ir, 2};
        cg(ir) = g.coef_short;  tg(ir) = g.t_valid;  fg(ir) = g.err_short;
        tp(ir) = R.pers{ir}.t_valid;  fp(ir) = R.pers{ir}.err_short;
        if o.built, co(ir) = o.coef_short; to(ir) = o.t_valid; fo(ir) = o.err_short; end
    end
    subplot(2, 3, c);  hold on;
    plot(R.r_list, fg, 's-', 'Color', cG, 'LineWidth', 1.5, 'DisplayName', 'POD-Galerkin');
    plot(R.r_list, fo, 'o-', 'Color', cO, 'LineWidth', 1.5, 'DisplayName', 'OpInf');
    plot(R.r_list, fp, ':', 'Color', cP, 'LineWidth', 1.5, 'DisplayName', 'persistence a(t)=a(t_0)');
    plot(R.r_list, R.proj_err_test, 'k--', 'LineWidth', 1, 'DisplayName', 'projection (best possible)');
    grid on;  xlabel('r');  ylabel('field error');  set(gca, 'YScale', 'log');
    title(sprintf('Re = %d: field error, first %g time units', Res(c), R.cfg.horizon));
    if c == 1, legend('Location', 'southwest'); end
    subplot(2, 3, 3 + c);  hold on;
    plot(R.r_list, tg, 's-', 'Color', cG, 'LineWidth', 1.5);
    plot(R.r_list, to, 'o-', 'Color', cO, 'LineWidth', 1.5);
    plot(R.r_list, tp, ':', 'Color', cP, 'LineWidth', 1.5);
    yline(R.cfg.Hv, 'k:', 'test horizon');
    grid on;  xlabel('r');  ylabel('valid prediction time');
    title(sprintf('Re = %d: time until coef. error > %.1f', Res(c), R.cfg.tv_thr));
end
save_fig(f, fig_dir, 'fig_r2_accuracy_vs_r');

%% ---------------------------------------------------------------- Fig 3: error vs time
f = figure('Position', [50 50 1300 380], 'Color', 'w');
rshow = [20 20 30];
for c = 1:3
    R = RR{c};  ir = find(R.r_list >= rshow(c), 1);
    dt = R.t(2) - R.t(1);  subplot(1, 3, c);  hold on;
    for im = 1:2
        x = R.res{ir, im};  if ~x.built, continue; end
        cc = x.coef_curves;  cc(isnan(cc)) = 2;  tt = (0:size(cc, 2)-1) * dt;
        col = cG;  if im == 2, col = cO; end
        plot(tt, cc', '-', 'Color', [col 0.25], 'HandleVisibility', 'off');
        plot(tt, mean(cc, 1), '-', 'Color', col, 'LineWidth', 2, 'DisplayName', x.rom.type);
    end
    pcv = R.pers{ir}.coef_curves;
    plot((0:size(pcv, 2)-1) * dt, mean(pcv, 1), ':', 'Color', cP, 'LineWidth', 2, 'DisplayName', 'persistence');
    yline(R.cfg.tv_thr, 'k:', 'HandleVisibility', 'off');
    grid on;  xlabel('time since start of prediction');  ylabel('coefficient error');
    title(sprintf('Re = %d, r = %d (test window, %d starts)', Res(c), R.r_list(ir), size(cc, 1)));
    ylim([0 2]);  legend('Location', 'southeast');
end
save_fig(f, fig_dir, 'fig_r3_error_vs_time');

%% ---------------------------------------------------------------- Fig 4: long free runs at Re = 40
R = RR{3};  kfree = find(R.t >= R.cfg.t_train(1));  tf = R.t(kfree);
rr = [10 30 121];                                            % 121: Galerkin only (converged)
f = figure('Position', [50 50 1300 800], 'Color', 'w');
for q = 1:numel(rr)
    ir = find(R.r_list >= rr(q), 1);
    subplot(numel(rr), 4, (q-1)*4 + (1:3));  hold on;
    plot(tf, R.D_fom(kfree), '-', 'Color', cF, 'DisplayName', 'FOM');
    for im = 1:2
        x = R.res{ir, im};  if ~x.built, continue; end
        col = cG;  if im == 2, col = cO; end
        plot(tf, x.D_free, '-', 'Color', col, 'DisplayName', sprintf('%s (%s)', x.rom.type, x.free_status));
    end
    xline(R.cfg.t_train(2), 'k--', 'end of training', 'HandleVisibility', 'off');
    grid on;  ylim([0 0.4]);  ylabel('D / D_{lam}');
    title(sprintf('Re = 40, r = %d: free run from t = %g', R.r_list(ir), tf(1)));
    if q == 1, legend('Location', 'northwest'); end
    if q == numel(rr), xlabel('time t'); end
    subplot(numel(rr), 4, q*4);  hold on;
    edges = linspace(0, 0.4, 41);
    histogram(R.D_fom(kfree), edges, 'Normalization', 'pdf', 'FaceColor', cF, 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'DisplayName', 'FOM');
    for im = 1:2
        x = R.res{ir, im};  if ~x.built, continue; end
        col = cG;  if im == 2, col = cO; end
        d = x.D_free(~isnan(x.D_free));
        histogram(d, edges, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', col, 'LineWidth', 1.5);
    end
    xline(0.15, 'k:');  grid on;  xlabel('D / D_{lam}');  ylabel('pdf');  title('distribution');
    set(gca, 'YScale', 'log');
end
save_fig(f, fig_dir, 'fig_r4_free_run_Re40');

%% ---------------------------------------------------------------- Fig 5: statistics vs r
f = figure('Position', [50 50 1300 380], 'Color', 'w');
for c = 1:3
    R = RR{c};  subplot(1, 3, c);  hold on;
    [mg, mo, sg, so] = deal(nan(1, numel(R.r_list)));
    for ir = 1:numel(R.r_list)
        g = R.res{ir, 1};  o = R.res{ir, 2};
        if strcmp(g.free_status, 'ok'), mg(ir) = g.free_stats.mean; sg(ir) = g.free_stats.std; end
        if o.built && strcmp(o.free_status, 'ok'), mo(ir) = o.free_stats.mean; so(ir) = o.free_stats.std; end
    end
    errorbar(R.r_list - 0.3, mg, sg, 's', 'Color', cG, 'LineWidth', 1.3, 'DisplayName', 'POD-Galerkin');
    errorbar(R.r_list + 0.3, mo, so, 'o', 'Color', cO, 'LineWidth', 1.3, 'DisplayName', 'OpInf');
    yline(R.fom_stats.mean, '-', 'Color', cF, 'DisplayName', 'FOM mean');
    yline(R.fom_stats.mean + R.fom_stats.std, ':', 'Color', cF, 'HandleVisibility', 'off');
    yline(R.fom_stats.mean - R.fom_stats.std, ':', 'Color', cF, 'HandleVisibility', 'off');
    grid on;  xlabel('r');  ylabel('<D>/D_{lam} \pm std (free run)');
    title(sprintf('Re = %d (missing marker = blow-up)', Res(c)));
    if c == 1, legend('Location', 'best'); end
end
save_fig(f, fig_dir, 'fig_r5_statistics_vs_r');

%% ---------------------------------------------------------------- Fig 6: energy budget
dfile = fullfile(res_dir, 'diag_Re40.mat');
if exist(dfile, 'file')
    L = load(dfile, 'G');  G = L.G;
    f = figure('Position', [50 50 1200 420], 'Color', 'w');
    subplot(1, 2, 1);  hold on;
    plot(G.cfg.r_list, G.proj_quiet, 'o-', 'Color', cO, 'LineWidth', 1.5, 'DisplayName', 'quiet snapshots');
    plot(G.cfg.r_list, G.proj_burst, 's-', 'Color', cG, 'LineWidth', 1.5, 'DisplayName', 'burst snapshots (D > 0.15 D_{lam})');
    grid on;  xlabel('r');  ylabel('projection error');  legend;  title('Re = 40: is the basis worse during bursts?');
    subplot(1, 2, 2);
    cl = G.closure;  cl = cl(arrayfun(@(d) ~isempty(d.r), cl));
    M = [[cl.Pt_q]; [cl.PG_q]; [cl.PO_q]; [cl.Pt_b]; [cl.PG_b]; [cl.PO_b]]';
    b = bar(categorical(arrayfun(@(d) sprintf('r = %d', d.r), cl, 'UniformOutput', false)), M);
    cols = [cF; cG; cO; cF; cG; cO];
    for k = 1:6, b(k).FaceColor = cols(k, :); if k > 3, b(k).FaceAlpha = 0.45; end, end
    set(gca, 'YScale', 'linear');  grid on;  ylabel('d(|a|^2/2)/dt');
    legend({'true (quiet)', 'Galerkin (quiet)', 'OpInf (quiet)', 'true (burst)', 'Galerkin (burst)', 'OpInf (burst)'}, 'Location', 'northwest');
    title('resolved energy budget on training snapshots');
    save_fig(f, fig_dir, 'fig_r6_energy_budget_Re40');
end

%% ---------------------------------------------------------------- Fig 7: unseen Re
pfile = fullfile(res_dir, 'rom_param.mat');
if exist(pfile, 'file')
    L = load(pfile, 'P');  PP = L.P;
    f = figure('Position', [50 50 1300 420], 'Color', 'w');
    ir = numel(PP.r_list);
    [~, o] = sort(PP.allRe);
    subplot(1, 2, 1);  hold on;
    for q = o
        x = PP.res{ir, q, 1};  y = PP.res{ir, q, 2};
        mk = 'o';  if ~x.seen, mk = 'p'; end
        errorbar(x.Re, x.fom_stats.mean, x.fom_stats.std, mk, 'Color', cF, 'MarkerFaceColor', cF, 'LineWidth', 1.2, 'HandleVisibility', 'off');
        if strcmp(x.free_status, 'ok'), errorbar(x.Re - 0.5, x.stats.mean, x.stats.std, mk, 'Color', cG, 'LineWidth', 1.2, 'HandleVisibility', 'off'); end
        if strcmp(y.free_status, 'ok'), errorbar(y.Re + 0.5, y.stats.mean, y.stats.std, mk, 'Color', cO, 'LineWidth', 1.2, 'HandleVisibility', 'off'); end
    end
    plot(nan, nan, 'o', 'Color', cF, 'MarkerFaceColor', cF, 'DisplayName', 'FOM');
    plot(nan, nan, 'o', 'Color', cG, 'DisplayName', 'POD-Galerkin');
    plot(nan, nan, 'o', 'Color', cO, 'DisplayName', 'OpInf');
    plot(nan, nan, 'pk', 'DisplayName', 'star = unseen Re');
    grid on;  xlabel('Re');  ylabel('<D>/D_{lam} \pm std (free run)');  legend('Location', 'northeast');
    title(sprintf('one parametric ROM, r = %d', PP.r_list(ir)));
    subplot(1, 2, 2);  hold on;
    for m = 1:2
        e = arrayfun(@(q) PP.res{ir, q, m}.coef_short, o);
        col = cG;  if m == 2, col = cO; end
        plot(PP.allRe(o), e, '-', 'Color', col, 'LineWidth', 1.5, 'DisplayName', PP.res{ir, 1, m}.method);
        un = ~arrayfun(@(q) PP.res{ir, q, m}.seen, o);
        plot(PP.allRe(o(un)), e(un), 'p', 'Color', col, 'MarkerSize', 12, 'MarkerFaceColor', col, 'HandleVisibility', 'off');
    end
    grid on;  xlabel('Re');  ylabel(sprintf('coef. error, first %g time units', PP.cfg.horizon));
    legend('Location', 'northwest');  title('short-horizon accuracy (stars: unseen Re)');
    save_fig(f, fig_dir, 'fig_r7_unseen_Re');
    P('\n==== parametric ROM (train Re = %s) ====\n', mat2str(PP.cfg.train_Re));
    for ir2 = 1:numel(PP.r_list)
        for q = o
            for m = 1:2
                x = PP.res{ir2, q, m};
                P(' r = %2d  Re = %2d (%-6s) %-12s coef err %.3f  t_v %5.1f  free %-6s <D> %.3f (FOM %.3f)  std %.3f (FOM %.3f)\n', ...
                    PP.r_list(ir2), x.Re, tern(x.seen, 'train', 'UNSEEN'), x.method, x.coef_short, x.t_valid, ...
                    x.free_status, x.stats.mean, x.fom_stats.mean, x.stats.std, x.fom_stats.std);
            end
        end
    end
end

%% ---------------------------------------------------------------- Fig 8: vorticity fields at Re = 40
R = RR{3};  ir = find(R.r_list >= 30, 1);  N = R.cfg.N;
S = matfile(R.cfg.train_file);
it = find(R.test);  k0 = it(1);  dt = R.t(2) - R.t(1);
lags = [0 2 5 10];
f = figure('Position', [50 50 1300 950], 'Color', 'w');
At = R.Aall_test(1:R.r_list(ir), :);                        % projected FOM in the test window
kk = k0 + round(lags / dt);                                  % snapshot indices of the lags
Ahs = cell(1, 2);
for im = 1:2                                                 % one integration per ROM
    x = R.res{ir, im};
    if x.built, Ahs{im} = rom_integrate(x.rom, At(:, 1), R.t(kk), 1e3 * norm(At(:, 1)), 1e-8); end
end
for q = 1:numel(lags)
    k = kk(q);
    wF = reshape(double(S.W(:, k)), N, N);
    wP = reshape(R.wbar + R.V(:, 1:R.r_list(ir)) * At(:, k - k0 + 1), N, N);
    fields = {wF, wP};  ttl = {'FOM', sprintf('projection, r = %d', R.r_list(ir))};
    for im = 1:2
        x = R.res{ir, im};
        if ~x.built, fields{end+1} = nan(N); ttl{end+1} = 'n/a'; continue; end   %#ok<SAGROW>
        fields{end+1} = reshape(R.wbar + R.V(:, 1:R.r_list(ir)) * Ahs{im}(:, q), N, N);  %#ok<SAGROW>
        ttl{end+1} = x.rom.type;                                                         %#ok<SAGROW>
    end
    lim = max(abs(wF(:)));
    for j = 1:4
        subplot(numel(lags), 4, (q-1)*4 + j);
        imagesc(fields{j});  axis xy equal tight off;
        caxis([-lim lim]);  colormap(gca, bwr(256));
        title(sprintf('%s, t_0 + %g', ttl{j}, lags(q)));
    end
end
save_fig(f, fig_dir, 'fig_r8_fields_Re40');

%% ---------------------------------------------------------------- Fig 9: lambda search
f = figure('Position', [50 50 1300 380], 'Color', 'w');
rlam = [10 10 20];
for c = 1:3
    R = RR{c};  ir = find(R.r_list >= rlam(c), 1);  o = R.res{ir, 2};
    subplot(1, 3, c);
    if ~o.built || ~isfield(o, 'lam_table'), axis off; continue; end
    T = o.lam_table;  S = T.score;  bl = isinf(S);  S(bl) = NaN;
    imagesc(log10(T.lam2s), log10(T.lam1s), log10(S));  axis xy;  colorbar;
    hold on;  [i1, i2] = find(bl);
    plot(log10(T.lam2s(i2)), log10(T.lam1s(i1)), 'kx', 'MarkerSize', 8, 'LineWidth', 1.2);
    plot(log10(o.info.lam(2)), log10(o.info.lam(1)), 'wp', 'MarkerSize', 14, 'MarkerFaceColor', 'w');
    xlabel('log_{10} \lambda_2 (quadratic)');  ylabel('log_{10} \lambda_1 (linear)');
    title(sprintf('Re = %d, r = %d: log_{10} score (x = blow-up)', Res(c), R.r_list(ir)));
    P('\n%s r = %d lambda search: %d of %d models blew up; chosen lam = [%g %g]\n', names{c}, ...
        R.r_list(ir), nnz(bl), numel(bl), o.info.lam);
    P('   fit residual range %.3f - %.3f; condition number range %.1e - %.1e\n', ...
        min(T.resid(:)), max(T.resid(:)), min(T.cond(:)), max(T.cond(:)));
end
save_fig(f, fig_dir, 'fig_r9_lambda_search');

fclose(sumf);
fprintf('figures written to %s\n', fig_dir);
type(fullfile(res_dir, 'summary_tables.txt'));

%% ---------------------------------------------------------------- helpers
function s = flag(x)
s = '';  if ~strcmp(x.free_status, 'ok'), s = '*'; end   % * = free run blew up
end
function o = tern(c, a, b)
if c, o = a; else, o = b; end
end
function save_fig(f, d, name)
try, theme(f, 'light'); catch, end
set(f, 'Color', 'w');
exportgraphics(f, fullfile(d, [name '.png']), 'Resolution', 200);
end
function c = bwr(m)
h = floor(m/2);  up = linspace(0, 1, h)';
c = [[up up ones(h,1)]; [ones(m-h,1) flipud(linspace(0,1,m-h)') flipud(linspace(0,1,m-h)')]];
end
