%% MAKE_ROM_FIGURES  Figures and summary tables of the ROM study.
% Reads ../rom_results/rom_Re20.mat, rom_Re30.mat, rom_Re40.mat, diag_Re40.mat,
% rom_param.mat (and baseline_standard_POD/ for the symmetry comparison) and
% writes PNG figures to ../rom_results/figures/ and a text summary.
%
% Every figure is drawn at the size it is printed (inches) with 8-9 pt text,
% so nothing has to be shrunk later:
%   paper, one column  (3.5 in):  fig_r1, fig_r7
%   paper, full width  (6.2 in):  fig_r2, fig_r4      (5.2 in: fig_r9)
%   slides only        (~7 in):   fig_r3, fig_r5, fig_r6, fig_r8

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

%% ---------------------------------------------------------------- Fig 1: symmetry (paper, one column)
[f, tl] = newfig(3.5, 2.4, 8, 1, 2);  ax0 = [];
for c = [1 3]
    ax = nexttile(tl);  hold on;  if isempty(ax0), ax0 = ax; end
    R = RR{c};
    bfile = fullfile(res_dir, 'baseline_standard_POD', ['rom_' names{c} '.mat']);
    if exist(bfile, 'file')
        B = load(bfile, 'R');  rb = cellfun(@(x) x.r, B.R.res(:, 1));
        plot(rb, B.R.proj_err_test, 's--', 'Color', cG, 'MarkerSize', 4, 'DisplayName', 'standard POD, test');
    end
    plot(R.r_list, R.proj_err_train, 'o:', 'Color', cO, 'MarkerSize', 4, 'DisplayName', 'invariant POD, train');
    plot(R.r_list, R.proj_err_test, 'o-', 'Color', cO, 'MarkerSize', 4, 'DisplayName', 'invariant POD, test');
    grid on;  box on;  xlabel('r');  ylim([0 1]);
    if c == 1, ylabel('projection error'); end
    title(sprintf('Re = %d', Res(c)));
end
shared_legend(ax0, 'north', 2);
save_fig(f, fig_dir, 'fig_r1_symmetry_projection');

%% ---------------------------------------------------------------- Fig 2: accuracy vs r (paper, full width)
[f, tl] = newfig(6.2, 3.9, 8, 2, 3);  axs = gobjects(2, 3);
for c = 1:3
    R = RR{c};  [tg, to, fg, fo, tp, fp] = deal(nan(1, numel(R.r_list)));
    for ir = 1:numel(R.r_list)
        g = R.res{ir, 1};  o = R.res{ir, 2};
        tg(ir) = g.t_valid;  fg(ir) = g.err_short;
        tp(ir) = R.pers{ir}.t_valid;  fp(ir) = R.pers{ir}.err_short;
        if o.built, to(ir) = o.t_valid; fo(ir) = o.err_short; end
    end
    axs(1, c) = nexttile(tl, c);  hold on;
    plot(R.r_list, fg, 's-', 'Color', cG, 'MarkerSize', 4, 'DisplayName', 'POD-Galerkin');
    plot(R.r_list, fo, 'o-', 'Color', cO, 'MarkerSize', 4, 'DisplayName', 'OpInf');
    plot(R.r_list, fp, ':', 'Color', cP, 'LineWidth', 1.5, 'DisplayName', 'persistence');
    plot(R.r_list, R.proj_err_test, 'k--', 'DisplayName', 'projection error (best possible)');
    grid on;  box on;  set(gca, 'YScale', 'log');
    yt = [0.05 0.1 0.2 0.3 0.5 1 2 3 5];  yl = ylim;  yt = yt(yt >= yl(1) & yt <= yl(2));
    set(gca, 'YTick', yt, 'YTickLabel', compose('%g', yt), 'YMinorGrid', 'off');
    if c == 1, ylabel('field error'); end
    title(sprintf('Re = %d', Res(c)));
    axs(2, c) = nexttile(tl, 3 + c);  hold on;
    plot(R.r_list, tg, 's-', 'Color', cG, 'MarkerSize', 4);
    plot(R.r_list, to, 'o-', 'Color', cO, 'MarkerSize', 4);
    plot(R.r_list, tp, ':', 'Color', cP, 'LineWidth', 1.5);
    yline(R.cfg.Hv, 'k:', 'test horizon', 'LabelVerticalAlignment', 'bottom', 'FontSize', 7);
    grid on;  box on;  xlabel('r');  ylim([0 1.15 * R.cfg.Hv]);
    if c == 1, ylabel('valid time t_v'); end
end
shared_legend(axs(1, 1), 'south', 4);
save_fig(f, fig_dir, 'fig_r2_accuracy_vs_r');

%% ---------------------------------------------------------------- Fig 3: error vs time (slides)
[f, tl] = newfig(7, 2.6, 9, 1, 3);  ax0 = [];
rshow = [20 20 30];
for c = 1:3
    R = RR{c};  ir = find(R.r_list >= rshow(c), 1);
    dt = R.t(2) - R.t(1);  ax = nexttile(tl);  hold on;  if isempty(ax0), ax0 = ax; end
    for im = 1:2
        x = R.res{ir, im};  if ~x.built, continue; end
        cc = x.coef_curves;  cc(isnan(cc)) = 2;  tt = (0:size(cc, 2)-1) * dt;
        col = cG;  if im == 2, col = cO; end
        plot(tt, cc', '-', 'Color', [col 0.25], 'LineWidth', 0.6, 'HandleVisibility', 'off');
        plot(tt, mean(cc, 1), '-', 'Color', col, 'LineWidth', 2, 'DisplayName', x.rom.type);
    end
    pcv = R.pers{ir}.coef_curves;
    plot((0:size(pcv, 2)-1) * dt, mean(pcv, 1), ':', 'Color', cP, 'LineWidth', 2, 'DisplayName', 'persistence');
    yline(R.cfg.tv_thr, 'k:', 'HandleVisibility', 'off');
    grid on;  box on;  xlabel('time since start');  ylim([0 2]);
    if c == 1, ylabel('coefficient error'); end
    title(sprintf('Re = %d, r = %d', Res(c), R.r_list(ir)));
end
shared_legend(ax0, 'south', 3);
save_fig(f, fig_dir, 'fig_r3_error_vs_time');

%% ---------------------------------------------------------------- Fig 4: long free runs at Re = 40 (paper, full width)
R = RR{3};  kfree = find(R.t >= R.cfg.t_train(1));  tf = R.t(kfree);
rr = [10 30 121];                                            % 121: Galerkin only (converged)
[f, tl] = newfig(6.2, 4.6, 8, numel(rr), 4);  ax0 = [];
edges = linspace(0, 0.4, 41);
for q = 1:numel(rr)
    ir = find(R.r_list >= rr(q), 1);
    ax = nexttile(tl, [1 3]);  hold on;  if isempty(ax0), ax0 = ax; end
    plot(tf, R.D_fom(kfree), '-', 'Color', cF, 'LineWidth', 0.8, 'DisplayName', 'FOM');
    for im = 1:2
        x = R.res{ir, im};  if ~x.built, continue; end
        col = cG;  if im == 2, col = cO; end
        plot(tf, x.D_free, '-', 'Color', col, 'LineWidth', 0.8, 'DisplayName', x.rom.type);
    end
    if q == 1
        xline(R.cfg.t_train(2), 'k--', 'end of training ', 'LabelOrientation', 'horizontal', ...
            'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top', 'FontSize', 7, 'HandleVisibility', 'off');
    else
        xline(R.cfg.t_train(2), 'k--', 'HandleVisibility', 'off');
    end
    grid on;  box on;  ylim([0 0.4]);  xlim([tf(1) tf(end)]);  ylabel('D / D_{lam}');
    title(sprintf('r = %d', R.r_list(ir)));
    if q == numel(rr), xlabel('time t'); end
    nexttile(tl);  hold on;
    histogram(R.D_fom(kfree), edges, 'Normalization', 'pdf', 'FaceColor', cF, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    for im = 1:2
        x = R.res{ir, im};  if ~x.built, continue; end
        col = cG;  if im == 2, col = cO; end
        d = x.D_free(~isnan(x.D_free));
        histogram(d, edges, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', col, 'LineWidth', 1.2);
    end
    xline(0.15, 'k:');  grid on;  box on;  ylabel('pdf');  set(gca, 'YScale', 'log');
    if q == 1, title('distribution'); end
    if q == numel(rr), xlabel('D / D_{lam}'); end
end
shared_legend(ax0, 'north', 3);
save_fig(f, fig_dir, 'fig_r4_free_run_Re40');

%% ---------------------------------------------------------------- Fig 5: statistics vs r (slides)
[f, tl] = newfig(7, 2.6, 9, 1, 3);  ax0 = [];
for c = 1:3
    R = RR{c};  ax = nexttile(tl);  hold on;  if isempty(ax0), ax0 = ax; end
    [mg, mo, sg, so] = deal(nan(1, numel(R.r_list)));
    for ir = 1:numel(R.r_list)
        g = R.res{ir, 1};  o = R.res{ir, 2};
        if strcmp(g.free_status, 'ok'), mg(ir) = g.free_stats.mean; sg(ir) = g.free_stats.std; end
        if o.built && strcmp(o.free_status, 'ok'), mo(ir) = o.free_stats.mean; so(ir) = o.free_stats.std; end
    end
    errorbar(R.r_list - 0.3, mg, sg, 's', 'Color', cG, 'MarkerSize', 4, 'DisplayName', 'POD-Galerkin');
    errorbar(R.r_list + 0.3, mo, so, 'o', 'Color', cO, 'MarkerSize', 4, 'DisplayName', 'OpInf');
    yline(R.fom_stats.mean, '-', 'Color', cF, 'DisplayName', 'FOM mean \pm std');
    yline(R.fom_stats.mean + R.fom_stats.std, ':', 'Color', cF, 'HandleVisibility', 'off');
    yline(R.fom_stats.mean - R.fom_stats.std, ':', 'Color', cF, 'HandleVisibility', 'off');
    grid on;  box on;  xlabel('r');
    if c == 1, ylabel('free-run \langleD\rangle / D_{lam}'); end
    title(sprintf('Re = %d', Res(c)));
end
shared_legend(ax0, 'south', 3);
save_fig(f, fig_dir, 'fig_r5_statistics_vs_r');

%% ---------------------------------------------------------------- Fig 6: energy budget (slides)
dfile = fullfile(res_dir, 'diag_Re40.mat');
if exist(dfile, 'file')
    L = load(dfile, 'G');  G = L.G;
    [f, tl] = newfig(7, 3.4, 9, 1, 2);
    nexttile(tl);  hold on;
    plot(G.cfg.r_list, G.proj_quiet, 'o-', 'Color', cO, 'MarkerSize', 5, 'DisplayName', 'quiet snapshots');
    plot(G.cfg.r_list, G.proj_burst, 's-', 'Color', cG, 'MarkerSize', 5, 'DisplayName', 'burst snapshots');
    grid on;  box on;  xlabel('r');  ylabel('projection error');  ylim([0 1]);
    legend('Location', 'southwest');  title('basis error, quiet vs burst');
    nexttile(tl);
    cl = G.closure;  cl = cl(arrayfun(@(d) ~isempty(d.r), cl));
    M = [[cl.Pt_q]; [cl.PG_q]; [cl.PO_q]; [cl.Pt_b]; [cl.PG_b]; [cl.PO_b]]';
    b = bar(categorical(arrayfun(@(d) sprintf('%d', d.r), cl, 'UniformOutput', false)), M);
    cols = [cF; cG; cO; cF; cG; cO];
    for k = 1:6, b(k).FaceColor = cols(k, :); if k > 3, b(k).FaceAlpha = 0.45; end, end
    grid on;  box on;  xlabel('r');  ylabel('\langle a \cdot da/dt \rangle');
    lg = legend({'true, quiet', 'Galerkin, quiet', 'OpInf, quiet', 'true, burst', 'Galerkin, burst', 'OpInf, burst'}, 'Location', 'southoutside');
    lg.NumColumns = 3;  lg.ItemTokenSize = [12 8];
    title('resolved energy budget');
    save_fig(f, fig_dir, 'fig_r6_energy_budget_Re40');
end

%% ---------------------------------------------------------------- Fig 7: unseen Re (paper, one column)
pfile = fullfile(res_dir, 'rom_param.mat');
if exist(pfile, 'file')
    L = load(pfile, 'P');  PP = L.P;
    [f, tl] = newfig(3.5, 2.5, 8, 1, 2);
    ir = numel(PP.r_list);
    [~, o] = sort(PP.allRe);
    ax0 = nexttile(tl);  hold on;
    for q = o
        x = PP.res{ir, q, 1};  y = PP.res{ir, q, 2};
        mk = 'o';  if ~x.seen, mk = 'p'; end
        errorbar(x.Re, x.fom_stats.mean, x.fom_stats.std, mk, 'Color', cF, 'MarkerFaceColor', cF, 'MarkerSize', 4, 'HandleVisibility', 'off');
        if strcmp(x.free_status, 'ok'), errorbar(x.Re - 0.8, x.stats.mean, x.stats.std, mk, 'Color', cG, 'MarkerSize', 4, 'HandleVisibility', 'off'); end
        if strcmp(y.free_status, 'ok'), errorbar(y.Re + 0.8, y.stats.mean, y.stats.std, mk, 'Color', cO, 'MarkerSize', 4, 'HandleVisibility', 'off'); end
    end
    plot(nan, nan, 'o', 'Color', cF, 'MarkerFaceColor', cF, 'MarkerSize', 4, 'DisplayName', 'FOM');
    plot(nan, nan, 'o-', 'Color', cG, 'MarkerSize', 4, 'DisplayName', 'POD-Galerkin');
    plot(nan, nan, 'o-', 'Color', cO, 'MarkerSize', 4, 'DisplayName', 'OpInf');
    plot(nan, nan, 'pk', 'MarkerFaceColor', 'k', 'MarkerSize', 7, 'DisplayName', 'unseen Re');
    grid on;  box on;  xlabel('Re');  ylabel('\langleD\rangle / D_{lam}');  xlim([17 43]);
    title('free-run dissipation');
    nexttile(tl);  hold on;
    for m = 1:2
        e = arrayfun(@(q) PP.res{ir, q, m}.coef_short, o);
        col = cG;  if m == 2, col = cO; end
        plot(PP.allRe(o), e, '-', 'Color', col, 'HandleVisibility', 'off');
        un = ~arrayfun(@(q) PP.res{ir, q, m}.seen, o);
        plot(PP.allRe(o(un)), e(un), 'p', 'Color', col, 'MarkerSize', 7, 'MarkerFaceColor', col, 'HandleVisibility', 'off');
    end
    grid on;  box on;  xlabel('Re');  ylabel('coefficient error');  xlim([17 43]);
    title(sprintf('first %g time units', PP.cfg.horizon));
    shared_legend(ax0, 'north', 2);
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

%% ---------------------------------------------------------------- Fig 8: vorticity fields at Re = 40 (slides)
R = RR{3};  ir = find(R.r_list >= 30, 1);  N = R.cfg.N;
S = matfile(R.cfg.train_file);
it = find(R.test);  k0 = it(1);  dt = R.t(2) - R.t(1);
lags = [0 2 5 10];
[f, tl] = newfig(6, 6.4, 9, numel(lags), 4);  tl.TileSpacing = 'tight';
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
    fields = {wF, wP};  ttl = {'FOM', 'projection'};
    for im = 1:2
        x = R.res{ir, im};
        if ~x.built, fields{end+1} = nan(N); ttl{end+1} = 'n/a'; continue; end   %#ok<SAGROW>
        fields{end+1} = reshape(R.wbar + R.V(:, 1:R.r_list(ir)) * Ahs{im}(:, q), N, N);  %#ok<SAGROW>
        ttl{end+1} = x.rom.type;                                                         %#ok<SAGROW>
    end
    lim = max(abs(wF(:)));
    for j = 1:4
        nexttile(tl, (q-1)*4 + j);
        imagesc(fields{j});  axis xy equal tight off;
        clim([-lim lim]);  colormap(gca, bwr(256));
        if q == 1, title(ttl{j}); end
        if j == 1
            text(-0.08, 0.5, sprintf('t_0 + %g', lags(q)), 'Units', 'normalized', 'Rotation', 90, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        end
    end
end
title(tl, sprintf('Re = 40, r = %d, test window', R.r_list(ir)), 'FontSize', 9);
save_fig(f, fig_dir, 'fig_r8_fields_Re40');

%% ---------------------------------------------------------------- Fig 9: lambda search (paper, 0.72 x full width)
[f, tl] = newfig(5.2, 2.0, 8, 1, 3);
rlam = [10 10 20];
for c = 1:3
    R = RR{c};  ir = find(R.r_list >= rlam(c), 1);  o = R.res{ir, 2};
    nexttile(tl);
    if ~o.built || ~isfield(o, 'lam_table'), axis off; continue; end
    T = o.lam_table;  S = T.score;  bl = isinf(S);  S(bl) = NaN;
    imagesc(log10(T.lam2s), log10(T.lam1s), log10(S), 'AlphaData', ~isnan(S));  axis xy;  colorbar;
    set(gca, 'Color', [0.92 0.92 0.92]);
    hold on;  [i1, i2] = find(bl);
    plot(log10(T.lam2s(i2)), log10(T.lam1s(i1)), 'kx', 'MarkerSize', 5, 'LineWidth', 1);
    plot(log10(o.info.lam(2)), log10(o.info.lam(1)), 'wp', 'MarkerSize', 9, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
    xlabel('log_{10} \lambda_2');
    if c == 1, ylabel('log_{10} \lambda_1'); end
    title(sprintf('Re = %d, r = %d', Res(c), R.r_list(ir)));
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
function [f, tl] = newfig(w, h, fs, rows, cols)
% Figure of w x h inches with fs-point text, holding a compact rows x cols tiled layout.
f = figure('Units', 'inches', 'Position', [1 1 w h], 'Color', 'w');
defs = {'DefaultAxesFontSize', fs; 'DefaultTextFontSize', fs; 'DefaultLegendFontSize', fs; ...
        'DefaultColorbarFontSize', fs; 'DefaultConstantLineFontSize', fs; ...
        'DefaultAxesTitleFontSizeMultiplier', 1; 'DefaultAxesLabelFontSizeMultiplier', 1; ...
        'DefaultAxesTitleFontWeight', 'normal'; 'DefaultLineLineWidth', 1.2};
for k = 1:size(defs, 1)
    try, set(f, defs{k, 1}, defs{k, 2}); catch, end      % skip any default this release lacks
end
tl = tiledlayout(f, rows, cols, 'TileSpacing', 'compact', 'Padding', 'compact');
end
function shared_legend(ax, where, ncol)
% One legend for the whole tiled layout, built from the DisplayNames in axes ax.
lg = legend(ax);  lg.NumColumns = ncol;  lg.ItemTokenSize = [18 10];
try, lg.Layout.Tile = where; catch, lg.Location = 'best'; end
end
function save_fig(f, d, name)
try, theme(f, 'light'); catch, end
set(f, 'Color', 'w');
exportgraphics(f, fullfile(d, [name '.png']), 'Resolution', 300);
end
function c = bwr(m)
h = floor(m/2);  up = linspace(0, 1, h)';
c = [[up up ones(h,1)]; [ones(m-h,1) flipud(linspace(0,1,m-h)') flipud(linspace(0,1,m-h)')]];
end
