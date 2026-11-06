%% AS5401 Project Proposal - Full-order model of 2D Kolmogorov flow
% Runs the pseudo-spectral FOM (kolmogorov_solver.m) for several Reynolds
% numbers and initial conditions, saves the snapshot data for the ROM
% phase, and produces the figures/numbers needed for Section 3.
%
% Outputs (folder kf_results/):
%   case_ReXX_sY.mat      snapshot matrix W (N^2 x m, single), times, diagnostics
%   fig1 ... fig7 .png    figures for the proposal
%   summary.txt           all printed numbers (diary)
%
% Runtime with the default settings: roughly 20-40 min on a laptop.
% Set quick_test = true first to check everything works (~1-2 min).

clear; clc; close all;

quick_test = false;
reuse      = false;   % true: reuse existing case_*.mat files instead of re-simulating

%% ------------------------------------------------------------ settings
base.n       = 4;          % forcing wavenumber
base.amp     = 0.1;        % initial perturbation (fraction of laminar max)
base.dt_diag = 0.05;
base.verbose = true;
if quick_test
    base.N = 64;  base.T = 60;  base.t_rec = 20;  base.dt_save = 0.5;
else
    base.N = 128; base.T = 400; base.t_rec = 100; base.dt_save = 0.5;
end

% Cases: Re, seed, role
%   train  : data used to build the ROMs (later phase)
%   testIC : same Re, different random initial condition
%   testRe : Reynolds number not seen in training (parametric test)
cases = struct('Re',  {8,       20,      30,      40,      20,       30,       40,       25,       35}, ...
               'seed',{1,       1,       1,       1,       2,        2,        2,        1,        1}, ...
               'role',{'train','train','train','train','testIC','testIC','testIC','testRe','testRe'});

outdir = 'kf_results';
if ~exist(outdir, 'dir'), mkdir(outdir); end
diary off; if exist(fullfile(outdir,'summary.txt'),'file'), delete(fullfile(outdir,'summary.txt')); end
diary(fullfile(outdir, 'summary.txt'));

fprintf('Kolmogorov flow FOM:  N = %d  (state dimension %d), n = %d\n', base.N, base.N^2, base.n);
fprintf('T = %g, recording from t = %g every %g -> %d snapshots per case\n\n', ...
        base.T, base.t_rec, base.dt_save, floor((base.T-base.t_rec)/base.dt_save)+1);

%% ------------------------------------------------------------ run all cases
R = cell(1, numel(cases));
for c = 1:numel(cases)
    p = base;
    p.Re = cases(c).Re;  p.seed = cases(c).seed;
    p.dt = 0.01;  if p.Re > 30, p.dt = 0.005; end    % CFL-safe step
    p.store = true;
    fprintf('Case %d/%d: Re = %g, seed = %d (%s), dt = %g\n', c, numel(cases), p.Re, p.seed, cases(c).role, p.dt);
    fname = fullfile(outdir, sprintf('case_Re%02d_s%d.mat', p.Re, p.seed));
    if reuse && exist(fname, 'file')
        S = load(fname, 'W', 'diagn');
        out = S.diagn;  out.W = S.W;
        out.w_final = reshape(double(S.W(:, end)), p.N, p.N);
        fprintf('   loaded existing %s\n', fname);
    else
        out = kolmogorov_solver(p);
    end
    out.role = cases(c).role;
    R{c} = out;

    % save snapshot data for the ROM phase (without the large W in R)
    W = out.W; t_snap = out.t_snap; params = out.p; role = out.role; x = out.x; %#ok<NASGU>
    diagn = rmfield(out, {'W','w_final'}); %#ok<NASGU>
    save(fname, 'W', 't_snap', 'params', 'role', 'x', 'diagn', '-v7.3');
    R{c}.W = [];                       % free memory, keep only what is plotted
    if strcmp(cases(c).role,'train') && p.Re == 40
        R{c}.W = out.W;                % keep Re = 40 snapshots for burst figure
    end
    fprintf('   done in %.1f s, max CFL = %.2f, saved %s (%d x %d)\n\n', ...
            out.runtime, out.cfl_max, fname, size(out.W,1), size(out.W,2));
end

train = find(strcmp({cases.role}, 'train'));
N = base.N; n = base.n;
cmap = bwr_colormap(256);

%% ------------------------------------------------------------ checks
fprintf('================ VALIDATION / PHYSICAL CHECKS ================\n');
% (1) laminar solution at Re = 8 (below Re_c ~ 9.9 for n = 4)
c8 = train([cases(train).Re] == 8);
[~, Yg] = meshgrid(R{c8}.x, R{c8}.x);
w_lam = -(8/n) * cos(n*Yg);
fprintf('Re = 8 : ||w - w_laminar|| / ||w_laminar|| at t = T : %.2e\n', ...
        norm(R{c8}.w_final(:) - w_lam(:)) / norm(w_lam(:)));

% (2) energy balance dE/dt = I - D, and long-time balance <I> = <D>
fprintf('\n  Re  seed | <D>/D_lam | min D/D_lam | max D/D_lam | <E>/E_lam | |<I>-<D>|/<D> | dE/dt residual | max CFL | E(k_max)/E_peak\n');
for c = 1:numel(cases)
    o = R{c}; Re = o.p.Re;
    D_lam = Re/(2*n^2);  E_lam = Re^2/(4*n^4);
    m  = o.t_diag >= o.p.t_rec;
    dEdt = gradient(o.E, o.t_diag);
    res  = norm(dEdt(m) - (o.I(m) - o.D(m))) / norm(o.D(m));
    kd   = floor(N/3) - 1;                     % last fully resolved shell
    Ek   = o.E_spec(2:end);
    fprintf('  %2g   %d   |  %7.4f  |   %7.4f   |   %7.4f   |  %7.4f  |   %8.2e    |    %8.2e    |  %5.2f  |   %8.2e\n', ...
        Re, o.p.seed, mean(o.D(m))/D_lam, min(o.D(m))/D_lam, max(o.D(m))/D_lam, ...
        mean(o.E(m))/E_lam, abs(mean(o.I(m))-mean(o.D(m)))/mean(o.D(m)), res, o.cfl_max, ...
        o.E_spec(kd+1)/max(Ek));
end

% (3) dominant frequency of D(t) (periodic vs. chaotic)
fprintf('\nDominant frequency of D(t) in the recorded window (share of fluctuation power):\n');
for c = train
    o = R{c}; m = o.t_diag >= o.p.t_rec;
    d = o.D(m) - mean(o.D(m));
    if std(d) < 1e-6*mean(o.D(m)), fprintf('  Re = %2g : steady (std(D)/<D> = %.1e)\n', o.p.Re, std(d)/mean(o.D(m))); continue; end
    P = abs(fft(d)).^2;  P = P(2:floor(end/2));
    fr = (1:numel(P)) / (numel(d)*o.p.dt_diag);
    [Pm, im] = max(P);
    fprintf('  Re = %2g : f = %.4f (period %.2f), %.1f %% of power\n', o.p.Re, fr(im), 1/fr(im), 100*Pm/sum(P));
end

%% ------------------------------------------------------------ Fig 1: vorticity fields
fig = figure('Position',[50 50 1300 340],'Color','w');
tl = tiledlayout(1, numel(train), 'TileSpacing','compact');
for j = 1:numel(train)
    o = R{train(j)}; nexttile;
    wf = o.w_final; lim = max(abs(wf(:)));
    imagesc(o.x, o.x, wf); axis xy equal tight; colormap(gca, cmap); caxis([-lim lim]);
    cb = colorbar; cb.Label.String = '\omega';
    xlabel('x'); ylabel('y'); title(sprintf('Re = %g, t = %g', o.p.Re, o.p.T));
end
title(tl, 'Vorticity field \omega(x,y) at the final time');
save_fig(fig, outdir, 'fig1_vorticity_fields');

%% ------------------------------------------------------------ Fig 2: E(t), D(t)
fig = figure('Position',[50 50 1100 800],'Color','w');
tl = tiledlayout(numel(train), 2, 'TileSpacing','compact');
for j = 1:numel(train)
    o = R{train(j)}; Re = o.p.Re;
    D_lam = Re/(2*n^2);  E_lam = Re^2/(4*n^4);
    nexttile; plot(o.t_diag, o.E/E_lam, 'b', 'LineWidth', 1); grid on;
    xline(o.p.t_rec, 'k--'); ylabel('E / E_{lam}'); title(sprintf('Re = %g : kinetic energy', Re));
    nexttile; plot(o.t_diag, o.D/D_lam, 'r', 'LineWidth', 1); grid on;
    xline(o.p.t_rec, 'k--'); ylabel('D / D_{lam}'); title(sprintf('Re = %g : dissipation rate', Re));
end
xlabel(tl, 'time t (non-dimensional)');
title(tl, 'Energy and dissipation histories (dashed line: start of recorded snapshots)');
save_fig(fig, outdir, 'fig2_energy_dissipation_time');

%% ------------------------------------------------------------ Fig 3: I-D phase portraits
unsteady = train([cases(train).Re] > 8);              % Re = 8 is steady laminar flow
fig = figure('Position',[50 50 1100 360],'Color','w');
tl = tiledlayout(1, numel(unsteady), 'TileSpacing','compact');
for j = 1:numel(unsteady)
    o = R{unsteady(j)}; Re = o.p.Re; D_lam = Re/(2*n^2);
    m = o.t_diag >= o.p.t_rec;
    nexttile; plot(o.I(m)/D_lam, o.D(m)/D_lam, 'k', 'LineWidth', 0.6); hold on;
    lims = [min([o.I(m) o.D(m)]) max([o.I(m) o.D(m)])]/D_lam;
    plot(lims, lims, 'r--'); grid on; axis tight;
    xlabel('I / D_{lam}'); ylabel('D / D_{lam}'); title(sprintf('Re = %g', Re));
end
title(tl, 'Energy input I vs. dissipation D (recorded window); red: I = D  (Re = 8 is steady: a single point)');
save_fig(fig, outdir, 'fig3_phase_portraits');

%% ------------------------------------------------------------ Fig 4: energy spectra
fig = figure('Position',[50 50 650 480],'Color','w');
cols = lines(numel(train));
for j = 1:numel(train)
    o = R{train(j)};
    loglog(o.k_spec(2:end), o.E_spec(2:end), 'Color', cols(j,:), 'LineWidth', 1.5, ...
           'DisplayName', sprintf('Re = %g', o.p.Re)); hold on;
end
xline(N/3, 'k--', '2/3 dealiasing cut-off', 'HandleVisibility','off');
grid on; xlabel('wavenumber |k|'); ylabel('E(|k|)');
legend('Location','southwest'); title('Time-averaged energy spectrum (resolution check)');
save_fig(fig, outdir, 'fig4_energy_spectra');

%% ------------------------------------------------------------ Fig 5: snapshot SVD (preliminary POD)
fprintf('\n================ PRELIMINARY POD (training cases) ================\n');
fprintf('  Re | r for 90%% | r for 99%% | r for 99.9%% energy  (mean-subtracted snapshots)\n');
fig = figure('Position',[50 50 1100 420],'Color','w');
tiledlayout(1, 2, 'TileSpacing','compact');
ax1 = nexttile; ax2 = nexttile;
for j = 1:numel(train)
    o = R{train(j)}; Re = o.p.Re;
    S = load(fullfile(outdir, sprintf('case_Re%02d_s%d.mat', Re, o.p.seed)), 'W');
    Wm = double(S.W); Wm = Wm - mean(Wm, 2);
    sv = svd(Wm, 'econ'); en = cumsum(sv.^2) / sum(sv.^2);
    m = o.t_diag >= o.p.t_rec;
    if std(o.D(m)) < 1e-6*mean(o.D(m))          % steady flow: POD of round-off is meaningless
        fprintf('  %2g | steady (std(D)/<D> = %.1e), skipped\n', Re, std(o.D(m))/mean(o.D(m))); continue;
    end
    fprintf('  %2g |   %4d    |   %4d    |   %4d\n', Re, find(en>=0.9,1), find(en>=0.99,1), find(en>=0.999,1));
    nr = min(100, numel(sv));
    semilogy(ax1, 1:nr, sv(1:nr)/sv(1), 'o-', 'Color', cols(j,:), 'MarkerSize', 3, 'DisplayName', sprintf('Re = %g', Re)); hold(ax1,'on');
    plot(ax2, 1:nr, en(1:nr), 'o-', 'Color', cols(j,:), 'MarkerSize', 3, 'DisplayName', sprintf('Re = %g', Re)); hold(ax2,'on');
end
grid(ax1,'on'); xlabel(ax1,'index i'); ylabel(ax1,'\sigma_i / \sigma_1'); title(ax1,'Normalised singular values'); legend(ax1,'Location','northeast');
grid(ax2,'on'); xlabel(ax2,'number of POD modes r'); ylabel(ax2,'cumulative energy'); title(ax2,'Captured fluctuation energy'); legend(ax2,'Location','southeast');
yline(ax2, 0.99, 'k--', '99%', 'HandleVisibility','off');
save_fig(fig, outdir, 'fig5_pod_singular_values');

%% ------------------------------------------------------------ Fig 6: burst sequence at Re = 40
c40 = train([cases(train).Re] == 40);
o = R{c40}; D_lam = 40/(2*n^2);
Ds = interp1(o.t_diag, o.D, o.t_snap);
Ds_in = Ds; Ds_in(o.t_snap < o.p.t_rec + 25 | o.t_snap > o.p.T - 15) = -Inf;
[~, ib] = max(Ds_in);                                   % strongest burst away from window edges
idx = max(1, min(numel(o.t_snap), ib + round([-20 -6 0 10]/o.p.dt_save)));
fig = figure('Position',[50 50 1300 650],'Color','w');
tl = tiledlayout(2, 4, 'TileSpacing','compact');
lim = max(abs(double(o.W(:, idx))), [], 'all');
for j = 1:4
    nexttile;
    imagesc(o.x, o.x, reshape(double(o.W(:, idx(j))), N, N)); axis xy equal tight;
    colormap(gca, cmap); caxis([-lim lim]); colorbar;
    title(sprintf('t = %.0f,  D/D_{lam} = %.2f', o.t_snap(idx(j)), Ds(idx(j))/D_lam));
end
nexttile([1 4]);
plot(o.t_diag, o.D/D_lam, 'r'); hold on; grid on;
plot(o.t_snap(idx), Ds(idx)/D_lam, 'ko', 'MarkerFaceColor', 'k');
xlim([o.p.t_rec o.p.T]); xlabel('time t'); ylabel('D / D_{lam}');
title('Re = 40: dissipation history; dots mark the vorticity snapshots shown above (before, rise, peak, decay of a burst)');
save_fig(fig, outdir, 'fig6_burst_sequence_Re40');

%% ------------------------------------------------------------ Fig 7: train vs test cases
fig = figure('Position',[50 50 1100 600],'Color','w');
tiledlayout(2, 1, 'TileSpacing','compact');
nexttile; hold on;
for c = 1:numel(cases)
    if ~strcmp(cases(c).role, 'testIC') && ~(strcmp(cases(c).role,'train') && cases(c).Re > 8), continue; end
    o = R{c}; Re = o.p.Re; D_lam = Re/(2*n^2);
    ls = '-'; if strcmp(cases(c).role,'testIC'), ls = '--'; end
    plot(o.t_diag, o.D/D_lam, ls, 'LineWidth', 1, 'DisplayName', sprintf('Re = %g, seed %d (%s)', Re, o.p.seed, o.role));
end
grid on; ylabel('D / D_{lam}'); legend('Location','eastoutside'); title('Training (solid) vs. different initial condition (dashed)');
nexttile; hold on;
for c = 1:numel(cases)
    o = R{c}; Re = o.p.Re; if Re == 8, continue; end
    D_lam = Re/(2*n^2); m = o.t_diag >= o.p.t_rec;
    mk = 'o'; if strcmp(o.role,'testRe'), mk = 's'; end; if strcmp(o.role,'testIC'), mk = '^'; end
    errorbar(Re, mean(o.D(m))/D_lam, std(o.D(m))/D_lam, mk, 'LineWidth', 1.2, 'MarkerSize', 8, ...
             'DisplayName', sprintf('Re = %g (%s)', Re, o.role));
end
grid on; xlabel('Re'); ylabel('<D> / D_{lam}  \pm  std'); legend('Location','eastoutside');
title('Mean dissipation vs. Re: circle = train, triangle = new initial condition, square = unseen Re');
save_fig(fig, outdir, 'fig7_train_test_cases');

%% ------------------------------------------------------------ data description
fprintf('\n================ DATASET ================\n');
fprintf('Each case_ReXX_sY.mat contains W of size %d x %d (single precision):\n', N^2, ...
        floor((base.T-base.t_rec)/base.dt_save)+1);
fprintf('  row i    = vorticity at grid point i of the %d x %d grid (column-major, reshape(W(:,j),N,N))\n', N, N);
fprintf('  column j = snapshot at time t_snap(j), t in [%g, %g], dt_save = %g\n', base.t_rec, base.T, base.dt_save);
fprintf('Total FOM wall-clock time: %.1f s\n', sum(cellfun(@(o) o.runtime, R)));
diary off;

%% ------------------------------------------------------------ helpers
function save_fig(fig, outdir, name)
    % force a light theme (MATLAB R2025a+ dark mode gives grey-on-black plots)
    try, theme(fig, 'light'); catch, end
    set(fig, 'Color', 'w');
    f = fullfile(outdir, [name '.png']);
    if exist('exportgraphics', 'file')
        exportgraphics(fig, f, 'Resolution', 200);
    else
        print(fig, f, '-dpng', '-r200');
    end
end

function c = bwr_colormap(m)
    % blue - white - red diverging colormap
    h = floor(m/2);
    up = linspace(0, 1, h)';
    c = [[up up ones(h,1)]; [ones(m-h,1) flipud(linspace(0,1,m-h)') flipud(linspace(0,1,m-h)')]];
end
