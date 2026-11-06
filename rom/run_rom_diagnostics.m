function G = run_rom_diagnostics(cfg)
% RUN_ROM_DIAGNOSTICS  Step 7: what do the ROMs do during dissipation bursts?
% (Re = 40; answers the reviewer's second comment)
%
% 1. Basis: projection error and projected dissipation, burst vs quiet snapshots.
% 2. Closure: the exact effect of the discarded POD modes on the retained ones,
%       tau(a) = V' rhs_FOM(w) - f_Galerkin(V'(w - wbar)),
%    computed from the FOM operators on training snapshots.  Its energy
%    transfer a.tau (the rate at which the discarded modes change |a|^2/2)
%    shows whether truncation removes a dissipation mechanism.
% 3. Energy budget: the true rate of change of the resolved fluctuation energy
%    a.(da/dt) versus the rate predicted by Galerkin and by OpInf, in quiet and
%    burst periods. The Galerkin deficit (true - Galerkin) is the energy drained
%    by the discarded modes; the 'share' is the fraction of it OpInf supplies.
%    (A correlation between f_OpInf - f_Galerkin and tau is NOT used: when the
%    closure dominates, both are close to -f_Galerkin and the correlation is
%    high regardless of what OpInf learned.)
% 4. Symmetric part of A: largest eigenvalue for Galerkin vs OpInf.
%
% cfg: N, n, Re, data_file, results_file, r_list, thr, stride, outfile

ops  = kf_ops(cfg.N, cfg.n);
Dlam = cfg.Re / (2*cfg.n^2);
L = load(cfg.results_file, 'R');  R = L.R;  clear L
S = load(cfg.data_file, 'W', 't_snap', 'diagn');
t = S.t_snap;
Dn = interp1(S.diagn.t_diag, S.diagn.D, t) / Dlam;
burst = Dn > cfg.thr;
V = R.V;  wbar = R.wbar;
say(cfg, '\n===== burst / closure diagnostics (Re = %g) =====\n', cfg.Re);
say(cfg, '%d of %d snapshots are bursts (D > %.2f D_lam)\n', nnz(burst), numel(t), cfg.thr);

% ---------------------------------------------------------------- 1. basis
[~, jp] = max(Dn);                                     % strongest burst
G = struct('cfg', cfg, 't', t, 'Dn', Dn, 'burst', burst, 't_peak', t(jp));
for ir = 1:numel(cfg.r_list)
    r = cfg.r_list(ir);  Vr = V(:, 1:r);
    eq = zeros(1, 0);  eb = eq;
    for j0 = 1:500:numel(t)
        j = j0:min(numel(t), j0+499);
        X = double(S.W(:, j)) - wbar;  Ap = Vr' * X;
        e = sqrt(max(0, sum(X.^2,1) - sum(Ap.^2,1)) ./ sum(X.^2,1));
        eq = [eq, e(~burst(j))];  eb = [eb, e(burst(j))];   %#ok<AGROW>
    end
    wp = wbar + Vr * (Vr' * (double(S.W(:, jp)) - wbar));
    G.proj_quiet(ir) = mean(eq);  G.proj_burst(ir) = mean(eb);
    G.Dpeak_proj(ir) = mean(wp.^2) / cfg.Re / Dlam;
end
G.Dpeak_true = Dn(jp);
say(cfg, '  r | proj. err quiet | proj. err burst | D/D_lam at t = %.1f: %.3f ->\n', t(jp), Dn(jp));
for ir = 1:numel(cfg.r_list)
    say(cfg, ' %2d |     %.3f       |     %.3f       |  %.3f\n', cfg.r_list(ir), ...
        G.proj_quiet(ir), G.proj_burst(ir), G.Dpeak_proj(ir));
end

% ---------------------------------------------------------------- 2-4. closure
tr = find(R.train);  js = tr(1:cfg.stride:end);        % training snapshots only
Wj = double(S.W(:, js));  bj = burst(js);
clear S
Rhs = zeros(size(Wj));
for j0 = 1:200:numel(js)
    j = j0:min(numel(js), j0+199);
    Rhs(:, j) = ops.rhs(Wj(:, j), cfg.Re);             % exact FOM time derivative
end
Xj = Wj - wbar;  clear Wj
r_all = cellfun(@(x) x.r, R.res(:, 1))';
say(cfg, '\n  Energy budget of the resolved fluctuations, d(|a|^2/2)/dt = a . da/dt  (mean over training snapshots)\n');
say(cfg, '  r | rhs error: Gal   OpInf | quiet: true    Gal      OpInf  | burst: true    Gal      OpInf  | OpInf share of closure drain q/b | max eig sym(A) Gal  OpInf\n');
for ir = 1:numel(cfg.r_list)
    r = min(r_all(r_all >= cfg.r_list(ir)));  k = find(r_all == r, 1);
    if isempty(k) || ~R.res{k, 2}.built, continue; end
    Gr = R.res{k, 1}.rom;   Or = R.res{k, 2}.rom;
    Vr = V(:, 1:r);  a = Vr' * Xj;
    dtrue = Vr' * Rhs;                                 % exact projected dynamics
    fG = rom_rhs(Gr, a);   fO = rom_rhs(Or, a);
    Pt = sum(a .* dtrue, 1);  PG = sum(a .* fG, 1);  PO = sum(a .* fO, 1);
    d.r = r;
    d.err_G = norm(fG - dtrue, 'fro') / norm(dtrue, 'fro');   % = |closure| / |da/dt|
    d.err_O = norm(fO - dtrue, 'fro') / norm(dtrue, 'fro');
    d.Pt_q = mean(Pt(~bj));  d.PG_q = mean(PG(~bj));  d.PO_q = mean(PO(~bj));
    d.Pt_b = mean(Pt(bj));   d.PG_b = mean(PG(bj));   d.PO_b = mean(PO(bj));
    % fraction of the closure energy drain (true - Galerkin) that OpInf supplies
    d.share_q = (d.PO_q - d.PG_q) / (d.Pt_q - d.PG_q);
    d.share_b = (d.PO_b - d.PG_b) / (d.Pt_b - d.PG_b);
    d.eigG = max(eig((Gr.A + Gr.A')/2));
    d.eigO = max(eig((Or.A + Or.A')/2));
    d.Pt = Pt;  d.PG = PG;  d.PO = PO;  d.bj = bj;  d.tj = t(js);
    G.closure(ir) = d;
    say(cfg, ' %2d |        %5.2f  %5.2f  | %8.1f %8.1f %8.1f | %8.1f %8.1f %8.1f |          %5.2f / %5.2f           |   %7.3f %7.3f\n', ...
        r, d.err_G, d.err_O, d.Pt_q, d.PG_q, d.PO_q, d.Pt_b, d.PG_b, d.PO_b, d.share_q, d.share_b, d.eigG, d.eigO);
end
save(cfg.outfile, 'G', '-v7.3');
say(cfg, 'saved %s\n', cfg.outfile);
end

function say(cfg, varargin)
% print to the console and append to cfg.logfile (live progress log)
s = sprintf(varargin{:});
fprintf('%s', s);
if isfield(cfg, 'logfile') && ~isempty(cfg.logfile)
    fid = fopen(cfg.logfile, 'a');  fprintf(fid, '%s', s);  fclose(fid);
end
end
