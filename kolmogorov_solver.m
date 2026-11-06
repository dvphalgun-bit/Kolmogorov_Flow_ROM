function out = kolmogorov_solver(p)
% KOLMOGOROV_SOLVER  Pseudo-spectral full-order model of 2D Kolmogorov flow.
%
% Governing equation (vorticity form, non-dimensional, doubly periodic
% domain [0, 2*pi]^2, Chandler & Kerswell 2013, J. Fluid Mech. 722):
%
%   dw/dt + u.grad(w) = (1/Re) lap(w) - n cos(n y)
%   lap(psi) = -w ,   u = d(psi)/dy ,   v = -d(psi)/dx
%
% The forcing term -n cos(n y) is the curl of the body force sin(n y) e_x.
% Laminar steady solution:  u = (Re/n^2) sin(n y),  w = -(Re/n) cos(n y).
%
% Discretisation:
%   - Fourier pseudo-spectral in x and y on an N x N grid
%     (state dimension = N^2 vorticity values)
%   - nonlinear term evaluated in physical space, 2/3-rule dealiasing
%   - time stepping: integrating-factor 4th-order Runge-Kutta (IFRK4);
%     the viscous term is integrated exactly, so only the advective
%     term limits the time step (CFL)
%
% Input  p : struct with fields (defaults in brackets)
%   Re        Reynolds number                                   [40]
%   N         grid points per direction                         [128]
%   n         forcing wavenumber                                [4]
%   dt        time step                                         [0.005]
%   T         final time                                        [400]
%   t_rec     time at which snapshot recording starts            [100]
%   dt_save   snapshot sampling interval (multiple of dt)        [0.5]
%   dt_diag   sampling interval of scalar diagnostics            [0.05]
%   seed      random seed for the initial perturbation          [1]
%   amp       perturbation amplitude (fraction of laminar max)  [0.1]
%   store     true -> keep snapshot matrix W                    [true]
%   verbose   print progress                                    [true]
%
% Output out : struct with
%   W        N^2 x m snapshot matrix (single), column j = w(:) at t_snap(j)
%   t_snap   1 x m snapshot times
%   t_diag, E, D, I   kinetic energy, dissipation, energy input vs time
%   k_spec, E_spec    time-averaged radial energy spectrum (recorded window)
%   w_final, x, cfl_max, runtime, p

% ------------------------------------------------------------ defaults
def = struct('Re',40,'N',128,'n',4,'dt',0.005,'T',400,'t_rec',100, ...
             'dt_save',0.5,'dt_diag',0.05,'seed',1,'amp',0.1, ...
             'store',true,'verbose',true);
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(p, f{i}), p.(f{i}) = def.(f{i}); end
end
Re = p.Re; N = p.N; n = p.n; dt = p.dt;

% ------------------------------------------------------------ grid
Lx = 2*pi;
x  = (0:N-1)' * Lx/N;
[X, Y] = meshgrid(x, x);          % rows = y, columns = x
dx = Lx/N;

k = [0:N/2-1, -N/2:-1]';          % integer wavenumbers (domain 2*pi)
[KX, KY] = meshgrid(k, k);
K2   = KX.^2 + KY.^2;
K2i  = 1 ./ K2;  K2i(1,1) = 0;    % inverse Laplacian, zero mean mode
mask = (abs(KX) < N/3) & (abs(KY) < N/3);   % 2/3-rule dealiasing
iKX = 1i*KX;  iKY = 1i*KY;

% ------------------------------------------------------------ forcing
F  = -n * cos(n*Y);
Fh = fft2(F);
fx = sin(n*Y);                    % body force (x-component) for energy input

% ------------------------------------------------------------ IFRK4 factors
Lin = -K2 / Re;                   % linear (viscous) operator in Fourier space
E1  = exp(Lin*dt/2);
E2  = E1.^2;

% ------------------------------------------------------------ initial condition
% laminar state + random large-scale perturbation (|k| <= 4)
rng(p.seed);
w_lam = -(Re/n) * cos(n*Y);
r  = randn(N);
rp = real(ifft2(fft2(r) .* (K2 <= 16 & K2 > 0)));
w0 = w_lam + p.amp * (Re/n) * rp / max(abs(rp(:)));
wh = fft2(w0);

% ------------------------------------------------------------ bookkeeping
nsteps     = round(p.T / dt);
save_every = round(p.dt_save / dt);
diag_every = round(p.dt_diag / dt);
rec_step   = round(p.t_rec / dt);

n_snap = floor((nsteps - rec_step)/save_every) + 1;
n_diag = floor(nsteps/diag_every) + 1;

if p.store, W = zeros(N*N, n_snap, 'single'); else, W = []; end
t_snap = zeros(1, n_snap);
t_diag = zeros(1, n_diag); E = t_diag; D = t_diag; I = t_diag;

kmax_bin = ceil(sqrt(2)*N/2) + 1;
kbin     = round(sqrt(K2)) + 1;
E_spec   = zeros(kmax_bin, 1);
n_spec   = 0;
cfl_max  = 0;

js = 0; jd = 0;
tic;
for s = 0:nsteps
    t = s*dt;

    % ---- diagnostics (energy, dissipation, input, CFL)
    if mod(s, diag_every) == 0
        ph = wh .* K2i;
        u  = real(ifft2(iKY .* ph));
        v  = real(ifft2(-iKX .* ph));
        w  = real(ifft2(wh));
        jd = jd + 1;
        t_diag(jd) = t;
        E(jd) = 0.5 * mean(u(:).^2 + v(:).^2);     % kinetic energy
        D(jd) = mean(w(:).^2) / Re;                % dissipation rate
        I(jd) = mean(u(:) .* fx(:));               % energy input rate
        cfl_max = max(cfl_max, max(abs(u(:)) + abs(v(:))) * dt / dx);
        if ~all(isfinite(wh(:)))
            error('Solution blew up at t = %.3f (reduce dt)', t);
        end
    end

    % ---- snapshots + spectrum (recording window only)
    if s >= rec_step && mod(s - rec_step, save_every) == 0
        js = js + 1;
        t_snap(js) = t;
        if p.store
            W(:, js) = single(reshape(real(ifft2(wh)), [], 1));
        end
        Eden   = 0.5 * abs(wh).^2 .* K2i / N^4;      % energy density per mode
        E_spec = E_spec + accumarray(kbin(:), Eden(:), [kmax_bin 1]);
        n_spec = n_spec + 1;
    end

    if s == nsteps, break; end

    % ---- IFRK4 step
    k1 = dt * rhs_nl(wh);
    k2 = dt * rhs_nl(E1 .* (wh + k1/2));
    k3 = dt * rhs_nl(E1 .* wh + k2/2);
    k4 = dt * rhs_nl(E2 .* wh + E1 .* k3);
    wh = E2 .* wh + (E2 .* k1 + 2*E1 .* (k2 + k3) + k4) / 6;

    if p.verbose && mod(s+1, round(nsteps/10)) == 0
        fprintf('   Re = %g : t = %6.1f / %g  (%.0f s)\n', Re, (s+1)*dt, p.T, toc);
    end
end
runtime = toc;

out.W = W(:, 1:js);           out.t_snap = t_snap(1:js);
out.t_diag = t_diag(1:jd);    out.E = E(1:jd);  out.D = D(1:jd);  out.I = I(1:jd);
out.k_spec = (0:kmax_bin-1)'; out.E_spec = E_spec / max(n_spec, 1);
out.w_final = real(ifft2(wh)); out.x = x;
out.cfl_max = cfl_max;        out.runtime = runtime;  out.p = p;

    % -------------------------------------------------------- nested RHS
    function Nh = rhs_nl(wh_in)
        % Fourier transform of  -(u.grad w) + forcing   (dealiased)
        psih = wh_in .* K2i;
        uu = real(ifft2( iKY .* psih));
        vv = real(ifft2(-iKX .* psih));
        wx = real(ifft2( iKX .* wh_in));
        wy = real(ifft2( iKY .* wh_in));
        Nh = -fft2(uu.*wx + vv.*wy) .* mask + Fh;
    end
end
