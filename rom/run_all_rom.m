%% RUN_ALL_ROM  Reduced-order modelling of 2D Kolmogorov flow: full pipeline.
%
%   Step 1  generate_rom_data.m      fine-sampled FOM data          (~20 min, run once)
%   Step 2  kf_ops.m                 FOM operators (used by Galerkin and diagnostics)
%   Steps 3-6  run_rom_case.m        POD + POD-Galerkin + OpInf + evaluation, per Re
%   Step 7  run_rom_diagnostics.m    burst / closure diagnostics at Re = 40
%   Step 8  run_rom_param.m          one parametric ROM, tested at unseen Re = 25, 35
%   Figures make_rom_figures.m
%
% Set 'steps' to run only part of the pipeline. Results: ../rom_results/

clear; clc;
here = fileparts(mfilename('fullpath'));  root = fileparts(here);
addpath(here);  addpath(root);
dat = fullfile(root, 'rom_data');   old = fullfile(root, 'kf_results');
out = fullfile(root, 'rom_results');
if ~exist(out, 'dir'), mkdir(out); end

steps = {'Re20', 'Re30', 'Re40', 'diag', 'param'};   % any subset

% settings shared by all cases
%   basis 'invariant' : translation-invariant POD (pod_basis_inv)
%   n_aug             : OpInf training data + 8 random x-shifted copies
base = struct('N', 128, 'n', 4, 'ratio_min', 3, 'stride', 5, ...
    'lam1s', 10.^(-8:2:2), 'lam2s', 10.^(-8:1:2), 'tv_thr', 0.5, 'tol', 1e-8, ...
    'basis', 'invariant', 'n_aug', 8, 'rows_max', 15000);

if any(strcmp(steps, 'Re20'))
    cfg = base;  cfg.name = 'Re20';  cfg.Re = 20;
    cfg.train_file = fullfile(dat, 'fine_Re20_s1.mat');
    cfg.t_train = [100 250];  cfg.t_test = [250 400];
    cfg.ic_file = fullfile(old, 'case_Re20_s2.mat');        % different attractor (steady)
    cfg.r_list = [4 6 8 10 14 20 30];
    cfg.horizon = 10;  cfg.Hv = 60;  cfg.nstarts_sel = 5;  cfg.nstarts_test = 5;
    cfg.outfile = fullfile(out, 'rom_Re20.mat');  cfg.logfile = fullfile(out, 'log_Re20.txt');
    run_rom_case(cfg);
end

if any(strcmp(steps, 'Re30'))
    cfg = base;  cfg.name = 'Re30';  cfg.Re = 30;
    cfg.train_file = fullfile(dat, 'fine_Re30_s1.mat');
    cfg.t_train = [100 250];  cfg.t_test = [250 400];
    cfg.ic_file = fullfile(old, 'case_Re30_s2.mat');
    cfg.r_list = [4 6 8 10 14 20 30];
    cfg.horizon = 10;  cfg.Hv = 60;  cfg.nstarts_sel = 5;  cfg.nstarts_test = 5;
    cfg.outfile = fullfile(out, 'rom_Re30.mat');  cfg.logfile = fullfile(out, 'log_Re30.txt');
    run_rom_case(cfg);
end

if any(strcmp(steps, 'Re40'))
    cfg = base;  cfg.name = 'Re40';  cfg.Re = 40;
    cfg.train_file = fullfile(dat, 'fine_Re40_s1.mat');
    cfg.t_train = [100 800];  cfg.t_test = [800 1100];      % both windows contain bursts
    cfg.ic_file = fullfile(dat, 'fine_Re40_s2.mat');
    cfg.r_list = [4 6 8 10 14 20 30 40 60 90 120];      % r > 60: Galerkin only (OpInf has too many unknowns)
    cfg.horizon = 5;  cfg.Hv = 40;  cfg.nstarts_sel = 8;  cfg.nstarts_test = 10;  cfg.rows_max = 20000;
    cfg.outfile = fullfile(out, 'rom_Re40.mat');  cfg.logfile = fullfile(out, 'log_Re40.txt');
    run_rom_case(cfg);
end

if any(strcmp(steps, 'diag'))
    cfg = struct('N', 128, 'n', 4, 'Re', 40, 'data_file', fullfile(dat, 'fine_Re40_s1.mat'), ...
        'results_file', fullfile(out, 'rom_Re40.mat'), 'r_list', [10 20 30 40 60], ...
        'thr', 0.15, 'stride', 5, 'outfile', fullfile(out, 'diag_Re40.mat'), ...
        'logfile', fullfile(out, 'log_diag.txt'));
    run_rom_diagnostics(cfg);
end

if any(strcmp(steps, 'param'))
    cfg = struct('N', 128, 'n', 4, ...
        'train_files', {{fullfile(dat, 'fine_Re20_s1.mat'), fullfile(dat, 'fine_Re30_s1.mat'), fullfile(dat, 'fine_Re40_s1.mat')}}, ...
        'train_Re', [20 30 40], 't_train', [100 250; 100 250; 100 800], ...
        'pod_counts', 300, 'deriv_rows', 3000, 'n_aug', 8, ...
        'test_files', {{fullfile(old, 'case_Re25_s1.mat'), fullfile(old, 'case_Re35_s1.mat')}}, ...
        'test_Re', [25 35], 'r_list', [10 20 30], ...
        'lam1s', 10.^(-8:2:2), 'lam2s', 10.^(-8:1:2), 'horizon', 5, 'Hv', 40, ...
        'nstarts', 5, 'tv_thr', 0.5, 'tol', 1e-8, ...
        'outfile', fullfile(out, 'rom_param.mat'), 'logfile', fullfile(out, 'log_param.txt'));
    run_rom_param(cfg);
end
fprintf('\nrun_all_rom finished: %s\n', datestr(now));
