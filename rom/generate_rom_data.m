%% GENERATE_ROM_DATA  Step 1 of the ROM study: fine-sampled training data.
% Uses the unchanged full-order solver kolmogorov_solver.m.
%
%   case            T     snapshots       dt_save   use
%   Re = 20, s1     400   t = 100..400    0.1       training (+ future-time test)
%   Re = 30, s1     400   t = 100..400    0.1       training (+ future-time test)
%   Re = 40, s1     1100  t = 100..1100   0.1       training t<=800, test t>800
%   Re = 40, s2     1100  t = 100..1100   0.5       new-initial-condition test
%
% Re = 25, 35 (unseen Re) and Re = 20, 30 seed 2 (new IC) are taken from the
% existing kf_results/ files (dt_save = 0.5), which is sufficient for testing.
%
% Output: ../rom_data/fine_ReXX_sY.mat  and  ../rom_data/generate_log.txt
% Runtime: about 20 minutes.

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(root);                                   % kolmogorov_solver.m
outdir = fullfile(root, 'rom_data');
if ~exist(outdir, 'dir'), mkdir(outdir); end
logf = fullfile(outdir, 'generate_log.txt');

cases = struct('Re',     {20,    30,    40,     40}, ...
               'seed',   {1,     1,     1,      2}, ...
               'T',      {400,   400,   1100,   1100}, ...
               'dt_save',{0.1,   0.1,   0.1,    0.5});

logmsg(logf, sprintf('START %s', datestr(now)));
for c = 1:numel(cases)
    p = struct('Re', cases(c).Re, 'N', 128, 'n', 4, 'T', cases(c).T, ...
               't_rec', 100, 'dt_save', cases(c).dt_save, 'dt_diag', 0.05, ...
               'seed', cases(c).seed, 'amp', 0.1, 'store', true, 'verbose', true);
    p.dt = 0.01;  if p.Re > 30, p.dt = 0.005; end
    logmsg(logf, sprintf('case %d/%d: Re=%g seed=%d T=%g dt_save=%g ...', ...
        c, numel(cases), p.Re, p.seed, p.T, p.dt_save));
    out = kolmogorov_solver(p);
    W = out.W; t_snap = out.t_snap; params = out.p; x = out.x;  %#ok<NASGU>
    diagn = rmfield(out, {'W', 'w_final'});                     %#ok<NASGU>
    fname = fullfile(outdir, sprintf('fine_Re%02d_s%d.mat', p.Re, p.seed));
    save(fname, 'W', 't_snap', 'params', 'x', 'diagn', '-v7.3');
    logmsg(logf, sprintf('   done in %.1f s, max CFL %.2f, %d snapshots -> %s', ...
        out.runtime, out.cfl_max, size(out.W, 2), fname));
    clear out W diagn
end
logmsg(logf, sprintf('ALL DONE %s', datestr(now)));

function logmsg(f, s)
    fprintf('%s\n', s);
    fid = fopen(f, 'a');  fprintf(fid, '%s\n', s);  fclose(fid);
end
