# Reduced-order models of 2D Kolmogorov flow: POD-Galerkin vs Operator Inference

AS5401 Data-Driven Modelling, final project, IIT Madras.

**Repository:** https://github.com/dvphalgun-bit/Kolmogorov_Flow_ROM

To get a copy of the code:

```bash
git clone https://github.com/dvphalgun-bit/Kolmogorov_Flow_ROM.git
```

This repository holds:

- A full-order model (FOM) of two-dimensional Kolmogorov flow: a pseudo-spectral DNS of the vorticity equation on a 128 × 128 doubly periodic grid, with forcing −4 cos(4y).
- Two quadratic reduced-order models (ROMs), compared at Re = 20, 30 and 40 and at the unseen values Re = 25 and 35:
  - **POD-Galerkin** (intrusive);
  - **Operator Inference, OpInf** (non-intrusive).

## File organisation

```
.
├── kolmogorov_solver.m        FOM: pseudo-spectral solver (IFRK4, 2/3 dealiasing)
├── run_kolmogorov_fom.m       FOM runs (Re = 8-40, two seeds) -> kf_results/
├── kf_results/                FOM figures fig1-fig7 and summary.txt   (case_*.mat not in git)
├── rom/                       ROM study (all scripts below)
├── rom_data/                  fine-sampled FOM data for the ROMs     (not in git, regenerate)
├── rom_results/               ROM results: logs, summary_tables.txt, figures/ (*.mat not in git)
│   └── baseline_standard_POD/ logs of the first run with ordinary POD (for Fig. r1)
└── README.md
```

### Scripts in `rom/`

The list below goes in pipeline order; helper functions are listed under the step that uses them.

| Step | File | Purpose |
|---|---|---|
| 1 | `generate_rom_data.m` | Re-runs the FOM with snapshot spacing 0.1 (Re = 20, 30, 40), writes `rom_data/fine_*.mat` |
| 2 | `kf_ops.m` | FOM operators (forcing, Laplacian, stream function, advection, rhs, E, D) on snapshot matrices |
| 3 | `pod_basis_inv.m` | Translation-invariant POD (x-Fourier × y-POD, cos/sin pairs); used for all final results |
|   | `pod_basis.m` | Standard snapshot POD (baseline only) |
|   | `shift_coeffs.m` | Exact x-shift of POD coefficients (rotation of cos/sin pairs), used for data augmentation |
| 4 | `galerkin_rom.m` | Builds POD-Galerkin operators c, A, H (affine in 1/Re) |
|   | `quad_pairs.m`, `rom_set_Re.m`, `rom_rhs.m` | Quadratic-term ordering, setting Re, ROM right-hand side |
| 5 | `fd_derivative.m` | 4th-order finite-difference time derivatives of the coefficients |
|   | `opinf_fit.m` | Regularised least-squares Operator Inference fit (QR) |
|   | `opinf_select.m` | λ grid search using training data only (stability + short-horizon error + statistics) |
| 6 | `run_rom_case.m` | One Re: basis, both ROMs, λ search, tests (test window, valid time, free run, new IC), persistence baseline |
|   | `rom_integrate.m` | ode45 integration with blow-up detection |
|   | `rom_error.m`, `rom_short_error.m` | Field, projection and coefficient errors |
|   | `rom_ED_pre.m`, `rom_ED.m`, `burst_stats.m` | Energy and dissipation from coefficients; burst statistics |
| 7 | `run_rom_diagnostics.m` | Re = 40: projection error in burst vs quiet periods; energy budget of the resolved modes |
| 8 | `run_rom_param.m` | One parametric ROM (train Re = 20, 30, 40), tested at the unseen Re = 25, 35 |
| – | `run_all_rom.m` | **Master script**: runs steps 3–8 with all settings |
| – | `make_rom_figures.m` | Figures `rom_results/figures/fig_r1`–`fig_r9` and `rom_results/summary_tables.txt` |

## Software dependencies

- **MATLAB** (developed with a recent release; R2020b or newer should work). Only base MATLAB is needed (`fft2`, `ode45`, `qr`, `histogram`); no toolboxes are required.
- About **2 GB of free disk space** for the generated data, and 8 GB of RAM or more.

## How to reproduce the results

Run everything from MATLAB with the repository root as the current folder. Approximate run times on a laptop are given for each step.

1. **FOM runs (about 30 min).** These create `kf_results/case_*.mat`, which the ROM study uses as the new-initial-condition and unseen-Re test data.
   ```matlab
   run_kolmogorov_fom          % set quick_test = true first for a 1-2 min check
   ```
2. **Fine-sampled training data (about 21 min).**
   ```matlab
   run('rom/generate_rom_data.m')
   ```
3. **ROM study (about 15 min).** This creates the `.mat` files and logs in `rom_results/`.
   ```matlab
   clear all; rehash; run('rom/run_all_rom.m')
   ```
   To run only part of the pipeline, edit `steps = {...}` at the top of the script.
4. **Figures and tables.**
   ```matlab
   run('rom/make_rom_figures.m')
   ```
   - The standard-POD curve in `fig_r1` needs `rom_results/baseline_standard_POD/rom_Re20.mat` and `rom_Re40.mat`.
   - To create them, set `base.basis = 'standard'` in `run_all_rom.m`, run the `Re20` and `Re40` steps, move the resulting `rom_Re*.mat` files into that folder, and restore `'invariant'`.

## Main results

- **Re = 20 (relative periodic orbit):**
  - OpInf stays accurate over the full 60-time-unit test horizon for r ≤ 20.
  - POD-Galerkin blows up or overpredicts dissipation by 6–7× for r ≥ 10.
- **Re = 40 (chaotic, with bursts):**
  - OpInf is no better than persistence for r ≤ 61. The resolved dynamics are not closed, and OpInf's fit residual stays at 0.90 or above.
  - POD-Galerkin becomes accurate only at r = 121: valid prediction time 7.6 vs 3.3 for persistence; mean dissipation 0.091 vs 0.090 D_lam; correct bursts.
- **Conditioning of OpInf was not the bottleneck** (κ ≤ 5.7 × 10³). Overfitting, the missing closure, and a POD basis that represents bursts poorly were.

The data files (`*.mat`, about 1.5 GB) are not stored in this repository because of GitHub's size limits. The steps above regenerate them.
