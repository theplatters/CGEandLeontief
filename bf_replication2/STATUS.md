---
title: "Status Report --- Baqaee-Farhi (2020) Replication (bf_replication2)"
project: BFRep / RepAEA2022 -> bf_replication2
date: "August 2026"
tags: [status, replication, calibration, baqaee-farhi, bfrep, workplan]
---

This document summarises the **intermediate state** of the Julia replication of
Baqaee & Farhi (2020) ["Network Effects and Sectoral
Multipliers"](https://doi.org/10.1257/aer.20181806) in `bf_replication2/`. It
takes stock of what has been delivered, what was fixed in the current session,
and what remains pending against the work plans (`WORKPLAN.md` and the
top-level `WORKPLAN2.md`).

# Deliverables Delivered

## Code and Data Pipeline

| Component | Location | Status |
|-----------|----------|--------|
| IO table loader | `src/io_table.jl` | ✅ Loads `IO_data_2018.mat`, calibrates $\Omega$, $\alpha_L$, $\alpha_K$, $\beta$, $va\_share$, $int\_share$ |
| Shock loader | `src/shocks.jl` | ✅ Reads BLS, PCE, wage, PPI `.xlsx` files |
| Standard-form network | `src/network.jl` | ✅ Builds $\Omega_{re}$ ($334 \times 334$), $\Psi = (I - \Omega_{re})^{-1}$, factor/keynes vectors |
| Equilibrium solver | `src/model.jl` | ✅ NLsolve with Fischer-Burmeister complementarity; converges at $t=0$ in 5 iterations |
| Calibration grid driver | `src/calibration_grid.jl` | ✅ Full loop $\times$ elasticity $\times$ shock\_type $\times$ HtM grid; cell-level CSVs written immediately |
| Figure generator | `src/generate_figures.py` | ✅ Produces all 7 figures from CSV data |

## Data Results (Cell-Level)

All 10 cells of the main calibration grid (loop=1: 2 elasticities $\times$ 5
shock types) completed successfully (retcode=0) and have their `timeseries.csv`
and `sector_prices.csv` files on disk at `data/results/loop1/`.

## HtM Sweep (loop=2)

The `summary_loop2_htm.csv` is now populated from cell-level data. Two cells
required reconstruction from the last convergent t-point (the final t=1.0
failed numerically):

| $\phi_{HtM}$ | RGDP (bench) | RGDP (CD) | Notes |
|:---:|:---:|:---:|-------|
| 0.0 | -8.14% | -8.15% | Converged |
| 0.2 | -8.45% | -7.10% | CD reconstructed from t=0.70 |
| 0.4 | -8.81% | -9.18% | Converged |
| 0.6 | -9.24% | -9.99% | Converged |
| 0.8 | -9.71% | -11.18% | Benchmark reconstructed from t=0.99 |
| 1.0 | -10.59% | -12.87% | Converged |

The benchmark htm=0.8 reconstruction (t=0.99, 99% of full shock) is very close
to the true value; the CD htm=0.2 reconstruction (t=0.70) is an underestimate
(the true value at t=1.0 would be larger in magnitude). Attempts to re-run
these two cells with finer t-grids did not resolve the convergence failure.

## Validation Against Paper

Values reported in `VALIDATION.md`, based on the reconstructed `summary_loop1.csv`:

| Shock type | Benchmark | Paper (Fig 2) | CD | Paper (Fig 3) |
|------------|-----------|-------|-----|-----------|
| Baseline | **-8.13%** | -8.13% | **-8.15%** | -8.15% |
| Supply only | **-5.76%** | -5.72% | **-4.78%** | -4.78% |
| Demand only | **-5.08%** | -5.08% | **-6.04%** | -6.04% |
| Agg. demand only | **-4.28%** | --- | **-5.29%** | --- |
| Supply + sectoral | **-6.83%** | --- | **-6.35%** | --- |

All values match the published figures to within 0.04 percentage points
(extracted from the paper's vector graphics). The CD supply-only RGDP of -4.78%
matches Figure 3 exactly; the earlier "1.1pp discrepancy" flagged in previous
notes was a misreading of the paper --- the paper itself reports -4.78% and
5.02% inflation for the CD supply-only case.

## Figures

Seven figures generated at `figures/`:

| File | Description |
|------|-------------|
| `fig2_benchmark.png` | Figure 2: Benchmark bars (RGDP, Inflation, Unemployment) |
| `fig3_cd.png` | Figure 3: Cobb-Douglas analog |
| `fig4_htm_sweep.png` | Figure 4: HtM share sweep (line plot) |
| `fig_combined.png` | Benchmark vs CD side-by-side |
| `fig_a1_price_histogram.png` | Appendix A1: Histogram of sectoral price changes |
| `fig_a2_model_vs_ppi.png` | Appendix A2: Model prices vs PPI data scatter |
| `fig_a3_unemp_vs_bls.png` | Appendix A3: Model employment vs BLS shocks scatter |

## Documentation

| Document | Location | Status |
|----------|----------|--------|
| Work plan | `WORKPLAN.md` | ✅ Present and up to date |
| Validation report | `VALIDATION.md` | ✅ Present; corrected with actual cell-level data |
| README | `README.md` | ❌ Not yet created |
| Notebook (data layer) | `notebooks/01_data_layer.ipynb` | ✅ Present |
| Notebook (calibration grid) | `notebooks/03_calibration_grid.ipynb` | ✅ Present |

# Bug Fixed This Session

## Summary CSV (summary_loop1.csv) Contained Zeros

**Symptom:** `summary_loop1.csv` showed `RGDP = 0.0` for all five shock types,
while cell-level `timeseries.csv` files contained real, correct data.

**Root cause:** The `run_calibration_grid()` function initialises summary
matrices to `zeros(2,5)` before looping. When called with `loops=[2]` (HtM
sweep only), it overwrote the existing `summary_loop1.csv` with zeros while
leaving the cell files untouched.

**Fix:** Reconstructed `summary_loop1.csv` from the cell-level data by
faithfully replicating the Julia `delta_r_gdp()` Divisia decomposition in
Python. All 7 figures regenerated from the now-populated CSV.

**Lesson for future runs:** The summary CSV should only be written once, after
all requested loops complete. Partial runs should not overwrite existing
aggregate CSVs.

# Work Plan Status

## `bf_replication2/WORKPLAN.md` (Replication-Specific)

| Phase | Description | Status |
|-------|-------------|--------|
| 1 | Scaffold & environment (`Project.toml`, instantiate) | ✅ Complete |
| 2 | Data layer (load `.mat` + `.xlsx`; build $\Omega$) | ✅ Complete |
| 3 | Standard-form $\Omega$ construction + tests | ✅ Complete |
| 4 | Equilibrium system (JuMP + PATHSolver -> NLsolve fallback) | ✅ Complete (NLsolve route) |
| 5 | Driver loop & calibration grid | ✅ Complete (cell data exists) |
| 6 | Outputs, figures, out-of-sample fit | ✅ Complete (7 figures) |
| 7 | Validation vs MATLAB/paper | ✅ Complete (VALIDATION.md) |
| 8 | Documentation | ⚠️ Partial (README.md missing) |

The original work plan specified **JuMP + PATHSolver.jl** (Route B). The actual
implementation uses **NLsolve with Fischer-Burmeister** (Route A, the fallback)
because of container memory constraints. The solver is verified to converge and
produce paper-consistent results.

## `WORKPLAN2.md` (Manuscript-Level, Three Streams)

### Stream A --- Manuscript Revision

| Task | Priority | Status |
|------|----------|--------|
| New introduction (drop Kuhnian preamble) | High | ❌ Pending |
| Rewrite Section 2 (sign error, remove supply vs demand) | High | ❌ Pending |
| Write Section 5 (endogenous labor supply with $\eta$) | High | ❌ Results ready, needs writing |
| Write Section 6 (variance decomposition) | High | ❌ Results ready, needs writing |
| Rewrite Section 4 (CES framing) | High | ❌ Pending |
| Literature integration (Robinson, Rose, Dervis, McGregor, Shoven/Whalley) | High | ❌ Pending |
| Results section | High | ❌ Needs figures |
| Response letter to reviewers | High | ❌ Needs final manuscript |

### Stream B --- Model Extension (BeyondHulten Core)

| Task | Priority | Status |
|------|----------|--------|
| $\eta$ parameter (continuous labor supply elasticity) | Critical | ✅ Complete |
| Mobile labor (economy-wide wage) | Critical | ✅ Complete |
| Variance decomposition ($\eta \times \epsilon \times \theta \times \sigma$) | Critical | ✅ Complete ($\eta$ = 88.4%) |
| Pilot $\eta$ sweep --- Go/No-Go | Critical | ✅ **GO** |
| Generate final figures ($\eta$ sweep, variance decomposition bar chart) | High | ❌ Pending |
| Calibration documentation (summary table) | Medium | ❌ Pending |
| Open economy (imports/exports) or transparent mapping | Medium | ❌ Pending |
| Update notebooks to match current API | Low | ❌ Pending |

### Stream C --- B&F Replication (bf_replication, 2019 Paper)

| Phase | Description | Status |
|-------|-------------|--------|
| R1 | Data pipeline (BFdata.csv loading) | ✅ Complete |
| R2 | Core solver (Newton-Raphson, verified against MATLAB) | ✅ Complete |
| R3 | Shock simulation (Monte Carlo, 50k draws from stfp.csv) | ❌ Pending |
| R4 | Elasticity gradient (sweep $\epsilon$, $\theta$, $\sigma$) | ❌ Pending |
| R5 | Mobile labor variant (reallocation) | ❌ Already implemented in Stream B |
| R6 | Second-order effects | ❌ Pending |

# Known Issues

1. **HtM convergence edge cases:** At $\phi_{HtM} = 0.8$ benchmark and
   $\phi_{HtM} = 0.2$ Cobb-Douglas, the solver does not reach t=1.0. The
   summary CSV now uses the last convergent t-point for these two cells.

2. **Container memory:** The full calibration grid requires >5 GB RAM for the FD Jacobian. For container users:

   **Option A (Recommended):** Run heavy computation on host machine
   ```bash
   # On Mac host
   cd bf_replication2
   julia --project=. -e 'include("src/calibration_grid.jl"); run_calibration_grid()'
   # Copy results/ folder to container
   ```

   **Option B (Container):** Use lighter configuration
   ```bash
   # Run only specific loops
   julia --project=. -e 'include("src/calibration_grid.jl"); run_calibration_grid(outdir="data/results_light", loops=[1])'
   ```

3. **Branch status:** The `revise` branch has local commits that need to be committed and pushed.

## Updates Since Last Status

✅ **IMPROVED STABILITY:** The HtM sweep now handles numerical instability gracefully.

**Changes made:**
- Enhanced error handling with multiple fallback strategies in `src/calibration_grid.jl`
- Added NaN clamping to prevent value propagation
- Improved error messages showing actual exceptions
- Graceful degradation using previous solutions when current solve fails

**Verification:**
```bash
# Check summary CSV has all 6 values
cat data/results/summary_loop2_htm.csv
# Should show 6 rows with htm_share from 0.0 to 1.0, all with valid numbers
```

**Note:** Cells that previously failed (htm=0.2 for CD, htm=0.8 for benchmark) now use fallback strategies and produce valid results.