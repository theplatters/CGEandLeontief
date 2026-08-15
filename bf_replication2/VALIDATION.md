---
title: "Validation Report: Baqaee-Farhi (2020) Replication in Julia"
author: "Jakob Kapeller"
project: "beyond-hulten"
date: "August 2026"
tags: [replication, calibration, validation, baqaee-farhi, hulten]
---

This document compares the Julia replication (bf_replication2) of Baqaee and Farhi
["Network Effects and Sectoral Multipliers"](https://doi.org/10.1257/aer.20181806)
(2020, *AER*) against the paper's published results and the original MATLAB codebase
(`Master_file_3.m`).

# Baseline Calibration

The model is calibrated to the 2018 US Input-Output table (66 sectors), using BLS
employment and PPI data. The benchmark calibration uses:

- $\sigma = 1.0$ (final demand elasticity)
- $\epsilon = 0.6$ (domestic supply elasticity)
- $\eta = 0.6$ (import substitution elasticity)
- $\theta_1 = 0.2$ (complementarity in final demand)
- Cobb-Douglas regime: all elasticities set to 1.0

## Real GDP Response (Figure 2/3)

Results from `summary_loop1.csv` versus the paper's published figures (extracted
from Figure 2 and Figure 3 vector graphics):

| Shock type | Julia (bench) | Paper (Fig 2) | Julia (CD) | Paper (Fig 3) |
|------------|:-------------:|:-------------:|:----------:|:-------------:|
| Baseline (all) | **-8.13%** | -8.13% | **-8.15%** | -8.15% |
| Supply only | **-5.76%** | -5.72% | **-4.78%** | -4.78% |
| Demand only | **-5.08%** | -5.08% | **-6.04%** | -6.04% |
| Agg. demand only | **-4.28%** | --- | **-5.29%** | --- |
| Supply + sectoral | **-6.83%** | --- | **-6.35%** | --- |

All values match the paper to within 0.04 percentage points. The CD supply-only
RGDP of -4.78% matches the paper's Figure 3 exactly (text reports "4.8
percent"). The replication therefore matches the published figures within
plotting precision.

## Unemployment and Inflation

Model versus paper (Figure 2 for benchmark, Figure 3 for Cobb-Douglas):

| Shock type | Infl (bench) | Infl paper | Unemp (bench) | Unemp paper | Infl (CD) | Infl paper | Unemp (CD) | Unemp paper |
|------------|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| Baseline | -1.49% | -1.49% | 6.40% | 6.4% | -1.47% | -1.47% | 6.91% | 6.91% |
| Supply only | 6.11% | 6.07% | 1.15% | 1.14% | 5.02% | 5.02% | 0.0% | 0.0% |
| Demand only | -4.66% | -4.66% | 9.47% | 9.47% | -3.69% | -3.69% | 11.27% | 11.27% |

All inflation and unemployment values match the paper's figures to within 0.04
percentage points. The CD supply-only inflation of 5.02% matches Figure 3
exactly. The near-zero unemployment in the CD supply-only case is theoretically
expected: without complementarity, all sectors adjust uniformly and labor
reallocation is minimal.

## Baseline Sector Fit (Appendix)

`baseline_fit.csv` contains model-implied sectoral price changes and unemployment
versus the data. The scatter plot shows a positive correlation for sticky-price
sectors, consistent with the paper's Appendix.

| Sector | Model price $\Delta$ | PPI data | BLS shock |
|--------|---------------------|----------|-----------|
| Air transportation | -18.0% | -9.1% | -39.8% |
| Petroleum & coal | -22.0% | -38.6% | -21.1% |
| Oil & gas extraction | -21.5% | -36.3% | -20.2% |

## HtM Sweep Status (With Stability Fixes)

✅ **IMPROVED STABILITY:** The HtM sweep now handles numerical instability gracefully with multiple fallback strategies.

### Current Results:

| φ_HtM | RGDP (bench) | RGDP (CD) | Status |
|-------|--------------|-----------|--------|
| 0.0   | -8.14%       | -8.15%    | ✅ Converged |
| 0.2   | -8.45%       | -7.10%    | ⚠ Reconstructed (last t) |
| 0.4   | -8.81%       | -9.18%    | ✅ Converged |
| 0.6   | -9.24%       | -9.99%    | ✅ Converged |
| 0.8   | -9.71%       | -11.18%   | ⚠ Reconstructed (last t) |
| 1.0   | -10.59%      | -12.87%   | ✅ Converged |

**Note:** All 6 HtM share values are now **populated**. Two cells did **not** converge at t=1.0 and are **reconstructed** from the last convergent t-point — htm=0.2 (CD, t=0.70) and htm=0.8 (benchmark, t=0.99). The CD htm=0.2 value is an **underestimate** of the true t=1.0 magnitude. The other 4 cells (htm=0.0, 0.4, 0.6, 1.0, both regimes) converged cleanly. Fallback strategies (continuation refinement, relaxed tolerance, previous-solution fallback, NaN clamping) used:
- Continuation refinement with relaxed tolerance
- Previous solution as fallback
- NaN clamping to prevent propagation

### Stability Improvements Added:

The code now includes enhanced error handling:

1. **Multiple fallback strategies** - If solver fails, tries:
   - Continuation refinement with midpoint t-value
   - Relaxed tolerance (1e-8 instead of 1e-10)
   - Uses previous solution as fallback
   - NaN clamping to prevent propagation

2. **Better error reporting** - Shows actual error messages instead of generic warnings

3. **Graceful degradation** - Uses previous values when current solve fails

Run the full HtM sweep:
```bash
julia --project=. -e 'include("src/calibration_grid.jl"); run_calibration_grid()'
```

## Model Implementation Notes

### Solver Architecture

The equilibrium is solved using NLsolve with a Fischer-Burmeister reformulation of the complementarity conditions. This follows "Route A" (as opposed to JuMP/PATHSolver), keeping dependencies minimal.

Key design choices:

- **Continuation method:** The shock is scaled by $t \in [0, 1]$ and the solver walks along the grid, using the previous solution as the next initial guess.
- **Continuation refinement:** If the solver fails at a grid point, a midpoint $t_{mid}$ is inserted and solved first.
- **Trunc_A persistence:** Following the MATLAB driver, the labor-demand Trunc_A matrix persists across all shock_type cells and determines demand-constrained sectors.

## Corrections Applied During Translation

| Issue | MATLAB behaviour | Initial Julia bug | Fix |
|-------|------------------|-------------------|-----|
| Sign of shock | `1 - t .* shock_A` | Typo `1 - t .* shock_A` sign | Changed to `1 + t .* shock_A` for absolute sign |
| $\lambda$ initialisation | Uses calibrated initial point | Self-referential $\lambda$ assignment | Fixed to `factor_clearing0` calculation |
| B renormalisation | Row vector $\Omega \cdot B$ | Transpose broadcasting | Removed transpose |
| $\chi$ parameter | Not used (blank parameter) | `rand()` seeded | Fixed to explicit `0.0` |
| **Trunc_A dimensions** | Persists across all cells | Only checked n_t dimension | Now checks both N and n_t dimensions |

## Dependency Status

```{.text}
Julia project: bf_replication2/src/
  model.jl          --- MCP model, solver, calibration
  network.jl         --- Network matrices (IO data loading)
  test_model.jl      --- Unit tests (data layer + equilibrium)
  calibration_grid.jl--- Main driver (Phase 5)
  generate_figures.py--- Figure generation (Python/matplotlib)

Notebooks:
  01_data_layer.ipynb    --- Data loading verification
  02_equilibrium.ipynb   --- Solver walkthrough + CSV export
  03_calibration_grid.ipynb --- Full calibration driver
  04_figures.ipynb        --- Figure generation
```
JuMP/PATHSolver), keeping dependencies minimal.

Key design choices:

- **Continuation method**: The shock is scaled by $t \in [0, 1]$ and the solver
  walks along the grid, using the previous solution as the next initial guess.
- **Continuation refinement**: If the solver fails at a grid point, a midpoint
  $t_{mid}$ is inserted and solved first.
- **Trunc_A persistence**: Following the MATLAB driver, the labor-demand Trunc_A
  matrix persists across all shock_type cells and determines demand-constrained
  sectors.

# Verification Checklist

- [x] Data layer loads correctly (all 5 checks pass)
- [x] Equilibrium solver matches initial calibration at $t=0$
- [x] Continuation converges for benchmark configuration
- [x] Baseline RGDP change: **-8.13%** (paper: -8.13%)
- [x] Supply-only: **-5.76%** (paper: -5.72%)
- [x] Demand-only: **-5.08%** (paper: -5.08%)
- [x] CD baseline: **-8.15%** (paper: -8.15%)
- [x] CD supply-only: **-4.78%** (paper: -4.78%)
- [x] CD demand-only: **-6.04%** (paper: -6.04%)
- [x] Loop-1 summary CSV populated (reconstructed from cell-level data)
- [x] HtM sweep summary CSV fully populated; 2 of 6 cells (htm=0.2 CD, htm=0.8 benchmark) reconstructed from last convergent t-point (NOT converged at t=1.0)
- [x] All figures generated from CSV data
- [ ] Fully automated `make all` pipeline

# Known Issues

1. **HtM convergence edge cases (reconstructed, not converged)**: Two cells — htm=0.2 (CD, t=0.70) and htm=0.8 (benchmark, t=0.99) — do not reach t=1.0. Both are **reconstructed** from the last convergent t-point and flagged in the tables; the CD htm=0.2 value is an underestimate. The earlier Trunc_A dimension-mismatch bounds error is patched.
2. **Reconstruction is not convergence**: The two reconstructed cells above are marked "⚠ Reconstructed", not "✅ Converged". All other HtM shares (0.0, 0.4, 0.6, 1.0) converge cleanly for both regimes.
3. **Container memory**: The full grid requires >5 GB RAM. Run on the host Mac.

# References

Baqaee, David Rezza, and Emmanuel Farhi. "Network Effects and Sectoral
Multipliers." *American Economic Review* 110, no. 12 (2020): 3712-58.
`https://doi.org/10.1257/aer.20181806`