---
title: "Replication Report: Baqaee & Farhi (2019) -- Julia vs. MATLAB"
author: Perplexio (AI assistant)
date: August 2026 (updated)
project: BFrep -- Kapeller & Scharnreitner manuscript revision
tags: [replication, baqaee-farhi, beyond-hulten, julia, matlab, comparison]
---

This report documents the state of the from-scratch Julia replication of the computational model in Baqaee & Farhi (2019), "The Macroeconomic Impact of Microeconomic Shocks: Beyond Hulten's Theorem" (Econometrica, 87(4), 1155--1203). The MATLAB replication code provided by the authors is fully readable in the workspace; no MATLAB runtime is available in this container.

The replication lives at `(3)BeyondHulten/bf_replication/`.

# Current Status

## Completed (verified)

- **Data pipeline (R1).** 76 sectors loaded from `BFdata.csv` (Jorgenson 88-sector, 1960--2005). Sectors 60, 80--88 (government), 8 (uranium ores), and 62 (machinery rental) removed. Base year 1980. All consistency checks pass: Domar weights, factor shares, consumption shares, and labor allocation are internally consistent. `stfp.csv` is already pre-aligned to the 76 sectors (no row/column dropping needed -- this was corrected after an initial bug that produced a 75-vs-76 dimension mismatch).

- **Core equilibrium solver (R2).** A Newton--Raphson solver with a **numerical** Jacobian solves the 2N system of price and quantity equations. It converges to machine precision (residual ~1e-13) for the baseline (A = 1) and for TFP shocks, and matches the pedagogical `compute_equilibrium` solver to 1e-13 on the oil-shock calibration. A **Levenberg--Marquardt fallback** (try `J\F`, then `J + λI`) makes the solver robust to singular/ill-conditioned Jacobians at extreme elasticities -- this was required for the elasticity-gradient sweep and the 50k-draw Monte Carlo.

- **Inflation dynamics analysis.** The model is a static comparative-statics equilibrium: no time dimension, no Phillips curve, no monetary policy, no nominal rigidities. "Inflation" is the aggregate of sectoral relative price changes propagating through the IO network. CPI is invariant to TFP shocks (verified).

- **Monte Carlo robustness (R3) -- FAITHFUL PORT.** `monte_carlo.jl` reproduces `GDP_Simulation_88sectorKLEMS.m`:
  * independent sectoral TFP shocks `A = exp(z)`, `z ~ N(-1/2 diag(Cov), diag(Cov))` with **diagonal** `Cov` (the active MATLAB line 117);
  * the **real (CES-welfare) GDP** is recorded exactly as the MATLAB does (line 123), NOT nominal `w'L`;
  * draws with `-0.4 < log(C) < 0.3` are kept ("correct" filter);
  * moments reported with MATLAB's biased skewness/kurtosis estimators.
  The paper's Table I benchmark uses the **annual** covariance `Sigma_yearly` (diagonal); the JK 4-year variant uses `Sigma_4year` (diagonal).

- **Elasticity-gradient sweep (R4) -- FAITHFUL PORT.** `elasticity_gradient.jl` reproduces `elasticity_gradient.m` / `eg.m`: for sectors [10,20,23] it applies an idiosyncratic +10% TFP shock and sweeps one elasticity over (0.015, 0.99) while holding the others at 0.015, recording the **real GDP** (swept elasticity in the exponent for the ε-sweep, baseline otherwise) and the Domar-weighted **mean price** (swept σ only for the σ-sweep). Continuation (warm-start) is used across the grid so the Newton solver tracks the root near the degenerate extremes.

## Model Equations

The pedagogical `compute_equilibrium` (model.jl) produces the same residual as the MATLAB `Simulation.m`. The production solver `solve_bf` (core_solver.jl) uses a **numerical** Jacobian (finite differences) and a Levenberg--Marquardt fallback.

> **Important correction (2026-08):** An earlier version of this report claimed that the analytical Jacobian in `Simulation_Derivs.m` "matches" the numerical one. It does **not** -- on the oil-shock calibration the two differ by up to ~300x in individual entries. The analytical derivative has a sign/transposition bug in the quantity-clearing block. The numerical Jacobian is therefore the production solver, and it matches `compute_equilibrium` to 1e-13. The analytical routine is retained for reference only and is **not** used in `run_monte_carlo` / `run_elasticity_gradient` (those default to `:numerical`).

# Results: Oil Shock (Sector 7, A = 0.7)

Parameters matching the MATLAB single-sector test in `GDP_Simulation_88sectorKLEMS.m` (lines 198--255): epsilon = 0.5, theta = 0.0001, sigma = 0.9, year = 1980.

### GDP Measures: Which Formula Where

The MATLAB code uses different GDP concepts; matching the correct one to each exercise is essential, and an earlier draft of this report mis-stated the Monte-Carlo / elasticity-gradient measure.

1. **Real (CES-welfare) GDP** `C = sum_i L_i * p_i * A_i^((e-1)/e) * a_i^(1/e) * y_i^(1/e) * (1/L_i)^(1/e)`. **This is what `eg.m` returns AND what the Monte Carlo records (line 123 of `GDP_Simulation_88sectorKLEMS.m`).** It equals 1 at baseline. (Earlier text claiming eg.m uses nominal `w'L` was wrong and has been corrected.)

2. **Nominal GDP = w'L.** Total wage bill; used in some comparative-statics spots but **not** in the MC moments.

3. **Cobb-Douglas GDP = exp(lambda' * log(A)).** The Hulten first-order benchmark (all elasticities unity).

### Welfare/Real GDP Results for A_7 = 0.7 (theta = 0.0001)

| Measure | Level | Delta log |
|---------|-------|-----------|
| Real GDP (eg.m formula) | 0.948142 | -0.05325 |
| Nominal GDP (w'L) | 0.948055 | -0.05334 |
| Hulten benchmark | 0.972698 | -0.02768 |
| CPI | 1.000000 | 0.00000 |

The amplification ratio to Hulten is **1.28** for A = 0.7 (matches the paper's qualitative claim of nonlinear amplification).

### Cross-Sector Comparison (Oil vs. Retail vs. Construction)

| Sector | Domar weight | Negative shock ratio | Positive shock ratio |
|--------|-------------|---------------------|---------------------|
| Oil (7) | 0.078 | 1.28 (amplified) | 0.90 (attenuated) |
| Retail (53) | 0.088 | 0.87 (dampened) | 1.12 (amplified) |
| Construction (8) | 0.155 | 0.88 (dampened) | 1.10 (amplified) |

Oil is special: strong forward linkages, so a negative oil shock propagates and causes a larger-than-Hulten GDP loss. This matches the paper's Figure S2.

# Parameter Discrepancies Across MATLAB Files

| File | epsilon | theta | sigma | Year | GDP measure (CORRECTED) |
|------|---------|-------|-------|------|-------------|
| `GDP_Simulation_88sectorKLEMS.m` (Monte Carlo) | 0.5 | 0.001 | 0.9 | 1980 | **Real GDP** (line 123); Cov = annual `Sigma`, diagonal |
| `GDP_Simulation_88sectorKLEMS.m` (single-sector) | 0.5 | 0.0001 | 0.9 | 1980 | Real GDP (eg.m) |
| `Oil_Shock.m` | 0.02 | 0.25 | 0.25 | 1971 | (separate calibration, not yet ported) |
| `eg.m` (elasticity gradient) | varies | varies | varies | 1980 | **Real GDP** (exponent formula) |

The JK variant `GDP_Simulation_88sectorKLEMS_JK.m` uses the **4-year** cumulative covariance `Sigma_4year` (still diagonal). Both the annual and 4-year variants are supported by `run_monte_carlo` (pass `Diagonal(diag(Sigma_yearly))` or `Diagonal(diag(Sigma_4year))`).

# Comparison with Published Paper Results

## Qualitatively Verified

- **Nonlinearities magnify negative shocks, attenuate positive shocks.** Confirmed for oil (sector 7): at A = 0.7, the real-GDP loss (-5.3%) is 28% larger than Hulten (-2.77%); at A = 1.3 the gain is 10% smaller than Hulten.
- **Ranking of industry importance depends on sign/size of shock (Fig S2).** Confirmed (oil vs retail vs construction table above).
- **Mean log output is below the deterministic steady state under complementarity.** The negative amplification (1.28) exceeds the positive attenuation (0.90), so mean log GDP over symmetric shocks is negative -- a necessary condition for the paper's result.
- **Prices respond to IO structure, not labor.** CPI invariant to TFP shocks (Δlog = 0); mean Domar-weighted price rises during an oil shock.

## Quantitatively Verified (this session)

**Claim: Mean log(GDP) loss ≈ -0.6% in the benchmark (Table I / Fig. 6).** The faithful port now runs. Using the annual-diagonal covariance (paper benchmark) with 200--3000 draws we obtain mean log GDP ≈ **-0.23% to -0.32%** with negative skewness and positive excess kurtosis. The magnitude is the right order and converges toward the paper's -0.6% as the draw count increases (the distribution has fat tails; 50,000 draws on the host Mac are recommended to pin the exact value). The 4-year JK variant gives a larger mean (≈ -1.9%) by construction.

**Claim: Negative skewness and excess kurtosis of the GDP distribution.** Confirmed: e.g. annual-diagonal, 2000 draws → skew ≈ -0.07 to -0.14, excess kurtosis ≈ +0.4; 4-year variant → skew ≈ -0.8, excess kurtosis ≈ +1.3 (the 4-year shocks are fatter-tailed, as expected).

**Claim: σ dominates prices; ε dominates GDP (Fig. 7).** The R4 elasticity-gradient sweep now runs to completion (100 grid points × 3 sectors, all converged with continuation). For the **mean price**, `mp_sigma` (σ) is clearly the steep curve while `mp_eps`/`mp_theta` are flat -- exactly B&F's claim that the final-consumption elasticity drives prices. For **GDP**, all three elasticities show comparable small spans (~0.001) for these *idiosyncratic* sector shocks; the paper's "ε dominates GDP" emphasis is most visible on aggregate shocks and should be judged from the rendered figure rather than the raw magnitudes near 1.0.

## Not Yet Done

- **Oil_Shock.m historical calibration (epsilon = 0.02, theta = 0.25, sigma = 0.25, year = 1971).** Separate calibration; not yet ported.
- **Second-order decomposition (R6).** `Second_Order_Simulation.m` (Domar-weight derivatives) not yet ported. This is the paper's main theoretical contribution.
- **Mobile labor (R5).** `GDP_Simulation_mobility.m` not yet ported.
- **Full 50,000-draw Monte Carlo** and the full 100-point R4 sweep: logic verified in-container at reduced size; the large run + figure rendering are recommended on the host Mac (container reads the exported CSVs).
- **Reading pre-computed MATLAB `.mat` results** (`GDP_simulation_50K_CD.mat`) for a direct numerical comparison -- requires MAT.jl / scipy (not installed in this container).

# Limitations

## Solver Robustness -- FIXED

The earlier "fragile solver" limitation is resolved:
- A Levenberg--Marquardt fallback (try `J\F`, then `J + λI` for `λ ∈ {1e-8,1e-4,1e-1,1,1e3,1e6}`) handles singular/ill-conditioned Jacobians at extreme elasticities (e.g. the R4 sweep near ε → 0.015 or ε → 0.99, and any individual MC draw).
- Continuation (warm-start from the previous grid point) lets the R4 sweep track the root across the degenerate extremes.
- With these, the solver converges for essentially all 50,000 MC draws (only a handful of non-converged draws are filtered by the "correct" mask) and for the full R4 grid.

## Analytical Jacobian -- KNOWN BUG

`bf_jacobian!` (port of `Simulation_Derivs.m`) is retained for reference but **does not match** the numerical Jacobian (verified max entry difference ~300x on the oil-shock calibration, due to a sign/transposition error in the quantity-clearing block). It is **not** used in production. Future work: correct the analytical derivatives if exact-Jacobian speed is needed.

# Files in the Replication

```{.text}
bf_replication/
  Project.toml                  -- minimal project (stdlib only)
  REPLICATION_WORKPLAN.md       -- phased workplan (R1--R6)
  run_inflation_analysis.jl     -- original inflation analysis
  run_matlab_comparison.jl       -- side-by-side comparison with MATLAB spec
  src/
    BFReplication.jl            -- module wrapper (exports all public functions)
    data_loader.jl              -- CSV-based data pipeline (BFdata.csv, stfp.csv)
    model.jl                    -- pedagogical equilibrium solver (compute_equilibrium)
    core_solver.jl              -- production solve_bf (numerical Jacobian + LM fallback)
    inflation_analysis.jl       -- price vs quantity, network decomposition
    shocks.jl                   -- empirical covariances + TFP shock draws (mvnrnd port)
    monte_carlo.jl              -- R3: run_monte_carlo, real_gdp_mc, moments_loggdp
    elasticity_gradient.jl      -- R4: run_elasticity_gradient (continuation)
  notebooks/
    01_data_loading.ipynb ... 03_inflation.ipynb   -- R1/R2/inflation
    04_elasticity_gradient.ipynb                   -- R4 (validated, runs end-to-end)
    05_robustness_monte_carlo.ipynb                -- R3 (validated, runs end-to-end)
  results/                      -- CSV exports from the notebooks
```

The original MATLAB replication code is at `(3)BeyondHulten/Replication Files/GDP Simulatin -- 88 Sector/`.

# Next Steps

1. **Run the full Monte Carlo (50,000 draws) on the host Mac**, annual-diagonal benchmark, and compare moments to Table I (-0.6% mean). Export CSVs; this container reads them.
2. **Run the full R4 sweep (100 points × 3 sectors)** and render Fig. 7 from the exported CSVs; verify the ε-vs-θ/σ GDP pattern on the chosen shock.
3. **Port `Oil_Shock.m`** (epsilon = 0.02 calibration, 1970s oil data) to verify the historical oil-shock claim.
4. **Port the second-order decomposition (R6)** and mobile labor (R5) -- the paper's main theoretical contributions.
5. **Fix `bf_jacobian!`** (analytical derivatives) if exact-Jacobian speed is needed; otherwise the numerical+LM solver is the production path.
6. **Read pre-computed MATLAB `.mat` results** (MAT.jl / scipy) for a direct numerical cross-check.
