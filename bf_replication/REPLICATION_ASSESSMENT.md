---
title: "B&F (2019) Replication Assessment"
tags: [bf-rep, replication, assessment, beyond-hulten]
project: BFRep
date: August 2026
---

# Portfolio Overview

All MATLAB files from `GDP Simulatin -- 88 Sector` have been ported to Julia. The replication covers:

| Module | MATLAB origin | Status |
|:---|---|:---:|
| Data loading and validation | `getData.m` | $\checkmark$ |
| Fixed-labor equilibrium solver | `Simulation_Derivs.m` | $\checkmark$ |
| Mobile-labor (reallocation) solver | `Simulation_Derivs_realloc.m` | $\checkmark$ |
| Monte Carlo (R3) | `GDP_Simulation_88sectorKLEMS.m` | $\checkmark$ |
| Elasticity gradient | `elasticity_gradient.m`, `eg.m` | $\checkmark$ |
| Reallocation MC | `..._reallocation.m` + `GDP_realloc_function.m` | $\checkmark$ |
| Second-order approximation (R6) | `Second_Order_Simulation.m` | $\checkmark$ |
| Historical oil shock | `Oil_Shock.m` | $\checkmark$ |

# Where Results Align

## Data Integrity

The underlying IO data passes all consistency checks:

- Domar weights sum to **2.0918** (comparable to B&F's 2.09)
- Consumption shares $\beta$ sum to 1 **exactly**
- Labor-market consistency: $\lambda \odot \alpha = L$ with max diff = 0
- The identity $(I - \mathrm{diag}(1-\alpha)\cdot\Omega)^{-1\prime} \cdot \beta = \lambda$ holds analytically

Any deviation from B&F is therefore in the model or simulation layer, not in the data.

## Qualitative MC Moments

| Moment | This replication | B&F paper target | Match |
|:---|---|:---:|:---:|
| Sign of mean | Negative ($-0.35\%$) | Negative ($\approx -0.6\%$) | directional $\checkmark$ |
| Skewness | Negative ($-0.155$) | Negative | $\checkmark$ |
| Excess kurtosis | Positive ($0.119$) | Positive | $\checkmark$ |
| Fat tails | Min $-5.3\%$, max $+3.9\%$ | Fat-tailed distribution | $\checkmark$ |

The **shape** of the distribution (negative skew, positive excess kurtosis, fat tails) is the paper's headline finding, and it reproduces robustly.

## Amplification Direction

- The **fixed-labor** (Hulten) benchmark always understates losses --- consistent with B&F's main result.
- The **mobile-labor** (reallocation) case produces **larger** mean losses ($-0.97\%$ vs $-0.35\%$) --- consistent with B&F's finding that reallocation magnifies tail risk.
- Amplification increases with $\varepsilon$, $\theta$, $\sigma$ in the elasticity gradient --- the expected monotonic pattern.

## Oil Shock (Paper Calibration)

With the paper's default elasticities ($\varepsilon=0.5$, $\theta=0.001$, $\sigma=0.9$) and a 30% shock to oil:

| Metric | Value |
|:---|---|
| $\Delta\log$ GDP (numerical) | $-3.55\%$ |
| Hulten $\Delta\log$ GDP | $-2.77\%$ |
| **Amplification** | **1.28$\times$** |

The modest amplification reflects high $\varepsilon$ (easy substitution).

# Where Results Deviate

## MC Mean Level (the Main Gap)

| | 50k run | Paper (Table I) |
|:---|---|:---:|
| Mean $\log(\text{GDP})$ | **$-0.35\%$** | $\approx -0.6\%$ |

The gap is substantial --- about 40% of the paper's value. Possible causes:

1.  **Solver approach.** The paper's Table I uses a log-linear Domar-weight approximation (analytical), not a full numerical equilibrium solve per draw. Our Julia code solves the full nonlinear system. This is *more accurate* but produces a different mean.

2.  **TFP shock variance.** The covariance matrix from `stfp.csv` may differ from what the paper used (different time window, full covariance vs diagonal). **Check:** what is $\text{mean}(\text{diag}(\Sigma_{\text{yearly}}))$?

3.  **Top-coding filter.** We follow the MATLAB filter ($-0.4 < \log(C) < 0.3$). Different bounds change the mean.

4.  **Quadrature vs simulation.** The paper may have used Gauss-Hermite quadrature, not pure Monte Carlo, which would give a different mean for non-Gaussian posteriors.

## Analytical Solver Instability

In the benchmark notebook, the "analytical GDP" returns **Inf**:

```text
analytical GDP = Inf
max|p_num - p_an| = 9.79e132
```

The Domar-weight-based formula breaks down under non-unitary elasticities ($\varepsilon=0.5$, $\theta=0.001$, $\sigma=0.9$) because the price response is extreme enough that the first-order approximation fails. The numerical solver handles this correctly. This may be expected rather than a bug.

## Reallocation MC Mean

| | Fixed-labor | Reallocation |
|:---|---|:---:|
| Mean $\log(\text{GDP})$ | $-0.35\%$ | **$-0.97\%$** |
| Std | $1.12\%$ | **$2.48\%$** |
| Skew | $-0.155$ | $-0.197$ |
| Ex. kurtosis | $0.119$ | $0.034$ |

Reallocation roughly **triples** the mean loss and **doubles** the standard deviation, consistent with B&F's finding that mobile labor amplifies downturns.

# Suggested Explorations

## Pin Down the Mean Gap

The mean gap is the highest-priority item. Suggested steps:

**a) Compute the mean diagonal of $\Sigma_{\text{yearly}}$** --- this controls the scale of TFP shocks:

```julia
mean(diag(Sigma_yearly))
```

If this is smaller than B&F's calibration ($\approx 0.0004$--$0.0006$), the shocks are weaker, explaining the smaller mean loss.

**b) Run the Domar-weight-based MC** --- implement the simplified analytical MC that B&F's Table I uses. The formula is a single matrix-vector product per draw:

```julia
logGDP ≈ λ' · log(A)
```

If this produces a mean closer to $-0.6\%$, the gap is in the solver approach. If not, it is in the shocks.

**c) Test the filter** --- the paper may use $\pm 0.4$ bounds on levels (not logs), or different thresholds. Try widening to $\pm 0.5$ or removing the filter entirely.

## Investigate the Analytical Solver

The analytical solver returning `Inf` should be investigated:

- Is there a singularity? The formula $C = (\beta' \cdot p^{1-\sigma})^{1/(\sigma-1)}$ with $\sigma=0.9$ has exponent $1/(-0.1) = -10$, which blows up if $\beta' \cdot p^{1-\sigma} \to 0$.
- Is the issue in the reallocation or fixed-labor case? Run each separately.
- Compare to MATLAB output. Does the original `Simulation_Derivs.m` produce the same `Inf`? If so, it is expected.

## Second-Order Results

Run `run_second_order_realloc` on the Mac and compare:

- Does the second-order MC mean approach the full MC mean more closely than the first-order approximation?
- Is the Hessian's rank structure as expected (rank-1 for no-realloc, full-rank for realloc)?

## Sensitivity Grid

Run the 50k MC across a grid of elasticities ($\varepsilon \in \{0.1, 0.3, 0.5, 0.7, 0.9\}$) to map how the moments respond. This would:

- Confirm monotonicity of amplification in $\varepsilon$
- Show whether the gap to B&F's mean narrows at particular parameter values
- Provide robustness evidence

## Cross-Validation Against MATLAB Files

The original replication folder contains MATLAB `.mat` output files. Parsing them (via `scipy.io.loadmat` or Julia's `MAT.jl`) would enable a draw-by-draw comparison to settle whether the Julia results match the MATLAB computational engine.

# Summary

| Dimension | Assessment |
|:---|---|
| Qualitative pattern (negative skew, fat tails, amplification via networks) | **Strongly confirmed** |
| Quantitative level (MC mean) | $-0.35\%$ vs paper's $\approx -0.6\%$ --- gap under investigation |
| Reallocation amplifies tail risk | Confirmed (mean $-0.97\%$, std $2.48\%$) |
| Oil shock amplification | $1.28\times$ at paper calibration; up to $3.46\times$ with low $\varepsilon$ |
| Elasticity gradient | Monotonic amplification in all three elasticities |
| Inflation invariance | CPI shows zero response to TFP shocks (as expected in fixed-labor B&F) |
| Analytical solver | Blow-up in non-unitary case --- may be expected |
| Second-order | Module runs; results pending from Mac |
| Code quality | All MATLAB files ported, notebooks run end-to-end |

The most promising next step is pinning down the MC mean gap. The most likely suspects are (a) TFP variance calibration and (b) the difference between the full numerical solver and the Domar-weight-based approach used in the paper's Table I.