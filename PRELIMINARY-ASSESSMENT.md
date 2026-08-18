---
title: "Preliminary Assessment: BeyondHulten Revision Plan and Code Plausibility"
author: "calculato (AI research assistant)"
date: "2026-08-18"
tags: [revision, assessment, beyondhulten, code-review]
project: "BFRep (3)BeyondHulten"
---

This document evaluates the revision plan laid out in `ROADMAP.md` against the
two referee reports in `(2)Reviews`, and assesses the plausibility of the code
that directly supports the revision. The two replications in `bf_replication/`
and `bf_replication2/` are explicitly **out of scope**. All empirical claims
below were reproduced by running the headless test harness under Julia 1.12.6;
`tests/minimal_test/` is byte-identical to `src/` for every core file.

# The ROADMAP Versus the Workplans

The `ROADMAP.md` is a **sound, review-aware plan**. It correctly reframes the
paper around "closure versus technology" (answering Rev2's objection that "the
bridge already exists -- endogenize labor supply"), and it lists non-negotiable
validation gates in sections 4.1 to 4.4: accounting consistency, a
Tornqvist/Divisia real-GDP index in *every* variant, an explicit financing
closure for the investment shock, and a clear taxonomy of labor-market closures
(voluntary supply, involuntary unemployment, mobility, perfectly elastic supply,
exogenous endowment expansion). It also states, correctly, that the manuscript
must not report a variance share until the underlying model has passed those
gates.

The `WORKPLAN.md`, `WORKPLAN2.md`, and `OUTLOOK_SUMMARY.md` documents
**contradict the ROADMAP**. They declare the four critical coding tasks
"complete and verified", the Go/No-Go criterion as **"GO"**, and report
*"$\eta$ = 88.4% of explained variance, 26/71 sectors with CV > 0.01"*. The
current code does not support either claim:

- The code produces **$\eta \approx 100\%$**, not 88.4%.
- That dominance is an **artifact of a broken model**, not evidence of the
  economic bridge.
- The "GO" criterion was satisfied precisely because $\eta$ collapses GDP so
  violently that sectoral variation appears large.

The planning documents are therefore ahead of the code in a dangerous way: they
certify success the code does not deliver.

# Code Plausibility Assessment

## What Works

- The base Baqaee-Farhi machinery (German IO calibration, CES cost functions,
  `conditional_input_shares`, `generate_data`) is plausible. The base CES model
  gives a baseline `real_gdp = 1.0` and a sensible **+3.26%** response to the
  standard construction demand shock.
- The modular structure (`ces.jl`, `leontief.jl`, `cobbdouglas.jl`,
  `interface.jl`, `solution.jl`) and the headless `tests/minimal_test/` harness
  (stubbed GLMakie/XLSX, identical to `src/`) are good engineering practice and
  made this verification possible.

## The Mobile-Labor Model Is Broken

The entire revision rests on a smooth "bridge": $\eta = 0$ (full-employment CGE)
should yield a GDP effect of approximately 0%, and $\eta \to \infty$ (IO
multiplier) should yield approximately +1.5%. The code produces the **opposite**.

- **Inverted response.** Raising $\eta$ *reduces* real GDP. At $\eta = 10$ the
  economy collapses to roughly 12% of baseline. Higher labor elasticity should
  make more labor available and raise GDP, not less.
- **Root cause (economic).** Labor supply is specified as
  $L = \bar{L} \cdot w^\eta$, which measures $w$ against 1.0. Under the demand
  shock the equilibrium wage settles at $w \approx 0.94 < 1$, so $w^\eta$ falls
  as $\eta$ rises, yielding *less* labor. A correct closure must anchor at the
  *baseline* wage, $L = \bar{L} \cdot (w/w_0)^\eta$, and must ensure the wage
  *rises* under excess labor demand. The wage currently *falls* even in the base
  CES model (to 0.989) under a positive demand shock -- a sign that the
  final-demand / income-expenditure closure is misspecified (the unnormalized
  `consumption_share .* demand_shock` in `final_demand` breaks the
  income-equals-expenditure identity).
- **Spurious equilibria / non-robustness.** At $\eta = 10$ the solver returns
  $w = 0.14$ and, depending on the starting vector, GDP of either 0.12 *or*
  1.13 for the same $\eta$. The system does not deliver a unique, economically
  meaningful equilibrium.
- **Cobb-Douglas limit fails.** `variance_decomposition` throws a
  `DomainError` (negative base to a non-integer power) at standard grid points
  with $\varepsilon = 0.99$ -- the model cannot even solve at the Cobb-Douglas
  case the paper discusses. Those points are silently dropped as `NaN`.

## Variance Decomposition Is Methodologically Fragile

- It computes a **one-way marginal** sum of squares $\text{SS}(f)/\text{SS}(\text{total})$,
  then the headline "% Share" **renormalizes** the four factors to sum to 100%.
  So "$\eta$ = 88.4%" really means *"$\eta$'s share of the sum of one-way factor
  variances"*, not an absolute share of output variance. The docstring's claim
  that the decomposition "sums to $\le 1$; residual = interactions" is wrong for
  this estimator: one-way marginal SS double-counts interactions, so the reported
  "Residual" can go negative.
- Because $\eta$ drives a degenerate GDP collapse, it trivially "dominates" --
  the decomposition **confirms the bug** rather than the bridge. The result is
  also grid-design-dependent (an $\eta$-range of $[0, 10]$ versus $[0, 0.5]$
  would flip the shares); that sensitivity is not reported.
- **No persisted results.** The decomposition numbers live only in workplan
  docstrings; no CSV is saved in the project, contrary to the stated
  reproducibility preference.

## ROADMAP-Compliance Gaps in the Code

- **Section 4.2 (real-GDP index).** `solve(Model{MobileLaborCES})` returns
  `real_gdp` via a **Laspeyres** index, not the mandated shared
  `tornqvist_quantity_index` used by the base CES model. The two variants use
  inconsistent real-GDP measures.
- **Section 4.2 (numeraire).** CPI = 1 is imposed by *dropping a price /
  zero-profit equation* and adding the CPI constraint. The ROADMAP explicitly
  warns against replacing a zero-profit equation merely to balance the equation
  count; this needs explicit verification, and empirically the system yields
  spurious equilibria.
- **Section 4.3 (financing).** The investment shock is an **unfinanced** demand
  add-on. There is no lump-sum-tax, expenditure-switching, or foreign-financing
  closure. Until one is added, the CGE experiment is not admissible by the
  ROADMAP's own rule.
- **Section 4.2 (homogeneity).** Not demonstrably verified.

# Recommended Next Steps

1. **Do not report 88.4% / "GO".** Re-derive after the model is fixed and the
   bridge runs in the correct direction.
2. **Fix the mobile-labor closure.** Anchor $L = \bar{L} \cdot (w/w_0)^\eta$ at
   the baseline wage; verify the wage *rises* under the demand shock; confirm
   $\eta = 0 \to \approx 0\%$ and $\eta \to \infty \to \approx +1.5\%$.
   Diagnose the income-expenditure imbalance in `final_demand`.
3. **Guard the Cobb-Douglas limit** (complex-safe powers) so
   $\varepsilon = 0.99$ solves.
4. **Use Tornqvist consistently** across all variants (section 4.2).
5. **Add a financing closure** (section 4.3) -- at minimum a budget-neutral
   preference reallocation or lump-sum tax.
6. **Report the variance decomposition** as absolute shares with the true
   residual, the exact grid, and a sensitivity to the $\eta$-range; save results
   to CSV.
7. **Verify homogeneity** (section 4.2).

# Evidence and Reproduction

The following were verified by executing the harness under Julia 1.12.6. The
shell command used was:

```{.bash}
/root/.julia/linux-julia-1.12.6/bin/julia --project=tests/minimal_test \
  tests/minimal_test/run_test.jl
```

Key results:

| Check | Result |
|-------|--------|
| Base CES model (full labor slack) | Baseline `real_gdp = 1.0`; +81% construction demand -> `real_gdp = +3.26%` |
| Mobile-labor $\eta$-sweep (current code) | $\eta=0\to0.982$, $\eta=1\to0.923$, $\eta=2\to0.867$, $\eta=5\to0.720$, **$\eta=10\to0.123$ (-86%)**, $\eta=50\to0.123$ |
| Mobile-labor solver at $\eta=10$ | Returns $w = 0.14$; result is initialization-dependent (GDP 0.12 or 1.13) |
| `variance_decomposition` (current code) | **$\eta = 100.0\%$** (not 88.4%); `DomainError` at $\varepsilon = 0.99$ |
| `pilot_eta_sweep` | Prints "GO" -- driven by the broken model |

Three reproducible diagnostic scripts remain in `tests/minimal_test/`:
`diag_corrected.jl`, `diag_vd.jl`, and `diag_baseces.jl`. They regenerate every
number above and may be removed once the model is fixed.
