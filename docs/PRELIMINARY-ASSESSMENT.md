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

> **Historical note (2026-09-03).** The `tests/minimal_test/` diagnostic scripts
> referenced throughout this document (`diag_review.jl`, `diag_milestoneA.jl`,
> `diag_B.jl`, `diag_vd.jl`, `diag_corrected.jl`, `diag_baseces.jl`,
> `run_test.jl`, `test_full_pilot.jl`) have since been **deleted** and replaced by
> the organised test suite (`tests/test_mobile_labor.jl`,
> `tests/test_fixed_closure.jl`, `tests/test_variance_decomposition.jl`,
> `tests/test_eta_sweep.jl`). The empirical claims here are retained as a historical
> assessment of the pre-fix code; they are not evidentially admissible for the
> corrected model. Regenerate any number via
> `julia --project=. -e 'using Pkg; Pkg.test()'`.

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

Four reproducible diagnostic scripts remain in `tests/minimal_test/`:
`diag_corrected.jl`, `diag_vd.jl`, `diag_baseces.jl`, and `diag_review.jl`. They
regenerate every number above and may be removed once the model is fixed.

# Addendum: Independent Verification of the Second Review (2026-08-18)

A second, independent review of the same codebase was received and saved at
`roadmaps/vertdict.md` (the filename is a typo for "verdict"). It diagnoses the
mobile-labor model as not an equilibrium and proposes a redesign. This addendum
records that I verified its code-line claims against the actual source and
reproduced its central diagnosis with an independent 71-sector computation.

## Verification of the review's code-line claims

Every specific code-line claim in the review was confirmed by reading the source:

| Review claim | Source line | Confirmed |
|---|---|---|
| Demand shock applied as unnormalized expenditure shifter | `mobile_labor.jl:148/150` | Yes |
| Wrong equation omitted: sector-1 zero-profit, not a goods-market equation | `mobile_labor.jl:158/159` | Yes |
| Solver return code not checked | `mobile_labor.jl:197/198` | Yes |
| Nominal `w` used in labor supply (should be `w/P`) | `mobile_labor.jl:165` | Yes (see Corrections) |
| High-$\eta$ limit misinterpreted: $w < 1 \Rightarrow L \to 0$ | mechanism | Yes |
| Fixed-base sum instead of shared Tornqvist index | `mobile_labor.jl:218` | Yes |
| Incomplete design silently treated as balanced | `variance_decomposition.jl:261` | Yes |
| "% Share" renormalized over main effects | `variance_decomposition.jl:308` | Yes |

## Reproduction of the central diagnosis

I rebuilt the equilibrium residuals directly from a solved 71-sector
mobile-labor model at $\eta = 0$ (standard construction shock). The only
violated equilibrium condition is the **omitted sector-1 zero-profit equation**:

| Check | Result |
|-------|--------|
| Household overspend (71-sector, $\eta=0$) | ratio 1.0449 (4.5% above income) |
| Sector-1 zero-profit residual | 10.924 (omitted equation; should be ~0) |
| Max zero-profit residual, sectors 2..N | 6.6e-5 (enforced) |
| Labor-market residual | 5.2e-5 (cleared) |
| Numeraire (CPI-1) residual | 0.0 (satisfied) |
| Max goods-market residual | 0.000236 (all N cleared) |
| Total residual norm | 10.924 (entirely the omitted sector-1 equation) |

The mechanism is therefore confirmed: the unnormalized `demand_shock` makes
household expenditure exceed income by 4.5%, so Walras' law cannot recover the
omitted condition, and the missing zero-profit equation permits a sector-1
profit residual of 10.924. This matches the review's two-sector finding
(expenditure 1.25485, sector-1 residual ~0.39); the magnitudes differ only
because the 71-sector `demand_shock` is a *level* vector (1.0 baseline plus a
bump), not a deviation vector.

The diagnostic is at `tests/minimal_test/diag_review.jl`.

## Relationship to the assessment above

The review **supports and extends** the assessment; there are no contradictions.
Both conclude that every current mobile-labor result -- including the "GO"
decision and the "88.4%" -- is invalid until regenerated. The review is sharper
on mechanism and redesign; this assessment adds independent empirical
reproduction and the workplan-versus-ROADMAP contradiction. Specifically, this
assessment adds:

- The empirical $\eta$-sweep breakdown ($\eta=10 \to$ GDP 0.123, i.e. -86%) and
  the initialization-dependent equilibria (GDP 0.12 or 1.13 for the same
  $\eta$).
- The Cobb-Douglas `DomainError` at $\varepsilon = 0.99$ behind the silent
  `NaN` dropping.
- The fact that the current code returns $\eta = 100.0\%$ in the "% Share"
  column, not 88.4% -- so the workplan's 88.4% is from an earlier or different
  code state.
- The observation that `run_test.jl` prints "Solved successfully" for a
  non-equilibrium, the practical face of the unchecked-return-code point.

## Corrections

Two clarifications to the discussion above and in chat:

1. **Overspend magnitude.** An earlier chat message implied an ~81% household
   overspend. That was wrong: `demand_shock` is a level vector (1.0 + bump),
   not a deviation, so the verified 71-sector overspend is **4.5%** (the
   review's two-sector figure is 25.5%). The mechanism -- an unnormalized shock
   breaks the income-expenditure identity -- is correct; only the magnitude was
   misstated.
2. **Nominal versus real wage (review point 4).** Because CPI = 1 is the
   numeraire, `w` is already the real wage numerically. The decisive fix is
   anchoring at the baseline wage,
   $L = \bar L [(w/P)/(w_0/P_0)]^\eta$, which alone repairs the
   $\eta \to \infty$ collapse. Do not over-rotate on the nominal/real
   distinction.

## Redesign and ROADMAP recommendations from the review

Condensed from `roadmaps/vertdict.md`:

- **Separate household demand from policy investment.** Budget-neutral
  preference reallocation uses renormalized CES weights
  $\tilde\beta_i = \beta_i d_i / \sum_j \beta_j d_j$ with
  $Z(p) = \sum_j \tilde\beta_j p_j^{1-\sigma}$ and
  $c_i^h = E_h \tilde\beta_i p_i^{-\sigma} / Z(p)$, guaranteeing
  $\sum_i p_i c_i^h = E_h$. Green investment uses a separate vector $g_i$ with
  $\sum_i p_i g_i = T + B + F$ (tax, borrowing, or foreign financing).
- **Correct square system.** Retain $N$ zero-profit, $N-1$ goods-market, one
  labor/closure, and one numeraire equation ($2N+1$ total). After solving,
  compute all $N$ goods-market residuals, including the omitted one; the
  solution is valid only if every residual, budget, and solver status is
  satisfactory and all prices, outputs, wages, and consumption are strictly
  positive.
- **Labor closures as separate types:** FixedLabor ($L = \bar L$),
  ElasticLaborSupply ($L = \bar L [(w/P)/(w_0/P_0)]^\eta$), FixedRealWage
  ($w/P = \bar\omega$), and UnemploymentComplementarity
  ($0 \le \bar L - L \;\perp\; w/P - \bar\omega \ge 0$).
- **Sensitivity redesign.** Replace the current estimator with a proper
  first-order index $S_f = \mathrm{Var}(\mathbb{E}[Y \mid f]) / \mathrm{Var}(Y)$
  under an explicit, documented probability measure; report total-effect and
  selected interaction indices, solver status, and residuals per evaluation;
  never silently drop failed cells; remove the second "% Share" normalization.
- **2019 versus 2022.** Keep a corrected 2019-style competitive CES core;
  import only a stripped-down 2022 complementarity closure for the unemployment
  extension. Switch entirely to the 2022 base only if the paper's question
  changes from "how closure choice moves results between CGE and IO endpoints"
  to "how sectoral wage rigidity and monetary conditions govern green-stimulus
  effects."
- **ROADMAP changes.** Mark the two recovery/audit tasks complete but record
  the audit outcome as rejection of their current results; state explicitly that
  no mobile-labor figures, GO decision, or sensitivity shares are evidentially
  admissible; distinguish preference reallocation, tax-financed investment,
  expenditure-switching investment, and debt/foreign-financed expenditure in
  Phase 2; add shared final-demand and institutional-budget components before
  labor closures in Phase 3; replace the finite $\eta=\infty$ approximation with
  the explicit fixed-real-wage closure; make omitted-equation invariance and
  household-expenditure exhaustion separate Phase 4 tests; require a documented
  product probability measure and complete designs in Phase 5.

## Disposition

`WORKPLAN2.md:16` ("complete and verified", GO decision, price-invariance
finding, 88.4% result) should be recorded as **superseded**: its claims are not
supportable by the current code. The critical path from the review -- accounting
reconciliation, shock definition and financing, shared equilibrium core, labor
closures, full residual validation, sensitivity analysis -- should replace the
optimistic "GO" framing in the workplans.
