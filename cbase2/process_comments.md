---
title: "Process Comments: Observations by Pipeline Chapter"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-09-15"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [cbase2, process-comments, pitfalls, data-audit]
version: 1
last-updated: "2026-09-15"
---

This file carries the *commentary*: findings, corrections, deviations, and
open questions, organized in chapters corresponding to the pipeline
notebooks. `documentation.md` stays a neutral, timeless description of the
pipeline as it stands; everything here is dated observation and may be
superseded as the pipeline evolves. Entries are additive; nothing is
silently rewritten -- later entries state what they replace.

# Environment (container)

- **2026-09-15.** The container's Julia setup had drifted: the image was
  built while the `release` channel was 1.13, so the baked depot carried an
  environment named `v1.13` and a kernel named `julia-1.13`, while the
  running binary is **1.12.7** (resolved through the sandbox-home juliaup,
  not the baked install). The 1.13-format Manifest also lost
  `git-tree-sha1` fields on some entries, so packages loaded as "not
  installed". Resolution: depot environment renamed to `v1.12` and
  regenerated natively (`Pkg.update()` + full precompile), stale `v1.11`
  environment and ~3.1 GB of unused sandbox-home depot data removed, kernel
  re-registered as `julia-1.12`, build ARG pinned to 1.12.7 in the image
  Dockerfile. Verified: plain `julia` loads the full stack (CSV 0.10.17,
  DataFrames 1.8.2, Plots 1.41.7, NonlinearSolve 4.32.0).

# Notebook 01 -- Data Wrangling

- **2026-09-15. Finding: silent NaN rows in the domestic share matrix.**
  Three sectors record imported intermediates exceeding their total
  intermediate use; after the domestic remainder is clamped, their
  $\Omega^{D}$ rows are 0/0 = NaN. The parent notebook carried these NaNs
  without notice (it asserted only the $Z^{D}$ column sums). Decision: such
  rows are encoded explicitly as **all-zero** ("no domestic content"),
  consistent with the fully-import-dependent interpretation, and the
  row-sum conventions ($\sum = 1$, or $0$ for these rows) are asserted.
- **2026-09-15. Clamps made audible.** The negative-remainder clamp in the
  domestic matrix and the zero-domestic-fraction clamp (final-demand
  bridge, notebook 02) previously passed silently; both now print a
  warning with the affected count/categories.
- **2026-09-15. Standing assumption (carried, documented).** The
  proportional allocation of imports across supplying sectors remains a
  stated assumption: the source table reports imports by *using* sector
  only. Any future Armington-type extension re-opens this.
- **2026-09-15. Porting note.** `DataFrame(A, :auto)` accepts a matrix
  directly; the intermediate `Matrix(A, :auto)` constructor form is not a
  valid call in the current stack (recorded because it cost one
  validation round).

# Notebook 02 -- Accounting Consistency

- **2026-09-15. Correction: residual documentation.** The parent's
  validation section stated the production-vs-expenditure residual as
  "< 1.5%" while the actual value in this data is **5.387%** and the hard
  gate is 10%. The text now states the real figure. The residual itself is
  a raw-table valuation discrepancy (documented, not fixed).
- **2026-09-15. Validation added.** The parent's bridge computed
  `cons_vec`/`cons_share` (the `consumption_share` calibration) without
  any checks; these are now asserted finite, non-negative, and summing to
  one. $\Omega^{raw}$ row sums are re-checked at the bridge.
- **2026-09-15. Open question: shock magnitude.** The defensive load of
  `impulses.csv` (27 rows x 75 columns) sums row 1's numeric columns to
  EUR 80.76 bn, while the paper states ~EUR 40.3 bn at 2019 prices. The
  factor ~2 discrepancy is unresolved: candidates are the row-1 year not
  being 2019, category double-counting in the column sum, or a units
  difference. To be checked before the financing notebooks rely on the
  impulse.
- **2026-09-15. Convention carried.** The expenditure side sums domestic
  final demand **including exports** (the parent pipeline's documented
  convention); it is part of what the 5.387% residual absorbs and is kept
  unchanged for comparability.

# Notebook 03 -- Financing Closures

- **2026-09-15. Design decisions (fixed before implementation).**
  Common programme frame: real bundle $g_i = m\,G_0\,\psi_i$ in numeraire
  units with incidence $\psi_i$ = year-1 (2024) sectoral shares of
  `data_raw/impulses.csv` and $G_0 = 40{,}300$ EUR m (the paper's stated
  programme); 1 model income unit = GDP at basic prices. F1 composes the
  household budget via preference shifters $d_i = 1 + G_0\psi^{F1}_i/(cs_i
  E_0)$ (no additive demand); F2 adds the bundle and levies an endogenous
  lump-sum tax $T(p) = \sum p_i g_i$; F3 adds the bundle untaxed and
  records the external balance $F = \sum p_i g_i$. F1's incidence is
  restricted to consumable categories ($cs_i > 0$; 10 zero-$cs$ sectors,
  zero impulse mass in this data) and renormalized.
- **2026-09-15. Recorded kernel DIFFs (hybrid drift rule).** The
  financing hook required two surgical edits to the core copies, flagged
  `DIFF` by `scripts/diff_kernel.jl`: `src/core/interface.jl`
  (AbstractFinancing supertype, `Model` gains `financing`, 3-arg
  constructor defaults to `NoFinancing`) and `src/core/mobile_labor.jl`
  (both residual systems route demand through the financing hooks; the
  $\eta=1$ fixed-wage scale-indeterminacy guard now keys on
  `has_additive_anchor`; post-solve consumption mirrors the financed
  expenditure). Deliberate design: ONE equilibrium system, so Stage 1.6
  validation exercises the system that produces the matrix. Backport to
  the parent src/ is pending and should follow once the Stage 2 numbers
  are stable.
- **2026-09-15. Residual-function lesson.** A positivity guard on the
  household expenditure base $E$ thrown *inside* the residual function
  kills the solve: the nonlinear solver legitimately explores
  negative-income trial points, and the legacy code tolerated that
  silently. Guards belong post-solve (headline assertions), not in
  residuals.
- **2026-09-15. Headline metric finding.** The B&F Tornqvist real GDP
  metric is a **household consumption** index. Under tax financing (F2)
  it contracts by roughly the programme share even when total final
  demand is flat, and under F3 it misses the externally financed bundle.
  The notebooks therefore record a supplementary **total-final-demand
  index** ($c + g$) next to the consumption metric. At $m = 1$
  (programme = 1.33% of GDP): F1 leaves both metrics flat (the
  $\eta$-invariance signature); F2 mobile shows full-employment
  crowding-out (consumption $-1.35\%$, total demand flat); F3 mobile
  keeps consumption $\approx$ flat with the external deficit absorbing
  the programme ($-0.06\%$ consumption, $+1.26\%$ total demand,
  $F = 1.33\%$); F2 under the sticky wage is nearly
  employment-neutral ($+0.008\%$) -- the retired +19.3 pp extensive-margin
  result was **unfinanced** manna, and the financed version must be
  reported against the pre-registration as a deviation finding (Stage 2).
- **2026-09-15. Open question narrowed: 80.76 vs 40.3 bn.** The 2024
  impulse row sums to EUR 80.76 bn (goods only; the "Wages and salaries"
  column is zero throughout). Two candidate explanations for the paper's
  ~EUR 40.3 bn: the **2046 impulse row** (40.372 bn, the decaying
  programme's mid-range year) or an **exact halving** (80.76/2 = 40.38,
  e.g. one of two funding sources or a net-of-something convention). The
  manuscript text must settle this; $G_0 = 40{,}300$ is used meanwhile
  per decision.
- **2026-09-15. Container env.** `ProgressMeter`/`ThreadsX` (kernel deps
  of util/variance_decomposition) added to the depot environment via the
  documented runtime `Pkg.add` route.

# Notebook 04 -- Labour Closures (BETA complete; DELTA blocked on a structural finding)

- **2026-09-15. BETA implemented and verified.** Elastic total labour supply
  enters as a labour-market-equation hook (`:beta` closure, `eta_s` on the
  elasticities struct, `solve_beta` η_s-continuation with fallback
  strategies). Verified: η_s = 0 reproduces ALPHA exactly (same financing);
  the implied elasticity $\hat\eta_s = \log L / \log w$ recovers
  0.5 / 1.0 / 2.0 (2.00001 at the extreme) from single equilibria -- no
  demand variation needed, because the numeraire pins the anchors at
  $w_0 = 1$, $\bar L = 1$. User question that triggered the simplification:
  the m-variation in the elasticity test was a test-design artifact, not a
  pipeline requirement; m-variation returns only in Stage 2 robustness
  sweeps.
- **2026-09-15. Cobb-Douglas limit guard implemented.** Analytic branch
  `cost = w^{fs}·ip^{1−fs}/A` at $|\;1-\epsilon\;| < 10^{-9}$ inside the
  shared `_ces_unit_cost` helper (both residual systems); smoke test pending
  re-run after the DELTA blocker below is resolved.
- **2026-09-15. STRUCTURAL FINDING (blocks DELTA, possibly GAMMA).** The
  analytic Leontief counterpart of the DELTA corner exposed a rank problem:
  under `:fixed` with demand-only shocks, prices are pinned by zero-profit
  ($p = 1$ exactly, demand-independent), intermediates are linear in $y$,
  and household consumption is proportional to wage income -- so the
  clearing system takes the form $y = (M + cs\cdot\kappa')y + g$ whose
  column sums are exactly 1 at the baseline calibration ($\kappa = fs$).
  The system is then **rank-deficient (unit root): the uniform-expansion
  direction is unresolved, and the employment level along the GAMMA/DELTA
  rows may be indeterminate -- chosen by the solver's warm-start path, not
  by the equilibrium**. Consistent symptom: the DELTA solve reports
  machine-zero residuals with employment FALLING (0.9957) under additive
  external demand, while the over-determined analytic check (all N clearing
  equations) diverges ($L = 10^{14}$). Economic reading: with zero savings,
  zero taxes at the margin, and fixed prices, the marginal propensity to
  spend is 1 -- the Keynesian multiplier is at the knife edge. Possible
  resolutions to investigate: (i) the model needs a leakage (saving, or a
  proportional -- not lump-sum -- tax rule) for the fixed-wage rows to be
  determinate; (ii) the CPI numeraire or the reallocation wedge must enter
  the fixed-wage system (the wedge is currently not applied under
  `:fixed`); (iii) an index error in the analytic check. This MUST be
  resolved before Stage 2: if (i)/(ii), the pre-registered "+19.3 pp
  extensive margin" (itself unfinanced and retired) and any GAMMA/DELTA
  employment levels are warm-start artifacts. The DELTA/CD smoke cells and
  notebook 04 are therefore NOT built yet.
- **2026-09-15. Solver lessons (BETA).** (i) Never gate solutions on the
  NonlinearSolve retcode: `Stalled` returns can already satisfy the 1e-6
  residual gate. (ii) Do not tighten tolerances beyond the proven
  `reltol = 1e-6` configuration -- tighter settings stall. (iii) Large
  parameter jumps stall where small continuation steps solve instantly;
  homotopize in $\eta_s$ (`solve_beta`) and in $m$ where needed.

# Notebooks 05--08 (preregistration, validation, matrix, sensitivity)

- No observations yet; chapters fill as the notebooks are built. Expected
  first entry for 03: the F1/F2/F3 equation design decisions and their
  rationale, recorded before implementation.
