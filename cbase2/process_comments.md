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

# Notebook 03b -- Open-absorption recalibration (v2, implemented 2026-09-16)

- **2026-09-16. v2 implemented and mobile rows verified.** Design per user
  decision (options 2+3, no saving rate, Omega_raw kept): government block
  gG (domestic government consumption), proportional balanced-budget tax
  (F2: tau = (sum gG + sum p g)/(w L), so F2 is a genuine balanced-budget
  multiplier), import margins m_i on household + programme demand (only the
  domestic content circulates; the import content leaks to the external
  account -- R2.4 now answered AT the margin AND at baseline). Household
  block stays the proportional CES (omega = clamped residual shares,
  sum omega = 1 keeps the numeraire consistent). IMPORTANT design lesson: a
  fixed-REAL baseline block breaks Walras add-up at p != 1 (first
  implementation attempt, rejected); the uniform-margin structure is the
  consistent one. Budget identity sum p*c_gross = E now holds EXACTLY
  (= 0.0) at every solved equilibrium. Mobile rows: F1 rel real_gdp ~ 0
  (composition-neutral, eta-invariance diagnostic), F2 -1.70%, F3 +0.057%
  with external deficit F = 0.00283 (import content). tau0 = 0.2068;
  household residual clamped to 0 in 10 government-heavy sectors (mass
  0.0208; parent convention).

- **2026-09-16. MAJOR bug found in the analytic Leontief multiplier (all
  versions).** `leontief_multiplier` solved G y = b instead of
  (I - G) y = b. The v1 "10^14 divergence" was this bug compounded with the
  unit root, not purely the unit root. Exposed by sector 59: an isolated
  node with (M y)_59 = 0 and zero final demand -- the wrong solve assigned
  it y = 6.66. Fixed; residual at the analytic point dropped 6.5 -> 6e-4.
  Finiteness gate (column sums of G < 1) verified: rho(G) = 0.888 (F2),
  0.820 (F3); Keynesian multipliers ~9 / ~5.5; cond(I-G) ~ 15. The v2
  fixed-wage rows are FINITE and well-conditioned -- the unit root is gone.

- **2026-09-16 (cont.). STRUCTURAL FINDING: the v2 import margin without a
  saving/export block has NO interior equilibrium.** Systematic arbitration
  (baseline_check / arbitrate / instrument scripts) established:
  (1) The `p[1]` scale pin OVERWROTE sector N's clearing equation
  (`out[2N]` collision) -- sector 71's market was silently unenforced and
  the numeric solver parked the import leakage there as overproduction
  (y_71 stuck at its init). The numeric "successes" were spurious.
  (2) With all N markets enforced, the Walras identity
  sum_i p_i * residual_i = mbar*E + programme import content (mbar =
  sum m_i cs_i) holds EXACTLY -- the N clearing equations are jointly
  INCONSISTENT for any E > 0: the import content is demand that VANISHES
  (no import supply channel, no export injection). The only solution is
  the corner E = 0, L = sum gG -- exactly what the corrected analytic
  returns (L_ana = 0.2068 = sum gG at g=0, verified residual 0.0).
  (3) Root cause: national accounting. Y = C + I + G + X - M with a
  fully-spending household and balanced government forces p'imp = 0: the
  import leak (mbar*E) has no offsetting injection, so income collapses
  until the leak vanishes. The household's saving rate is NOT optional --
  it is the accounting partner of the export/investment injections:
  S = I + X - M (with T = G). The 2026-09-15 decision to skip the saving
  rate was wrong; the data already implies it (private domestic
  consumption = 1 - gG - X_dom, so s = X_dom/(1 - tau0) ~ 0.55 for
  Germany's open economy).
  DESIGN CONSEQUENCE (v3, needs user sign-off): full open-economy
  Keynesian structure -- exogenous investment + exports (injections,
  import content included), household saving rate s (leakage), government
  (balanced), import margins. Multipliers become SMALL (marginal leakage
  s + mbar + tau ~ 0.75) -- the honest German open-economy answer, and
  the full answer to R2.4. The (I - G) analytic fix stands; the
  leontief_multiplier then needs the X/I/S terms.

- **2026-09-16 (cont.). v3 IMPLEMENTED (open-economy Keynesian structure).**
  User approved: exports exogenous and NO import margin (domestic sales
  abroad); investment I + government G exogenous with margins; household
  consumes (1-s)E, saving sE leaks; calibrated s from the data residual.
  Calibration VERIFIED: s = 0.3979 (in the predicted 0.3-0.45 band),
  tau0 = 0.2141, export share = 0.4219, investment share = 0.1624 --
  Germany's open economy, exactly as predicted. Finiteness gate holds.
  Mobile rows: F1 (composition-only) resid 9e-10, budget identity
  sum p*c - (1-s)E = 0.0 EXACTLY, L = 1.0; F2 (balanced-budget) resid
  5.6e-7, budget exact. Solver hardening: quality gate = ACTUAL residual
  (never the retcode), bounded Levenberg-Marquardt polish (2 attempts,
  maxiters 2000) after a stalled Newton; maxiters 20000 on the main solve.

- **2026-09-16 (cont.). v3 OPEN ITEMS (next session, in this order).**
  (1) Tornqvist real-GDP base is still the v2 form: under v3 the baseline
  household consumption is (1-s)E_h0 * omega = 0.473 * omega, but the
  metric's base vector is the raw omega -- so real_gdp reports 0.199 at
  the v3 baseline. FIX: the Tornqvist base must be
  data.household_baseline (c0_gross) in both _solve_mobile and
  _solve_fixed. Every v3 headline number is contaminated until this is
  done (F2's "-24.6%" is an artifact of the broken base).
  (2) F3 mobile solve stalls at max|resid| = 2e-4 (LM polish insufficient);
  F2 needed no polish (5.6e-7). Investigate: possibly the demand-side
  E(w) coupling stiffened by the (1-s) factor; consider a wage-anchored
  homotopy or a damped Newton.
  (3) The S = I + X - M canary printed diff = -0.16 at the ref: the
  canary code passes w = 1.0 to sectoral_labor_demand (the equilibrium
  wage is sol.wages_raw[1], not 1) and E_ref inherits the error -- the
  canary needs the equilibrium wage; recheck whether the residual identity
  then closes. Also re-derive whether the identity needs p = 1 (it is a
  NOMINAL identity: p-weighted residuals + zero-profit; valid at any p,
  but re-verify with the actual p from the solve).
  (4) The v3 equilibrium real wage: F2 mobile implies w* ~ 0.48 under
  CPI = 1 -- a ~50% real-wage drop at full employment. Either the solve
  found a spurious second equilibrium, or the numeraire/CPI interaction
  with the new demand block needs re-examination. Do NOT trust any v3
  aggregate until (1)-(4) are resolved; the F1/F2 budget identities and
  the calibration are the only fully-trusted v3 numbers so far.
- **2026-09-16 (evening). Mobile reduced-formulation round: ALL SECTIONS
  THROUGH DELTA NOW PASS.** The unified fix: (a) tau_rate for F1/F3
  reverted to the FIXED NOMINAL tax T = sum gG (the "constant rate" idea
  was wrong -- it made the system unsolvable by every algorithm at a
  5.6e-4 residual floor; the assessment's F3 is exactly "baseline tax
  unchanged, programme externally financed", and the saving + import
  leaks close the system without a marginal tax); (b) the analytic F3
  aligned to the same rule (E = L - sum gG) -- DELTA equivalence now
  EXACT for BOTH F2 and F3 (L num/ana = 1.165584/1.165584 and
  1.176709/1.176709, rel y error 0.0); (c) the mobile `problem` keeps the
  FULL 2N+1 formulation (p1..pN, y, w) with ALL N zero-profit + ALL N
  clearing + labour -- no numeraire, no pin (the earlier p1-pin silently
  broke sector 1's zero-profit, letting the wage drift to 0.54); (d)
  injection continuation (exo_scale, bisection start at s = 0) runs in
  5.4 s (was 517 s) with max|p-1| = 0.46 at the endpoint -- the price
  explosion is GONE; (e) F1/F2/F3 mobile rows: resid <= 3.7e-7, budget
  identities exact, F3 external balance F = 0.001472 recorded; (f)
  solver hardening: quality gate on the ACTUAL residual + bounded
  Levenberg-Marquardt polish; q clamped >= 0 post-solve (Tornqvist needs
  nonneg; beta rungs have solver dust); real_gdp returns NaN when
  consumption < 0 (beta rungs with E < 0) instead of throwing.
  REMAINING OPEN (single item): the BETA verification (solve_beta
  eta-continuation under v3) runs > 25 min without completing -- each
  rung now potentially triggers the LM polish (2000 iters, 142-dim FD
  Jacobian) -- needs either a longer timeout run, rung-specific solver
  tuning, or the exo-continuation applied to the beta ladder. Everything
  else in v2_verify passes. The v3 aggregate numbers (w* = 0.54,
  real_gdp = 0.77 at baseline) still need the economic interpretation
  pass -- in p1-units the wage is 0.54 while prices stay near 1; whether
  this is the correct new-normalization reading or a remaining
  inconsistency is the first interpretive question for notebook 04.


- **2026-09-16 (cont.). Continuation diagnostic: the price explosion is a
  CONTINUOUS model pathology, not a solver artifact.** Injection
  continuation (exo_scale 0.47 -> 1, s recalibrated endogenously per step,
  warm-starting each solve; start scale found by bisection on s = 0 since
  s < 0 makes the round-gain matrix non-contractive) tracked the good root
  smoothly -- every step solved to resid <= 5e-7, no stalls. The
  trajectory: max|p-1| grows EXPONENTIALLY (0.97 -> 2.1 -> 5.6 -> 16.6 ->
  52 -> 187 -> 991) while real_gdp (fixed c0_gross base) declines smoothly
  0.976 -> 0.421 and w* falls 0.98 -> 0.55 at constant full employment
  L = 1. So the direct solve's wild prices are the endpoint of a continuous
  branch: as the exogenous injections (I + X = 0.58 of GDP) crowd the
  household's marginal demand share down to (1-s) = 0.60, a price cluster
  in the CES production side (eps = theta = 0.5, complements) explodes and
  the whole economy reprices around it. NEXT SESSION (top priority):
  identify the exploding sector cluster (top-10 prices + their c0_gross /
  exogenous demands / Omega column sums) and choose the structural fix.
  Candidates: (a) isolated/clamped sectors whose price is economically
  irrelevant -- pin p_i = 1 and drop their equations; (b) unit-elastic
  bound on theta (gross-substitute production network, as in BF 2019's
  calibration); (c) re-examine the v3 numeraire interaction. ALSO OPEN as
  before: Tornqvist base fixed this session (base = c0_gross; verify
  real_gdp = 1 at step k=0 of the continuation -- it printed 0.976, close
  but not exactly 1: the clamped-sector baseline is not exactly the
  c0_gross vector at exo_scale > 0 -- recheck); F3 mobile stall at 2e-4;
  canary diff -0.094 with the equilibrium wage (still open -- the
  identity may genuinely fail at the found solution because the dropped
  N-th market absorbs the import content in the MOBILE path too: the
  mobile `problem` still drops the N-th clearing equation, which under v3
  is NOT Walras-redundant -- the reduced all-N formulation exists only in
  problem_fixed).

- **2026-09-16 (cont.). Solver-formulation fix for `:fixed` (done for the
  fixed path; the MOBILE path still drops the N-th market -- apply the same
  reduced formulation there next session).**
  The reduced formulation (p1 = 1 removed from the unknowns; 2N-1
  unknowns; zero-profit 2..N + clearing 1..N) is the correct square system
  once the v3 accounting is in place; sector 1's zero-profit becomes a
  post-solve check. NOT yet applied -- the v2 demand system must be fixed
  first, since the overdetermination it exposes is real, not a solver
  artifact. CD-guard NaN (0/0 for factor_share ~ 1) still open.


- **2026-09-16. `:fixed` solver closure fix (superseded by the structural
  finding below).** With the import margin, the
  N-th market is no longer Walras-redundant, so `problem_fixed` now retains
  ALL N clearing equations and pins the price scale with p[1] = 1 (the CPI
  numeraire is dropped as an equation; checked post-solve). Employment along
  GAMMA/DELTA is determinate. Remaining open (next session): the corrected
  analytic still differs from the numeric solve (F2: 0.209 vs 0.748; F3:
  0.537 vs 0.762) with a 6e-4 residual at the analytic point. Prime
  suspect: the F2 analytic solution has NEGATIVE income (E = -0.01) -- an
  economically invalid branch the numeric's positivity floors never visit;
  the equivalence test may need to be F3-only, or the analytic needs the
  positivity-respecting branch. Also open: CD-guard Diagonal NaN (0/0 for
  factor_share ~ 1 sectors at some solver iterates) -- robustness edge.



- **2026-09-15. Design decisions (user-approved: options 2+3, skip 1, Ω_raw kept).**
  Purpose: give the fixed-wage rows a marginal leakage so GAMMA/DELTA are
  determinate, answer R2.4 *at the margin*, and make F2 a genuine
  balanced-budget multiplier. Specification:
  - **Technology unchanged**: Ω_raw stays the production structure (BF/Domar
    comparability). The economy opens only in *absorption*.
  - **Government baseline**: $gG_i$ = domestic government consumption
    (AC artifact), exogenous fixed real vector, always financed by a
    proportional income tax; baseline rate $\tau_0 = \sum_i gG_i / E_0$.
  - **Exogenous baseline demand**: $X_i$ = domestic equipment investment +
    construction investment + inventories + exports (all exogenous by
    construction, non-income-responsive).
  - **Household**: residual private consumption $c^0 = \lambda - M\lambda -
    gG - X$ (exact baseline clearing, parent's residual convention retained),
    $cs = c^0/E_0$; proportional tax $E = (1-\tau)\,w\sum L$.
  - **Import margins (marginal only)**: household and programme demand carry
    sector import shares $m_i$ (from the §4.1 proportional import split --
    the final-demand domestic fractions by category, sector-composition
    weighted). Marginal domestic content: $dc_i = (1-m_i)\,cs_i\,dE$; the
    programme bundle enters as $(1-m_i)\,g_i$ domestically. The import
    content is recorded on the external account. Baseline levels are the
    calibrated residual (margin applies to changes); documented semantics.
  - **F2 rule**: the proportional tax finances ALL government purchases
    (baseline + programme): $\tau = (\sum p_i gG_i + \sum p_i g_i)/(w\sum L)$.
    F1/F3 keep $\tau_0$. F3's external deficit now records the import
    content of programme + induced consumption, as the assessment promised.
  - **Finiteness gate**: the baseline round-gain matrix
    $M + diag(1-m)\,cs\cdot fs'$ must have all column sums strictly below 1
    (spectral radius < 1) -- the acceptance test for the whole fix.
  - **Consequences**: baseline real GDP is re-normalised (new reference
    solve); notebook 03 smoke numbers must be regenerated under the v2
    calibration before pre-registration; notebooks 01--02 unchanged.

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
