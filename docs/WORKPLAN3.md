---
title: "Workplan: Repairing the BeyondHulten Mobile-Labor Core"
date: "2026-08-18"
---

> **Historical note (2026-09-03).** The `tests/minimal_test/` diagnostic scripts
> referenced below (`diag_review.jl`, `diag_milestoneA.jl`, `diag_B.jl`,
> `diag_vd.jl`) have since been **deleted** and replaced by the organised test
> suite (`tests/test_mobile_labor.jl`, `tests/test_fixed_closure.jl`,
> `tests/test_variance_decomposition.jl`, `tests/test_eta_sweep.jl`). The
> `output/*.csv` files mentioned are **locally generated, gitignored artifacts**,
> not committed evidence. The numerical results recorded here (Milestones A–F)
> predate the 2026-09-03 fixes (log-space mobile allocation, `Ω_raw`
> calibration/equilibrium alignment, Sobol grid/output validation, `:fixed`
> scale-indeterminacy guard). The **2026-09-03 re-run is now complete** under the
> fixed code/data; its current numbers supersede the historical diagnostics below
> and are summarised in the addendum at the end of this document. No fabricated
> numbers appear here.

# Scope and Verdict

This workplan repairs the core issues documented in `PRELIMINARY-ASSESSMENT.md`
and corroborated by `roadmaps/vertdict.md`. It targets the revision-support code
(`src/`, `tests/test_mobile_labor.jl` and the other `tests/test_*.jl` files) only.
The two `bf_replication/` directories are excluded by agreement.

Verdict adopted from the second review: build a corrected, budget-complete
2019-style competitive CES core; make labor and financing closures modular; use
a stripped-down 2022-style complementarity closure only for the unemployment
extension. Do not replace the model wholesale with Baqaee-Farhi (2022).

# Core Issues to Fix

1. Mobile-labor "bridge" is inverted: $\eta \to \infty$ yields a GDP collapse
   instead of the expected $+1.5\%$ input-output multiplier. Root cause is the
   non-equilibrium model, not the economics.
2. Real GDP in `MobileLaborCES` uses a Laspeyres index
   (`sum(consumption)/sum(consumption_share)`) instead of the mandated
   Tornqvist/Divisia index.
3. `solve()` passes `complex`-typed CES powers to `nlsolve` with
   `autodiff=:central`, producing a `DomainError` at the Cobb-Douglas limit
   ($\varepsilon \to 1$).
4. `partial_r2` is a one-way marginal $R^2$ with no interaction or conditional
   handling; `summary_table` renormalizes over main effects (double counts). The
   cited "$88.4\%$" is a renormalized main-effect share, not a variance share.
   The current code's $100\%$ ($\eta$) result is an artifact of the broken GDP
   collapse.
5. Sector-1 zero-profit equation is omitted (`out[1:N-1]` never sets `out[1]`);
   `solve()` ignores `retcode`; household demand is unnormalized, so the model
   is not an equilibrium (Walras' law cannot recover the omitted condition).

# Milestones

## Milestone A - Complete and Anchored Equilibrium (DONE 2026-08-18)

Changes in `src/mobile_labor.jl`:

- Restore the sector-1 zero-profit equation: `out[1] = p[1] - cost[1]`
  (currently only `out[1:N-1] .= p[2:N] .- cost[2:N]`).
- Specify labor supply with a baseline anchor:
  $L = \bar L\,(w/w_0)^\eta$ with $w_0 = 1.0$ (numeraire baseline), so the
  competitive limit ($\eta \to \infty$) keeps $w$ at the numeraire and lets $L$
  absorb shocks, while $\eta = 0$ pins $L = \bar L$ and forces $w$ to clear.
- Normalize household final demand to income:
  `final_demand = (consumption_share ./ sum(consumption_share)) .* E_h ./ p`,
  with `E_h = sum(consumption_share .* p .* x)` evaluated at baseline.
- (Deferred to Milestone E.) The scale/financing channel: the code currently has
  NO investment or financing term, so the "debt-financed investment" assumption is
  not yet implemented. This is where the bridge must be built (see Milestone E).
- In `solve()`, assert or return `retcode`; fail loudly on divergence instead
  of silently returning a degenerate solution.

Validation gate A (historical — `diag_review.jl` was deleted; use
`tests/test_mobile_labor.jl`):

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expect: every residual (including the sector-1 zero-profit residual) at $\sim
10^{-5}$; total residual norm $\approx 0$. Baseline solve: real GDP $\approx
1$, `sum(labor_share) == 1`.

Status (2026-08-18): VERIFIED (historical — `diag_milestoneA.jl` / `diag_B.jl`
were deleted; the same checks now live in `tests/test_mobile_labor.jl`). Residuals
$\sim 10^{-15}$ (machine precision) for $\eta = 0 \dots 50$; sector-1 zero-profit
satisfied; `solve()` fails loudly on non-convergence.
NOTE: with the current composition-only demand shock the model is $\eta$-invariant
(real GDP $\Delta = +0.0667\%$, wage $= 1.0$ for all $\eta$) -- the bridge is
absent because no scale/financing channel exists yet. That is Milestone E.

## Milestone B - Correct Real-GDP Index (DONE 2026-08-18)

- Replace the Laspeyres real-GDP expression in `MobileLaborCES` with the already
  implemented `tornqvist_quantity_index` from `src/solution.jl` (consumption-based,
  matching `ces.jl`).
- Validation gate B: baseline real GDP $= 1.0$; a $+81\%$ construction demand
  shock gives a small, correct-sign real-GDP response ($\approx +0.07\%$ under the
  current composition-only shock -- see the historical `diag_B.jl`, now
  `tests/test_mobile_labor.jl`). The base-CES $+3.26\%$
  replication (gate as originally written) requires the expansive shock from
  Milestone E, not a composition shock.

Status (2026-08-18): VERIFIED. `real_gdp(sol)` is now the Tornqvist index;
baseline $= 1.000000$.

## Milestone C - Solver Stability

- Remove `complex`-typed powers from the CES aggregator; use a sign-safe real
  power (for example `exp(y * log(max(x, eps())))` with proper bracketing), or
  reformulate the aggregator to avoid complex intermediates. Set
  `autodiff = :finite` in `nlsolve`.
- Validation gate C (historical — `diag_vd.jl` was deleted; use
  `tests/test_variance_decomposition.jl` and `tests/test_eta_sweep.jl`): the full
  $\eta \times \varepsilon \times \theta \times \sigma$ sweep completes without a
  `DomainError`; the $\eta$-sweep is monotonic and bounded.

## Milestone D - Honest Variance Decomposition (VERIFIED 2026-08-19)

Sobol first-order (S_f) and total-order (ST_f) indices on a 4×2×2×2 = 32-point
full factorial grid over (η, ε, θ, σ). Method: ANOVA sum-of-squares on a
balanced full factorial design. S_f = SS_f / SS_total (first-order). ST_f =
1 − SS_{-f} / SS_total (total-order, captures all interactions involving f).
Absolute shares reported (no renormalization over main effects). Results saved
to CSV at `output/variance_decomposition_sobol.csv` (locally generated, gitignored
artifact — not committed; regenerate via `rerun_results.jl` / the test suite).

RESULT (sector-35 construction supply +30%, ε ∈ {0.5, 2.0}, θ ∈ {0.5, 0.99},
σ ∈ {0.5, 0.99}, η ∈ {0, 0.5, 1, 2}):

  Factor   S_f (first-order)  ST_f (total-order)  ST_f−S_f (interactions)
  η        0.0001             0.0050              0.0049
  ε        0.0990             0.1091              0.0101
  θ        0.5186             0.5239              0.0053
  σ        0.3720             0.3746              0.0026
  Sum S_f  0.9897
  1−Sum S_f = 0.0103 (interaction share)

FINDINGS:
1. θ (intermediate-good substitution) dominates at 51.9% of variance.
2. σ (consumption substitution) is second at 37.2%.
3. ε (labor/intermediate substitution) at 9.9%.
4. η (labor mobility) at 0.01% — the reallocation bridge is so small that
   η barely contributes to GDP variance in the competitive one-factor model.
5. Interactions are modest (1% of variance — 1−Sum S_f = 0.0103).
   η's total-order index (0.005) is only marginally above its first-order,
   meaning η's tiny effect is mostly a direct main effect.

Comparison to earlier artifacts:
- The original "η=100%" result was from the broken model (GDP collapse at high η).
- The "η=88.4%" was a renormalized main-effect share that masked interactions.
- The honest Sobol S_f gives η ≈ 0.01%, ruling out η as a driver of GDP variance
  in this competitive one-factor specification.

Validation gate D: PASS. η share is sensible (essentially zero), residual is
reported (0.0103 = interaction share), and the result is not the 100% artifact.
All 32 grid points solved successfully (0 failures, maxiters=1000, reltol=1e-6).

## Milestone F — Binary Wage Regime (VERIFIED 2026-08-19)

The `:fixed` (sticky-wage) closure pins the real wage at 1.0 and lets
employment absorb the shock — a different economic mechanism (extensive
margin) from the η mobility channel (intensive margin). Previously the
`:fixed` closure failed for moderate-to-large shocks because the wage was
treated as a 2N+1 unknown with a linear constraint, creating a near-singular
Jacobian. Fixed with a specialized 2N solver (`problem_fixed` + `_solve_fixed`)
that hard-codes w = 1.0 and solves for prices and quantities only.

RESULT (autonomous demand shock to construction, ε=0.5, η=0.5):

  mult    η=0 Δ%    :fixed Δ%   fix_bridge(pp)  employment  resolution
  0.1     +0.002    +0.251      +0.249          1.0025     1.4e-11
  0.2     +0.009    +0.635      +0.626          1.0064     9.9e-07
  0.5     +0.047    +1.203      +1.156          1.0120     9.4e-09
  1.0     +0.159    +2.287      +2.129          1.0229     2.0e-07
  2.0     +0.467    +4.228      +3.761          1.0423     1.0e-07
  5.0     +1.525   +11.362      +9.837          1.1136     6.6e-11
 10.0     +3.146   +21.661     +18.514          1.2166     7.0e-10

FINDINGS:
1. The binary wage regime (sticky vs flexible) produces an aggregate GDP
   effect 74–112× LARGER than the continuous η gradient for small shocks.
2. At large shocks (mult=10), :fixed gives +21.7% GDP vs :mobile's +3.1%
   — a 7× amplification through employment expansion (emp→1.22).
3. The η gradient (η=0 → η=1 within the :mobile closure) never exceeds
   0.01 pp — negligible by comparison.
4. Labour-market closure matters for aggregate GDP, but through the WRONG
   margin relative to the original hypothesis: wage regime, not mobility.

CONCLUSION: The paper's central result must shift from "η determines the
aggregate multiplier" to "the binary wage regime (sticky vs flexible) is the
relevant labour-market margin; η is a second-order correction." Results
saved to `output/wage_regime_comparison.csv` (locally generated, gitignored
artifact — not committed; regenerate via `rerun_results.jl` / `tests/test_fixed_closure.jl`).

Validation gate F: :fixed solves for all tested multipliers (0.1—10.0) with
max residuals ~1e-7 to 1e-10. The employment and GDP responses are monotonic
in shock size and economically interpretable.

# Verification Protocol

Each milestone is validated by the organized test suite (`tests/test_*.jl`),
superseding the deleted `tests/minimal_test/` diagnostics, and the gate
listed above. No "GO" is asserted until all gates pass and the equilibrium
residuals are $\sim 10^{-5}$. The false certifications in `WORKPLAN.md` and
`WORKPLAN2.md` (line 16 of each) are superseded by this plan.

## 2026-09-03 Current-Rerun Addendum (supersedes the historical Milestone D/F numbers above)

The empirical re-run executed on **2026-09-03** under the fixed code/data
(`rerun_results.jl` + `tests/`) completed successfully. Its results **supersede**
the 2026-08-19 Milestone D (Sobol) and Milestone F (wage regime) diagnostics
recorded above, which were produced before the 2026-09-03 fixes. The reproducible
local artifacts (`rerun_results.log`, `output/variance_decomposition_sobol.csv`)
are **gitignored and not committed**; regenerate with
`julia --project=. rerun_results.jl`.

**Canonical accounting result:** GDP P = 3,027,818; I = 3,027,818; E = 2,864,724;
residual 5.387% (documented raw-table valuation gap); `sum(λ) = 2.1099`;
`sum(labor_share) = 1`; reference `ref_gdp = 0.9998546537504605`.

**§3 wage-regime table (autonomous demand to construction; `:fixed` residuals ≤ 1.9e-7):**

| mult | η=0 Δ% | η=1 Δ% | :fixed Δ% | η_bridge (pp) | fix_bridge (pp) | fixed_emp |
|:----:|:------:|:------:|:---------:|:-------------:|:---------------:|:---------:|
| 0.1 | -0.045697 | +0.014512 | +0.271460 | +0.060209 | +0.317157 | 1.0026 |
| 0.2 | -0.051978 | +0.014496 | +0.518147 | +0.066473 | +0.570124 | 1.0050 |
| 0.5 | -0.090681 | NaN (stalls) | +1.227407 | NaN | +1.318088 | 1.0121 |
| 1.0 | -0.203714 | NaN (stalls) | +2.318170 | NaN | +2.521884 | 1.0230 |
| 2.0 | -0.542053 | NaN (stalls) | +4.254413 | NaN | +4.796466 | 1.0424 |
| 5.0 | -2.009894 | NaN (stalls) | +9.957729 | NaN | +11.967623 | 1.0994 |
| 10.0 | -5.406077 | NaN (stalls) | +19.344778 | NaN | +24.750854 | 1.1933 |

The flexible η=0 closure turns **increasingly negative** with shock size; the sticky
`:fixed` closure gives a large **positive** extensive-margin response (employment up
to 1.19×). Because η=0 is negative, the old "5–7×" / "74–112×" amplification *ratios*
are **not used**; the comparison is the **percentage-point gap** (fix_bridge above).
η=1 (mobile) solves only at mult 0.1/0.2 and **stalls** from 0.5 upward (no
autonomous/investment anchor — the `:fixed` scale-indeterminacy guard).

**§2 Sobol (32/32 grid points, 0 failures):** θ dominant (S_f = 0.3951, ST_f = 0.3972);
σ second (S_f = 0.2734, ST_f = 0.2744); ε (S_f = 0.1650, ST_f = 0.1719) and η
(S_f = 0.1571, ST_f = 0.1662) each ~16%. ΣS_f = 0.9907; interaction share = 0.0093.
**η is material (~15.7% first-order), not negligible and not dominant.** This corrects
both the original "η = 88.4%/100%" artifact and the 2026-09-02 "~2% negligible" reading.

# Out of Scope

`bf_replication/` and `bf_replication2/` are excluded. The narrative
reframing in `ROADMAP.md` (two-referee, 2019-vs-2022) is retained; only the
numerical core is repaired here.
