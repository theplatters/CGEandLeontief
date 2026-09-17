# BETA — elastic total labour supply

Implementation: `src/closures/labor/types.jl` (`ElasticLaborClosure`),
`src/closures/labor/labor.jl` (`labor_market_residual`, `solve_beta`,
`solve_verified`), hook consumed by `src/core/equilibrium.jl` (`problem`).
Status lives only in `registry/closures.toml` (`[labor.BETA]`) and renders on
`docs/status.md` (ADR-0003); this page carries no status.

## Formulation (ADR-0002)

`Σ L_i = L̄ · ((w/P)/(w₀/P₀))^η_s`. `η_s` is always the supply elasticity;
never write it as `η` (that is the BF reallocation parameter).

## Design notes (preserved from the frozen `cbase2/src/closures.jl` header)

BETA is the referee-facing principal addition: `L^s = L̄ · [(w/P)/(w₀/P₀)]^{η_s}`,
interpreted as labour–leisure. `η_s = 0` reproduces ALPHA (vertical supply);
`η_s → ∞` approaches the fixed-real-wage regime. It is implemented as a
labour-market-equation hook: the `:beta` closure solves the SAME 2N+1
flexible-wage system as ALPHA with the market-clearing equation `Σ L_i = L̄`
replaced by `Σ L_i = L̄ · (w/w0)^{η_s}`. With the CPI numeraire (P = 1) the
real wage IS w; the anchor is w0 = 1 (baseline numeraire wage).

## What was promoted (Phase 2, ADR-0005)

- `ElasticLaborClosure(η_s; w0 = 1.0)` with unchanged validation (η_s finite
  and non-negative, w0 positive) and the BETA residual hook.
- `MobileLaborCESElasticities` gains `eta_s` (4-arg constructor defaults to
  0.0); `:beta` closure symbol; `labor_closure` returns
  `ElasticLaborClosure(η_s)` for `:beta`.
- `solve_beta` η_s-continuation from the ALPHA equilibrium (ladder, LM/jitter/
  halving fallbacks, fold jump) and `solve_verified`/`_solve_rung_verified`
  (quality gate on the actual residual, never the retcode).
- `mobile_labor_model` gains `financing` and `eta_s` kwargs (`eta_s` forces
  `:beta`); all existing positional/keyword forms keep working.
- `cbase2/src/calibration.jl` was NOT promoted (Phase 3); `solve_beta` here
  runs on the root `Data` with its compatibility defaults.

## Remaining gates and caveats

Source: `registry/closures.toml` (`[labor.BETA]` `open_gates`) and
`cbase2/review.md` §§1, 3.4.

- Contract tests for the promoted code land in the next step (kept in
  `open_gates`; status stays `implemented`).
- The labour–leisure interpretation is not earned by the implementation (no
  income effect; review §1).
- The w0 = 1 anchor vs the baseline wage question stands (review §3.4; the
  v3 equilibrium real wage read ≈ 0.48 under CPI = 1 in cbase2).
- The v3 `solve_beta` η_s-continuation ran > 25 min without completing
  (`cbase2/process_comments.md`, 2026-09-16); each rung can trigger the
  2000-iter LM polish.
- The v2-stage elasticity identification recovers η_s by construction and is
  circular (review §3.4).
