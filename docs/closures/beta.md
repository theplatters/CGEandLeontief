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

## The sectoral form (ADR-0022, 2026-09-19)

`SectoralElasticLaborClosure(η_s_vec)` — the N-market generalisation of the
scalar closure. The single aggregate equation is replaced by N sectoral ones,

    L^cm_i(p, y, w_i) = Lbar_i · (w_i / P)^{η_s,i},   i = 1..N,   P = CPI,

solved as a `3N+1` system `[p; y; w(1:N); F]`. The closure is selected by
`MobileLaborCESElasticities.eta_s_vec` (`nothing` = the scalar 2N+2 form), so
the registry id stays `BETA`: the sectoral form is BETA's N-market
generalisation, not a new closure.

**Endpoints.** `η_s,i = 0` for every sector is *exactly* the `η = 0` endpoint of
ADR-0020 option C — the allocation frozen by a vertical supply curve instead of
by a datum — proved analytically in ADR-0022 and asserted in
`tests/test_sectoral_labour.jl` (`0.00e+00` when warm-started from the endpoint
solution, `< 1e-10` cold). The `η_s,i → ∞` limit is the GAMMA corner (the real
wage pinned, employment absorbing); it is a limit, not a cell — the direct
formulation stiffens and eventually fails to solve.

**Why it is the demand-sensitive form.** With one wage, zero profit (N
equations) plus the CPI numeraire pin `p` and `w` under demand-only shocks, so
the price block is demand-free (this is the `BETA ≡ ALPHA` degeneracy). With N
sectoral labour markets the N wages are tied to N employment levels, so demand
enters the price block through the wage structure: `max|p−1|` differs across
F1/F2/F3, which no single-wage closure achieves. The measured ladder
(`matrix_5x3-v9-BETA-F2-etas*`): 0.152897 (0.25), 0.105666 (0.5), 0.065422 (1),
0.037171 (2). Rigid-group variants: programme sectors rigid 0.272033, largest
half rigid 0.109583 (F1).

`η_s,i = 0` for every i is exactly the ADR-0020 option C endpoint (the executed
`matrix_5x3_v6` BF row) — the promotion's nesting canary; `η_s,i → ∞` is
GAMMA-like in prices, approached only asymptotically (the supply slope
diverges).

Unlike the scalar form, the sectoral form moves prices with the demand
composition: on full-71 a uniform `η_s = 0.5` gives `max|p-1|` ≈ 0.10 across the
financing columns, against `2.7e-15` for the scalar closure at the same
elasticity. The scalar row stays degenerate under demand-only shocks *by
design*; the sectoral vector is the labour door of `docs/ETAs.md` candidate 1 /
`paper/equivalence.tex` Section 6.

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
