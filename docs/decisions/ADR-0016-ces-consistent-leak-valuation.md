# ADR-0016 — CES-consistent valuation of the intermediate-bill leaks

- **Status:** accepted (user instruction, 2026-09-18)
- **Date:** 2026-09-18
- **Related:** ADR-0010 (the canary is the asserted acceptance identity),
  ADR-0012 (A-bill), ADR-0013 (row 75 as the third leak),
  `src/core/equilibrium.jl` (`_mobile_market_demand`, `problem_fixed`,
  `external_balance_canary`), `tests/test_promoted_closures.jl`,
  `tests/test_calibration.jl`; found in the 2026-09-18 review of commit range
  `5ea54465..901ec07`

## Context

ADR-0012/0013 split the per-user purchaser-price intermediate bill as
`(1−fs_u)·λ_u ≡ A_bill_u + M_int_u + T_int_u`, charge only `A_bill` to
domestic demand, and book `M_int` (row 74) and `T_int` (row 75) as
external-account leaks.

The domestic demand system already values the bill with the CES expenditure
factor `k_u = p_u^ϵ · a_u^(ϵ−1) · P_u^(1−ϵ)` (collapsed from the conditional
CES demand and the CES price index; `a` = supply shock,
`P` = `_intermediate_price(Ω_raw, p, θ)`). `external_balance_canary`, however,
valued the two leaks at `p_u · (bill_u/λ_u) · y_u`, i.e. with the factor
effectively frozen at `k_u ≡ 1`. Because both leaks are constant shares of
the SAME per-user CES intermediate bundle, the identity from ADR-0010
(review finding 2.1) could then only close in the baseline/demand-only gauge
(`p ≡ 1`, `a ≡ 1`).

The measured consequence, on the full-71 A-bill calibration with a +20 %
sector-1 supply shock: equilibrium residual 5.26e-13, omitted-market residual
7.966554152441511e-4 vs canary 9.380591541865158e-4 — a gap of
−1.4140373894236467e-4, above the 1e-6 `assert_external_canary` gate
(`experiments/run.jl`), so a correctly solved equilibrium was rejected.
Supply-side scenarios — the regime the BETA `η_s` identification needs —
would misreport the external position.

A one-sector exact root makes the mechanism transparent: with fs = 0.5,
λ = 2, A_bill = 0.5, M_int = T_int = 0.25, supply shock a = 4, ϵ = θ = 0.5,
p = 1, y = 12, w = 9, all imposed residuals are exactly zero, the
omitted-market residual is 1.5, and the pre-fix canary reads 3.0 because
`k = a^(ϵ−1) = 0.5`.

## Decision

`external_balance_canary` computes
`k_bill = p .^ ϵ .* shocks.supply_shock .^ (ϵ - 1) .* _intermediate_price(data.Ω_raw, p, θ) .^ (1 - ϵ)`
and values both leaks as `dot(k_bill .* (data.M_int ./ data.λ), y)` and
`dot(k_bill .* (data.T_int ./ data.λ), y)`.

Nothing else changes: the demand system, solutions, acceptance gates, and
the returned fields (`M` excludes `T`; `diff = S − (IX − (M + T))`) are
untouched.

## Consequences

- Committed generations are invariant at reported precision: every matrix run
  uses unit supply shocks (`experiments/run.jl` builds
  `Shocks(ones(N), ones(N), zeros(N))`) and the `matrix_5x3_v3` manifests
  record `max|p−1| ≤ 5.5e-13`; with `a = 1`, `k_u` is a weighted power mean
  of prices, so `|k_u−1| ≤ max|p−1|` and the leak total (~0.247 of GDP) moves
  by at most ~2e-13. Recorded `canary_diff` values and the paper flow tables
  are unchanged at every reported precision; no `-v4` generation is minted
  and `matrix_5x3_v3` remains the cited generation (ADR-0004: v1/v2 stay as
  history).
- Supply-side scenarios now close at machine precision (full-71, +20 %
  shock: 2.05e-14 after the fix).
- Regression tests added: exact one-sector root in
  `tests/test_promoted_closures.jl` (pre-fix value 3.0 vs exact 1.5) and a
  guarded full-71 +20 % shock check in `tests/test_calibration.jl` asserting
  the identity closes (`atol = 1e-9`) and that the pre-fix valuation would be
  off by more than 1e-6.
- Open items unchanged: recycling of `T_int` into the government budget
  (ADR-0013); `experiments/run.jl` records `canary_s`/`canary_ixm`/
  `canary_diff` but not a separate `T` component.
