# ADR-0018 — Real GDP measurement in the open economy

- **Status:** accepted (user instruction, 2026-09-18)
- **Date:** 2026-09-18
- **Supersedes:** —
- **Related:** ADR-0003 (registry single source of truth), ADR-0004 (runs are immutable manifests), ADR-0010 (omitted N-th market residual and the canary), ADR-0013 (row-75 `T_int` leak), ADR-0016 (CES-consistent leak valuation), ADR-0017; ROADMAP.md §4.1 ("Use the shared Törnqvist/Divisia helper for real GDP in every new model variant"; "Report real GDP separately from welfare"); `src/core/equilibrium.jl` (`Solution.real_gdp`, both solve paths), `src/core/diagnostics.jl` (new measurement API), `experiments/run.jl` (manifest schema v2), `tests/test_gdp_measurement.jl`; `paper/main.tex` line ~419 (Törnqvist approximation to the Divisia quantity index of final demand)

## Context

The kernel's `real_gdp` accessor and the `Solution.real_gdp` field hold the
Törnqvist quantity index of the gross household consumption block against
`data.household_baseline` (`src/core/equilibrium.jl`, both the fixed-wage and
the mobile solve paths; ported from cbase2). That is the Baqaee–Farhi
constant-returns final-demand/welfare aggregator: in the closed B&F economy it
coincides with real GDP, but the BeyondHulten open economy also carries
government demand, investment, exports, final-import margins, intermediate
imports `M_int` (ADR-0012), product taxes on intermediate use `T_int`
(ADR-0013), and external financing (F3). ROADMAP.md §4.1 therefore forbids
calling a welfare index GDP without proving equivalence.

Verified at HEAD on the full-71 A-bill calibration (scratch checks, not
committed): at any solved equilibrium the model's value added satisfies
`w·ΣL = Σ_j V_j + p·ρ`, where the seven signed aggregate components are
`V = [C_gross, G+programme, I, X, −M_final, −M_int, −T_int]`, with
`C_gross = c_dom/(1−m)`, `M_final` the final-import margin content of C+G+I
(including the programme), and `M_int`, `T_int` valued with the ADR-0016 CES
bill factor `k = p^ε·a^(ε−1)·P^(1−ε)`. At the calibration baseline
`V = [0.691675, 0.214101, 0.162367, 0.421945, −0.242941, −0.221403, −0.025745]`
and `ΣV = w·ΣL = 1.000000` to 2.2e-16. `p·ρ` is the omitted N-th market
residual, i.e. the canary external position (ADR-0010): zero at the baseline
and at fixed-wage cells, nonzero at mobile η = 1 cells (measured: +0.060 % of
GDP for ALPHA-F2, +0.850 % for ALPHA-F3).

Under the demand-only 5×3 matrix all prices are pinned at one, so the GDP
deflator is exactly one and real GDP growth equals the growth of `w·ΣL`.

## Decision

1. Real GDP is the income-side measure: `real GDP = (w·ΣL) / P^GDP`,
   chain-linked to the reference solution, where `P^GDP` is the Törnqvist price
   index of the seven components above (unit values `|V_j|/Q_j`, signed value
   shares `V_j/ΣV`). At the matrix prices `P^GDP = 1`, so real GDP is the
   `w·ΣL` index (in the demand-only matrix, the employment index at the sticky
   wage).
2. The expenditure-side Divisia index `(ΣV_t/ΣV_0)/P^GDP` is the dual and is
   recorded as a cross-check: equal to the income measure at the baseline and
   at fixed-wage cells, differing by the external wedge `ΣV − w·ΣL = −p·ρ` at
   mobile η = 1 cells. The wedge is reported as a diagnostic and is never
   gated to zero.
3. Conventions: GDP at basic prices, so `T_int` remains a booked external leak
   (ADR-0013); imports are valued at domestic basic prices with share margins
   as the model does; the seven quantities
   `Q = [ΣC_gross, Σ(G+programme), ΣI, ΣX, ΣM_final, ΣM_int, ΣT_int]` are all
   strictly positive on the real calibration, and the programme is aggregated
   into government demand so no zero-base component arises. For components
   whose value is zero in one of the two equilibria, the index treats the
   component's unit value as unchanged (no measured price change) and drops
   components that are zero in both; this keeps the smoke fixtures finite.
4. Kernel API (additive, no silent semantic change): new exported functions
   `gdp_components`, `gdp_deflator`, `gdp_income`, `gdp_expenditure`,
   `gdp_wedge`; `real_consumption(sol)` is the canonical name of the old
   consumption index. The legacy `real_gdp`/`Solution.real_gdp` accessor and
   field keep their current value (household-consumption index, the welfare
   aggregator) and are documented as legacy; the paper and manifests no longer
   use them as GDP.
5. Manifests (`runs/<id>/manifest.toml`) move to schema version 2: `[metrics]`
   records `gdp`, `gdp_rel`, `gdp_expenditure`, `gdp_expenditure_rel`,
   `gdp_deflator`, `gdp_wedge`, `consumption`, `consumption_rel`, `employment`,
   `wage`, `nominal_gdp`, `max_abs_price_dev`; `real_gdp*` keys are retired.
   Existing v1 manifests are immutable history.
6. Because this changes `src/`, the provenance of the v3 generation is
   invalidated: the paper-facing generation is minted as `matrix_5x3_v4` (same
   design cells and parameters; v1/v2/v3 stay as history). Paper tables and
   DOCS_ASSESSMENT report real GDP and welfare as separate columns and cite
   the v4 run ids.

## Consequences

- The F2 mobile reading changes interpretation: income GDP is flat (0.000 %)
  while the consumption/welfare index falls 1.69 %; the pre-registered
  "aggregate ~ 0" signature is recovered for GDP. GAMMA/DELTA-F3: income GDP
  +1.78 % versus consumption +2.26 %. The v3 table's "Real GDP rel." column is
  a welfare column and must not be cited as GDP.
- Provenance: the v4 generation is required before any paper claim; v3 stays
  visible (ADR-0004).

## Enforcement

- `tests/test_gdp_measurement.jl` pins the identity `|ΣV − w·ΣL| ≤ 1e-10` on
  the real calibration, the wedge/canary identity at mobile η = 1, the
  deflator inertness under demand-only shocks, and a supply-shock case that
  exercises the deflator.
- `experiments/README.md` documents the v2 manifest metric contract;
  `tests/test_run_manifest.jl` requires the new keys.
