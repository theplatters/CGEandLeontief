# ADR-0005 — Phase 2 layout and cbase2 backport

- **Status:** accepted
- **Date:** 2026-09-17
- **Supersedes:** —
- **Related:** ADR-0001, ADR-0002, ADR-0003, `registry/closures.toml`,
  `registry/freeze.toml`, `docs/closures/beta.md`, `docs/closures/delta.md`,
  `docs/closures/financing.md`, `src/core/`, `src/closures/`

## Context

Phase 1 left the root kernel in flat `src/*.jl` files while BETA, DELTA and
F1–F3 lived only in the frozen `cbase2/` snapshot, with two recorded kernel
DIFFs (`src/core/interface.jl`, `src/core/mobile_labor.jl`) and the v3 open
items (Törnqvist base, F3 stall, omitted N-th market, θ continuation). The
characterization suite (`tests/test_kernel_regression.jl`, commit `8dd9907`)
pins pre-migration mobile behavior with golden numbers.

## Decision

**Layout.** The canonical kernel is now:

- `src/core/accounting.jl` (types, `Data`, `Shocks`, `Model`, data I/O),
- `src/core/technology.jl` (CES + Leontief + Cobb-Douglas),
- `src/core/equilibrium.jl` (`Solution` + mobile-labor equilibrium),
- `src/core/diagnostics.jl` (utilities, impulses, variance decomposition),
- `src/closures/labor/types.jl` (closure descriptions),
- `src/closures/labor/labor.jl` (BETA/DELTA mechanics),
- `src/closures/financing/financing.jl` (F1/F2/F3),
- `src/closures/registry.jl` (id → constructor; `ZETA` maps to `nothing`),
- `src/plots.jl` (unchanged stubs).

Include order is accounting → technology → labor/types → equilibrium →
labor/labor → financing → registry → diagnostics → plots. Status lives only
in `registry/closures.toml`; `src/closures/registry.jl` carries no status.

**Backport.** The recorded cbase2 DIFFs are ported onto the relocated kernel
with these compatibility guarantees:

- `Data` gains the six v3 absorption fields with cbase2's names/order; all
  legacy constructors default them to zeros EXCEPT `household_baseline`,
  which defaults to the legacy Törnqvist base
  (`consumption_share .* Σ labor_share`) — NOT cbase2's zeros — so legacy
  `real_gdp` behavior is preserved; the Phase 3 calibration overwrites it.
- `Model` gains `financing::AbstractFinancing`; the 3-arg constructor
  defaults to `NoFinancing()`.
- `MobileLaborCESElasticities` gains `eta_s` (4-arg defaults to 0.0);
  `:beta` symbol; `labor_closure` returns `ElasticLaborClosure(η_s)`.
- Demand in `problem`/`problem_fixed` is financed household demand +
  programme additive demand + baseline `gov`/`exo`/`exports` demand PLUS the
  legacy manna terms A/G. cbase2 retired manna; the root kernel keeps it as
  an explicitly documented compatibility path so the characterization
  goldens and `rerun_results.jl` stay reproducible. With `NoFinancing` and
  zero v3 fields the demand is numerically identical to the pre-Phase-2 code.
- `_intermediate_price` (θ = 1 CD limit) and `_ces_unit_cost` (ϵ = 1 CD
  limit) guards ported; residual-gated bounded LM polish ported (mobile
  1e-5, fixed 1e-6) with the "already exact at init" fast path and the
  `"did not converge"` error surface kept; solver-internal failures are
  reported as non-convergence.
- `equilibrium_residuals`: fixed accepts 2N or 2N+1, mobile requires 2N+1.
- `real_gdp` base is `data.household_baseline` (legacy default above); the
  non-negative-consumption NaN guard is kept.
- `mobile_labor_model` gains `financing`/`eta_s` kwargs (`eta_s` forces
  `:beta`); all existing forms keep working.

**Intentional fixed-wage change.** `problem_fixed` adopts the cbase2 system:
all N market-clearing equations enforced, `w = 1` numeraire, no CPI pin; the
η ≈ 1 anchor guard keys on `has_additive_anchor(financing)` OR legacy
manna. The guard keeps both legacy substrings (`"scale-indeterminate"` and
`"autonomous or investment"`). This resolves the omitted-N open item; fixed
levels move by design while the closure contracts stay green.

**Retirements and scope.** `cbase2/scripts/diff_kernel.jl` retires as a
checker (frozen snapshots keep their recorded parent commits as history; the
file stays untouched in `cbase2/`). `cbase2/src/calibration.jl` and the
cbase2 scripts are NOT promoted in Phase 2 (Phase 3).

## Consequences

- One kernel with closure plugins; BETA/DELTA/F1–F3 run in the root package.
- The Phase 1 open items Törnqvist base, F3 stall (gates/polish), omitted N
  (fixed formulation), and θ continuation (CD guards) are resolved as kernel
  changes; genuine open gates (income effect, w0 anchor, Type-I/II,
  budget-pricing, external-balance scope, contract tests, calibration/
  experiments) stay in `registry/closures.toml`.
- Promoted closures stay `implemented` until the next step's contract tests
  land (then `tested`).

## Enforcement

- `registry/closures.toml` `symbols` must appear in the listed `files`
  (`scripts/status.jl` greps them); `files` point at the new
  `src/closures/...` paths.
- `registry/freeze.toml` records the backport note under `[frozen.cbase2]`
  without changing hashes/commits.
- `tests/test_kernel_regression.jl` stays green and untouched; any
  `tests/test_fixed_closure.jl` expectation change requires justification
  here and in the log (none was needed).
