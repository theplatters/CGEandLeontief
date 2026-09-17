# DELTA — IO-type endpoint (Leontief corner)

Implementation: `src/closures/labor/labor.jl` (`delta_elasticities`,
`delta_model`, `leontief_multiplier`); the corner solves through the GAMMA
(`:fixed`) system in `src/core/equilibrium.jl`. Status lives only in
`registry/closures.toml` (`[labor.DELTA]`) and renders on `docs/status.md`
(ADR-0003); this page carries no status.

## Formulation (ADR-0002, incl. the 2026-09-17 amendment)

GAMMA (`w/P = 1`) + Leontief limit of the CES core (`θ, ϵ, σ → 0⁺`) at
`η = 1`. DELTA is a **corner** reached as GAMMA + Leontief technology — the
endpoint row of the 5×3 matrix, not an independent mechanism.

## Design notes (preserved from the frozen `cbase2/src/closures.jl` header)

Under demand-only shocks (A = 1) the zero-profit system pins prices at p = 1
exactly, so relative prices are fixed and quantities absorb demand — the
Robinson (2006) fixed-price multiplier structure. The analytic counterpart
(exact linear Leontief system with the endogenous consumption feedback) is
`leontief_multiplier` and is used as the equivalence test. `delta_elasticities`
solves the limit at a small `ε` for numerical robustness; results are
insensitive to `ε` up to O(ε) under `:fixed` + demand-only shocks.

## What was promoted (Phase 2, ADR-0005)

- `delta_elasticities`, `delta_model` (`:fixed` + Leontief elasticities,
  `financing` kwarg), and `leontief_multiplier` (exact `(I − G)y = b` solve
  with the endogenous consumption feedback, finiteness gate on the gain
  column sums) verbatim, modulo the root-module paths.
- The fixed-wage system it solves is the Phase 2 formulation: all N
  market-clearing equations enforced, `w = 1` numeraire, no CPI pin
  (resolves the omitted-N open item).
- `cbase2/src/calibration.jl` was NOT promoted (Phase 3); the analytic
  counterpart reads the root `Data` v3 absorption fields at their
  compatibility defaults until the Phase 3 calibration lands.

## Remaining gates and caveats

Source: `registry/closures.toml` (`[labor.DELTA]` `open_gates`) and
`cbase2/review.md` §§1, 3.5, 2.9.

- Contract tests for the promoted code land in the next step (kept in
  `open_gates`; status stays `implemented`).
- Type I vs Type II multiplier identification unresolved: the analytic
  counterpart includes endogenous consumption, saving and import leaks (a
  Type II/SAM-style multiplier, not `(I−A)⁻¹`; review §1).
- Exactness rests on p = 1 under demand-only shocks, not on the epsilon
  limit (review §3.5); near-unit round-gain column sums make it fragile.
- The recorded equivalence asserts relative agreement below 5e-3, not
  machine zero (review §2.9).
