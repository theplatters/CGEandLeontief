# ADR-0015 — monotone polish to ~1e-10 (tolerance-independent metrics)

- **Status:** accepted (user instruction, 2026-09-17)
- **Date:** 2026-09-17
- **Related:** ADR-0014 (which measured the problem), DE-0010 (ladder retired),
  `designs/matrix_5x3_v3.toml`

## Context

The solver accepted a cell as soon as its residual fell inside the gate
(1e-6 fixed / 1e-5 mobile) and only polished *beyond* that when it did not.
Five matrix cells therefore stopped at 4.9e-7 … 9.96e-7, and their reported
metrics were functions of the solver's stopping point rather than of the model:
`GAMMA-F2` (resid 4.85e-7) moved **8.3e-7** in `real_gdp_rel` under a 1e-16
perturbation of the residual (ADR-0014, measured).

## Decision

The residual-gated Levenberg–Marquardt polish in both solve paths
(`_solve_fixed` and the mobile `solve`) now triggers below **1e-10** instead of
the acceptance gate, runs with `reltol = abstol = 1e-12`, and is **monotone**: a
polish step is kept only if it strictly improves the residual, so polishing can
never turn a passing cell into a failing one. The acceptance gates are
unchanged.

## Consequences

- The five affected cells now sit at machine precision: `GAMMA-F2` 4.85e-7 →
  **1.04e-13**, `GAMMA-F1` 9.96e-7 → **1.95e-13**, `DELTA-F2` 8.71e-7 →
  **2.2e-16**, `BF-F1` 9.27e-7 → **2.2e-15**; cells already at machine
  precision are untouched (`ALPHA-F2` reproduces its v2 `real_gdp` to the last
  digit).
- Their *metrics* move, though: `GAMMA-F1`'s `real_gdp` shifts by 1.8e-5 and
  `GAMMA-F2`'s by 5.4e-6. The v2 generation's numbers for those cells are
  therefore superseded — hence `designs/matrix_5x3_v3.toml`, which the paper
  cites. v1 and v2 stay as history (ADR-0004).
- A side effect worth having: the GAMMA ≡ DELTA limit relation now holds to
  **5e-12** in `real_gdp` (it was 1e-9, i.e. limited by the solver, not the
  model).
- The earlier lesson "do not tighten tolerances" was about a calibration that
  embedded an inconsistency (ADR-0012/0013 removed it); with a consistent
  calibration the polish converges, so the lesson no longer binds.
