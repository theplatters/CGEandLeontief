# DE-0010 — The cbase2 multi-stage solver ladder

- **Status:** dead end — superseded by the root kernel's single `solve()` path
- **Recorded:** 2026-09-17
- **Origin:** `cbase2/src/solvers.jl` (frozen snapshot, ADR-0001), Stages 1–2
- **Related:** ADR-0006 (experiment entry point), ADR-0010 (N−1 + CPI mobile
  form), ADR-0012 (A-bill calibration), `registry/freeze.toml` open items
- **Evidence:** `runs/matrix_5x3-*/manifest.toml` (13 `executed` cells,
  2026-09-17, commit `848a9cb`), `runs/index.csv`

## What was tried

A staged solver ladder for the mobile system: an `exo_scale` bisection to the
smallest scale with `s ≥ 0`, then a ladder of `exo_scale` steps × a θ ladder
(2.0 → 0.5), with per-rung fallbacks — residual-gated Levenberg–Marquardt
polish on the 142-dimensional system, jitter restarts, trust-region retries,
and an `η_s`-homotopy for BETA. The design assumed the equilibrium could only
be reached by creeping along a continuation from a well-behaved corner.

## Why it seemed necessary

Three solver generations failed to get below a residual floor of
2.1e-4 / 3.6e-4 / 5.2e-4, the θ = 1, k = 0 continuation stalled at 3.1e-4, the
F3 mobile solve stalled at ~2e-4, and the BETA `η_s`-continuation ran >25 min
without completing. A ladder with fallbacks was the natural response.

## How it failed

The floor was not a solver property. ADR-0012 showed it was the calibration's
own clamp hole spread over the markets (clamp mass / N = 4.3e-4), and ADR-0013
removed the last accounting gap. With the A-bill calibration the standing
`solve()` reaches machine tolerance directly: 11 of the 13 executed matrix
cells land between 2.2e-15 and 8.9e-16, the remaining three
(`BF-F1` 9.27e-7, `DELTA-F2` 8.71e-7, `GAMMA-F2` 4.85e-7) pass the 1e-6 gate
without any fallback, and the whole reference continuation plus 15 cells runs
in about two minutes. No jitter, trust-region retry or homotopy rung was
needed anywhere.

## What replaced it

The root kernel's `solve()` — Newton with a residual-gated LM polish — plus the
residual external-account canary as the acceptance test. The `exo_scale`
bisection and the θ ladder remain in `experiments/run.jl` as a *warm-start
continuation* (they produce the reference solution and the `real_gdp` anchor),
which is a different thing from a solver ladder: the ladder's claim was that
the *target* economy needs the staging, and that claim is dead.

## Conditions for revival

A calibration that genuinely embeds an inconsistency (e.g. a retained-economy
variant with an unresolved identity gap larger than the gate) would bring the
floor back. In that case the ladder is the wrong remedy — fix the calibration,
as ADR-0012 and ADR-0013 did.
