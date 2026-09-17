# DE-0003 — `η = 10⁶` as a fixed-real-wage approximation

- **Status:** dead end — the approximation does not have the claimed limit
- **Recorded:** 2026-09-17
- **Origin:** high-elasticity runs of the mobile-labour model, 2026-09
- **Related:** ADR-0002, closures `BETA` and `GAMMA`; `ROADMAP.md` §5 Closure C
- **Evidence:** `roadmaps/vertdict.md` ("the high-η limit was misinterpreted as
  unlimited labor instead of L → 0"; employment collapses as η grows),
  `ROADMAP.md` §3

## What was tried

Approximate the fixed-real-wage closure by letting the supply elasticity
`η_s → 10⁶` in `L = L̄ (w/P)^{η_s}`, expecting the quantity to become
unconstrained.

## Why it seemed promising

The textbook limit of a very elastic supply curve is a horizontal supply
curve; numerically one parameter goes to a large number.

## How it failed

At the equilibrium used, `w/P < 1`, so `L = L̄ (w/P)^η → 0` as `η` grows, not
`∞`. Employment collapses instead of being unconstrained. The reported
"high-η" results mixed two different regimes and produced an artificial
"price invariance with changing employment" signature — the red flag recorded
in `ROADMAP.md` Phase 4.

## What replaced it

An explicit **GAMMA** fixed-real-wage closure (`w/P = w̄`, employment
endogenous), and, for a bounded/unemployment version, a genuine
complementarity formulation (`ZETA`:
`0 ≤ L̄ − L ⊥ w/P − ω̄ ≥ 0`). `ROADMAP.md` Phase 3 explicitly forbids
`η = 10⁶` as the endpoint.

## Revival conditions

Only as a numerical continuation device inside a closure whose limit is
proven analytically, and never as the definition of the fixed-wage closure.
