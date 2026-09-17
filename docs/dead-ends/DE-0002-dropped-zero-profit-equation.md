# DE-0002 — Dropping a zero-profit equation to close the system

- **Status:** dead end — do not retry
- **Recorded:** 2026-09-17
- **Origin:** rejected mobile-labour system (`src/mobile_labor.jl` pre-fix), 2026-09
- **Related:** ADR-0002, DE-0001, closure `BF`, `ROADMAP.md` §4.2
- **Evidence:** `src/mobile_labor.jl:158` (pre-fix), `roadmaps/vertdict.md`;
  sector-1 profit residual ≈ 0.39 at the "solution"

## What was tried

After adding the labour-market equation and a numeraire, drop a **zero-profit**
equation to restore a square system.

## Why it seemed promising

The equation count balances, the solver converges, and the missing condition
is invisible unless residuals are checked explicitly.

## How it failed

Walras' law permits dropping a redundant **goods-market** equation, not a
zero-profit condition. The system therefore solved a different model, with a
large profit residual in the "dropped" sector, while reporting convergence
(the solver return code was not checked either).

## What replaced it

The square system documented in `roadmaps/vertdict.md` and `ROADMAP.md`:
`N` zero-profit equations, `N−1` goods-market equations, one labour-closure
equation, one numeraire equation — plus the post-solve requirement to compute
**all** `N` goods-market residuals and assert them, repeat the solve while
omitting each market in turn (omitted-equation invariance), and check the
solver termination status and all residuals.

## Revival conditions

None for the original form. If an equation is ever omitted for numerical
reasons, the omitted-equation invariance gate of `ROADMAP.md` Phase 4 must
show real allocations are invariant to which market is omitted and that every
omitted residual is small.
