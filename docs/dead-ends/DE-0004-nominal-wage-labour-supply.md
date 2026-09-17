# DE-0004 — Nominal wage in the labour-supply function

- **Status:** dead end — violates numeraire homogeneity
- **Recorded:** 2026-09-17
- **Origin:** rejected mobile-labour system, 2026-09
- **Related:** ADR-0002, closures `BF` and `BETA`; `ROADMAP.md` §4.2
- **Evidence:** `roadmaps/vertdict.md` item 4 ("the model uses nominal `w`
  rather than `w/P` in labor supply"), `ROADMAP.md` §3 and §4.2

## What was tried

Labour supply as a function of the nominal wage, `L = L̄ w^{η}`.

## Why it seemed promising

The estimated elasticity literature often regresses employment on nominal
wages, and the numeraire was pinned anyway, so the distinction looked
harmless.

## How it failed

Changing the nominal unit changes real allocations: the model fails the
numeraire-homogeneity requirement. It also confounds the wage regime with the
price level and invalidates comparisons across closures with different
numeraires.

## What replaced it

Real-wage supply everywhere: `L = L̄ ((w/P)/(w₀/P₀))^{η_s}` for `BETA`
(`cbase2/src/closures.jl`), and the fixed-*real*-wage `GAMMA` peg.
`ROADMAP.md` §4.2 requires a homogeneity test ("changing the nominal unit
must not alter real allocations") for every closure.

## Revival conditions

None. A nominal-wage rule can only appear as an explicitly named nominal
rigidity (for example a one-sided nominal floor) with its own closure id,
calibration, and homogeneity discussion.
