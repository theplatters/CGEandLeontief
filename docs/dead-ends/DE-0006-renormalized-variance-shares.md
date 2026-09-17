# DE-0006 — Variance shares renormalized over main effects / silent incomplete designs

- **Status:** dead end — the "88.4%" and "85.9%" claims are not variance shares
- **Recorded:** 2026-09-17
- **Origin:** `src/variance_decomposition.jl` (pre-fix), 2026-09
- **Related:** ADR-0003, ADR-0004, `ROADMAP.md` §5, `roadmaps/vertdict.md`
- **Evidence:** `src/variance_decomposition.jl:261` (incomplete design silently
  converted), `src/variance_decomposition.jl:308` (share renormalized over
  main effects), `roadmaps/vertdict.md` ("0.859 would mean 85.9% of total
  variance … whereas 88.4% is merely η's share of the sum of reported main
  effects")

## What was tried

Compute elasticity "variance shares" from a balanced grid of solves, then
report each factor's percentage after renormalizing over the reported main
effects; incomplete designs were completed implicitly.

## Why it seemed promising

The headline "η explains 88.4%" is exactly the kind of quantitative result
the paper wants, and the balanced grid makes the arithmetic look like an
ANOVA decomposition.

## How it failed

1. The percentage was normalized over the **sum of main effects**, not over
   total variance, so it cannot be quoted as a share of variance.
2. Missing factorial cells were silently converted into a nonorthogonal
   design instead of raising an error.
3. Solver failures were discarded, so the design may be incomplete without
   notice.
4. Independently, the underlying mobile-labour results were invalid
   (DE-0001–DE-0003), so even a correct decomposition of them was moot.

## What replaced it

`SobolResult` (with deprecated `VarianceDecompositionResult` alias) and
`eta_sweep_diagnostics`: explicit factor levels and probability weights,
first-order and total-effect indices, selected interactions, solver status
and residuals for every evaluation, and the assumed product measure. A
missing required factorial cell aborts the decomposition; failed regions are
reported as part of the feasible parameter domain. `ROADMAP.md` §5 forbids
renormalizing over main effects.

## Revival conditions

Only as a separately named unbalanced regression/ANOVA with an explicitly
chosen sum-of-squares convention, and never quoted as a variance share.
