# ADR-0038 — The wage-pinned closure is an employment-adjustment closure

- **Status:** accepted
- **Date:** 2026-10-23
- **Supersedes:** the closure name and the real-wage formulation recorded for
  GAMMA and DELTA (ADR-0002's taxonomy wording; `registry/closures.toml`
  `[labor.GAMMA]`, `[labor.DELTA]`). Stable ids, kernel code, recorded runs and
  every real outcome are unchanged.
- **Related:** ADR-0037 (consumption-unit price reporting), ADR-0017, ADR-0019,
  ADR-0020, ADR-0029; `OpenCGE@src/core/equilibrium.jl` (`problem_fixed`,
  `_solve_fixed`), `OpenCGE@src/closures/labor/types.jl`.

## Context

The manuscript described GAMMA as a two-sided **real**-wage peg,
`w/P(p) = wbar = 1`, with employment endogenous and uncapped, and DELTA as its
Leontief corner. The kernel implements something else: `problem_fixed` solves a
2N system in `(p, y)` with `w = 1` hard-coded as the numeraire and `F ≡ 0`, and
records employment as a post-solve outcome. No real-wage equation is imposed,
because none can be.

Corollary `cor:realwage` already gave the reason: with a single common wage the
zero-profit block is homogeneous of degree one in `(p, w)` and contains no demand
variable, so `p = w·π(A)` and

    w/P(p) = w/(w·P(π)) = 1/P(π(A)),

a function of technology alone. The real wage is therefore a solution object,
not an instrument. A peg on it holds if and only if `P(π(A)) = 1` — the
calibration point — and nowhere else; appended to the wage-pinned system it is
one equation on a scale that the price block has already fixed. The old wording
was thus not merely a description that drifted from the code: as written, the
closure was not implementable in this model class.

The two defects surfaced together because they share one cause. Under
demand-only shocks relative prices do not move, `P(π) = 1`, and both readings
coincide — which is why the demand-only matrix was unaffected. Open the pricing
doors (ADR-0029) and the CPI moves between -0.9 % and +1.9 % across the executed
cells; the
reported deflators then differed between wage and consumption units (ADR-0037),
and the real wage visibly left its "peg" in exactly the same cells. The units
error and the real-wage deviation are the same number to first order:
`ΔP_wage − ΔP_cons ≈ 100(w/P − 1)`.

## Decision

1. **Rename the closure.** GAMMA is the **employment-adjustment closure**
   (wage-pinned): the labor-supply equation is dropped, the money wage is held at
   its baseline value as the numeraire, the net external transfer is fixed at
   zero, and employment is endogenous, uncapped and demand-determined. DELTA is
   its fixed-coefficient corner, not an independent labor-market mechanism. The
   ids GAMMA and DELTA and the label `eq:closure-gamma`/`eq:closure-delta` are
   retained so no cross-reference or run id moves.
2. **Record the formulation as solved**, not as described: `w = 1` (numeraire),
   `F = 0`, `L_i = L_i^cm(p, y_i, w)`, `η = 1`, all N goods markets clearing, no
   aggregate employment constraint, no wage-setting rule.
3. **State the identity in the manuscript** as Corollary `cor:peg`, with Remark
   `rem:empadj` naming what the closure implements. `w = 1` is a unit
   convention, not a wage rigidity: a static real model has no monetary anchor
   for a sticky-money-wage mechanism (Corollary `cor:nogo`), and it is not a
   wage floor either, since no employment constraint is imposed.
4. **Do not interpret any result as holding the real wage constant.** Under the
   pricing doors GAMMA employment responses come with a moving real wage; the
   demand-only coincidence with a real-wage-anchored closure is a property of the
   slice, not of the closure.
5. **Scope the nesting claims to the slice.** GAMMA is the `η_s → ∞` endpoint of
   the elastic family on the demand-only branch only. In the single-wage member
   the flat limit is empty off that branch (the anchor is technological); the
   sectoral family owns the N wages that make an anchor attainable in principle,
   so there the limit is unattained rather than empty. Figure 1 and Table A1
   carry this asymmetry.
6. **Keep the two corrections separate.** ADR-0037's consumption-unit reporting
   is necessary and unchanged: renormalization cannot repair a moving real wage,
   and the rename cannot repair mixed units. All GAMMA/DELTA pricing rows stay in
   consumption units, where the real-wage movement is visible.

## Consequences

No number moves: real outcomes are invariant to the nominal scale, so no run
regeneration, pin bump or kernel edit is involved. The label change reaches the
manuscript (§1, §3, §4, §5, §6, App A, App B, App H), nine table files and two
generators, and the R1.5 / R2 responses. The interpretation change is confined to
what the closure is *called* and to the claim that workers' purchasing power is
held fixed — which the revision no longer makes. The FIDELIO/E3ME mapping
survives as a wage-setting-with-quantity-adjustment reading, not as a fixed
real wage.

The load on Corollary `cor:nogo` increases: with the peg gone, the categorical
difference between the IO and CGE poles is which equation is absent, not a
parameter value, and that is the corollary's claim.

## Enforcement

- `grep -rn "fixed real wage\|fixed-wage\|real wage pegged" revision/` returns
  hits only in comments and in the retired-name history.
- `julia --project=. tests/runtests.jl` in the kernel repo is unaffected (no code
  change); `scripts/check_repo.jl` and `scripts/status.jl` in the paper repo.
- `python3 revision/tables/make_table4_results.py` and
  `python3 revision/tables/make_pricing_door.py` regenerate their tables from the
  same inputs; the only diff against the previous generation is the row label.
- `python3 scripts/check-manuscript-refs.py`; `latexmk -pdf` in a temporary copy.
