# ADR-0014 — real-wage elastic supply and a verified scale-determinacy guard

- **Status:** accepted (user rulings, 2026-09-17)
- **Date:** 2026-09-17
- **Supersedes:** the `has_additive_anchor` admissibility heuristic for the
  fixed-wage η = 1 system; the raw-wage form of the elastic-supply residual
- **Related:** DE-0004 (nominal-wage supply is a dead end; the real wage is
  required), DE-0010 (the solver ladder is retired), ADR-0010 (N−1 + CPI mobile
  form), ADR-0012/ADR-0013 (A-bill calibration), `designs/matrix_5x3_v2.toml`,
  `docs/log/2026-09.md` (2026-09-17 entries)

## Context

Two findings from the post-matrix probes (full-71 A-bill data, no code change):

**1. The elastic-supply residual was nominal in form.** DE-0004 records
nominal-wage supply as a dead end and requires the real wage,
`L^s = L̄ · [(w/P)/(w₀/P₀)]^{η_s}`. The implementation used the raw `w`, on the
stated grounds that the CPI numeraire makes `P = 1` — true in the current gauge,
but it leaves the closure's *meaning* numeraire-dependent and re-opens DE-0004
the moment the numeraire changes.

**2. The fixed-wage η = 1 guard was a heuristic.** `_solve_fixed` rejected the
cell whenever `!has_additive_anchor(financing)`, inherited from the cbase2
**closed-economy** unit root (column sums exactly 1 at MPS = 1, κ = fs, no
margins). It never computed a rank. Measured on the A-bill open economy:

| check | result |
| --- | --- |
| Jacobian at the F1 point (fixed, η = 1) | full rank, σ_min/σ_max = **0.1586**, no singular value < 1e-8·σ_max |
| F1 system solved from λ, 2λ, λ/2 (guard bypassed) | all three converge to the **same** root: resid 4.4e-16, L = 0.9990271532561, max\|y − λ\| = 0.0063279836, agreeing to 1e-15 |

So `GAMMA-F1` and `DELTA-F1` were executable; the guard's premise did not hold.

## Decision

1. **The elastic-supply residual is deflated by the CPI.**
   `labor_market_residual(c::ElasticLaborClosure, model, L_sum, w, cpi)` computes
   `L̄ · [ (w/P) / (w₀/P₀) ]^{η_s}` with `P = cpi`; the anchor `c.w0` is the
   baseline **real** wage. The ALPHA method takes the same argument (unused) for
   a uniform interface, the mobile system passes its own CPI, and the
   experiment gate deflates identically so the gate tests the same equation.
2. **Scale determinacy is verified, not assumed.** The fixed-wage η = 1 guard
   now tests the round-gain criterion: the clearing block is `y = G·y + const`
   with column sums

   ```
   colsum_u = A_bill_u/λ_u + (1 − m_u)·(1 − s)·fs_u,
   ```

   so `(I − G)` is singular — a unit root, hence a continuum of solutions —
   exactly when `max(colsums)` reaches 1. The guard fires only then. This is the
   same quantity as the finiteness gate of `recalibrate_open`.
3. **A new generation of the matrix is run** as `designs/matrix_5x3_v2.toml`
   (15 cells, `matrix_5x3-v2-<L>-<F>`). The v1 runs stay as history (ADR-0004);
   paper tables cite the v2 ids.

## Consequences

- **The criterion discriminates exactly as intended**: on the closed
  `tiny_fixture` (m = s = 0, A_bill = (1−fs)λ) `max(colsums)` = **1.0** and the
  guard fires (the contract test's expected failure is unchanged); on full-71 it
  is **0.8559345674224728** and the system is determinate.
- **`GAMMA-F1` and `DELTA-F1` now execute** (resid 9.96e-7 and 9.41e-7, L =
  0.9990416 / 0.9990154) — the matrix is complete at 15 cells, with no locked
  cells and no dropped scenarios.
- **The CPI deflation is numerically inert where the solution is at machine
  precision**: `BETA-F2` reproduces its recorded `real_gdp_rel` to 5.6e-16 and
  `ALPHA-F2` to 2.2e-16.
- **Open item (new, measured):** cells that sit *at* the solver's 1e-6
  acceptance threshold have tolerance-dependent metrics. `GAMMA-F2` (resid
  4.85e-7) moves by **8.3e-7** in `real_gdp_rel` under a 1e-16 perturbation of
  the residual. Five cells are affected (BF-F1, GAMMA-F1, DELTA-F1, GAMMA-F2,
  DELTA-F2, all 4.9e-7…9.96e-7). Tightening the polish target from 1e-6 to
  ~1e-10 would push them to machine precision; it is **proposed, not done**
  (the v1 lesson "don't tighten tolerances" was about a calibration that
  embedded an inconsistency, which ADR-0012/0013 removed — so it no longer
  binds, but the change is a separate decision).
- **BETA stays uninformative for demand-only shocks.** The deflation does not
  change that: `η_s = 0.5 / 2 / 5` give solutions identical to ~1e-15 under the
  programme, because the price block is demand-invariant (`p ≡ 1`) so the real
  wage sits at its anchor. `η_s` bites only with a supply-side scenario
  (+20 % sector-1 shock: L = 1.00245 at `η_s = 0.5` vs 1.00983 at `η_s = 2.0`).
  Recorded in `registry/closures.toml` (BETA open gates).
