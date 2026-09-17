# ADR-0012 — A-bill calibration: the domestic intermediate bill replaces the purchaser-price bill

- **Status:** accepted
- **Date:** 2026-09-17
- **Supersedes:** the clamp-based household-residual calibration inside the v3
  open-economy block (the "clamped in 31 sectors, mass 0.0858" behaviour)
- **Related:** ADR-0010 (retained-sector rebuild; sanctioned N−1+CPI mobile
  form), ADR-0011 (notebook conformance), review.md findings 2.4/2.5 and §3.2,
  `registry/preregistration.toml` (`designs.calibration_abill`), the
  external-balance canary (`src/core/equilibrium.jl`)

## Context

The v3 open-economy calibration derived the household baseline block `c0` as a
**residual** of the sector accounts:

```
c0_dom[i] = λ[i] − Σ_u Ω_raw[u,i]·(1−fs[u])·λ[u] − (1−m[i])·(gG+inv)[i] − expo[i]
```

The intermediate-demand term multiplied the domestic technology shares Ω_raw
(which are built from the *domestic* product flows, rows 1–71 of the raw use
table) by the sector bill `(1−fs[u])·λ[u]`. Measured against the raw table,
that bill is exactly the **purchaser-price intermediate bill**:

```
(1−fs[u])·λ[u] ≡ A[u] + Imp[u] + Tx[u]        (row 73 + row 74 + row 75;
                                               identity verified to 1.4e-17)
```

with `A[u]` the **domestic** intermediate flow (row 73, "Gesamte Verwendung
der inländischen Produktion"), `Imp[u]` the imported intermediates (row 74,
"Verwendung der Importe", 0.2214 of GDP on full-71), and `Tx[u]` the product
taxes (row 75).

Charging the total bill against domestic production asks the domestic sectors
to supply goods they never produced: the imported and taxed input content is
neither domestic output nor household consumption. The shortfall fell on the
household residual, which went **negative in 31 sectors (0.0858 GDP)**; the
clamp zeroed it and thereby

1. destroyed the exactness of the baseline (a uniform clearing hole of
   clamp-mass/N ≈ 4.3e-4 per market — the measured "solver floor"), and
2. inflated the calibrated saving rate to s ≈ 0.40 (the clamp *removed*
   household demand mass, so Σc0 fell relative to the exact-bill value).

The per-market hole also explains the historical convergence floors
(2.1e-4 / 3.6e-4 / 5.2e-4 across formulations): no solver can remove an
inconsistency that the calibration itself embeds.

## Decision

1. **Intermediate demand is charged with the domestic bill.** In
   `_mobile_market_demand`, `problem_fixed`, and `leontief_multiplier` the
   coefficient `(1−fs[u])` on `y[u]` in the intermediary-demand term is
   replaced by `a[u] = A_bill[u]/λ[u]`, with `A_bill[u] = row73[u]/GDP_P` read
   from the (possibly rebuilt) raw table. The imported+taxed content leaves
   the demand system entirely.
2. **The imported intermediates become an explicit external-account leak.**
   `M_int[u] = row74[u]/GDP_P` scales with sectoral output and enters the
   external-balance canary on the M side: `S = I + X − M` with
   `M = M_cons + M_inj + M_prog + M_int`.
3. **The household block is the table's own domestic household final demand.**
   By the row identity, the residual `c0` under the A-bill equals the observed
   household column; there is nothing to clamp. Only floating-point dust is
   cleaned (threshold 1e-5; full-71 dust ≈ 4e-19; 70s carries one microscopic
   −2.4e-6 from the retained-economy identity gap).
4. **Ω_dom is retired** from the calibration path: Ω_raw is already the
   domestic technology (the shares are built from domestic flows); the old
   Ω_dom construction was a mis-built copy.
5. **The saving rate is re-anchored** to the identity-implied value
   (pre-registered in `registry/preregistration.toml`):
   **s = 0.1199 (full-71), 0.1285 (70s)**. The household expenditure share of
   income rises to 1 − s ≈ 0.88; its import content is carried by the
   final-demand margins and the external account.

## Consequences

- **Baseline root test passes at machine zero**: max|resid| = 3.5e-18 at
  (p, y, w) = (1, λ, 1) on full-71 — the calibration IS an equilibrium.
- **The solver floor is gone through the standing `solve()` path**: exo=1
  3.5e-18 (0.8 s), exo=0.8 4.4e-16, F2 financed programme 2.2e-15 (2.5 s).
  No special solver machinery is required; the exo-continuation's bottom rung
  is obsolete (the target economy solves directly).
- **Negative implied saving rates inside the exo range** (s ≈ −0.07 at
  exo = 0.8) are admissible: shrinking injections raise the residual household
  demand above income, and the excess is externally financed. Never clamped.
- **The "reduced" variant fails its identity gap loudly** (mass 7.7e-3 in 10
  sectors ≈ the 18 bn retained-economy gap) and is deferred together with the
  dataset-variant decision; the 70s variant passes with the documented
  microscopic clamp.
- **The external-identity canary still shows −2.6e-2 at baseline** on full-71:
  this is the raw table's own income-vs-expenditure reconciliation gap (the
  ~3.3 % item), now *visible* instead of hidden inside the clamp. It is
  recorded, not threshold-fitted; booking it as an explicit discrepancy item
  is the pending decision.
- **All downstream headline numbers are re-anchored**: w* ≈ 1 at the
  numeraire, household consumption 0.88 of income (0.69 GDP, of which ≈ 0.24
  import content), employment neutrality of the financed programme confirmed
  (F2: L = 1.000000 at machine residuals).
- The wedge remains retired (ADR-0010): only the endpoints η ∈ {0, 1} carry
  the storyline.
