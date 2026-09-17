# ADR-0013 — the intermediate-bill tax term: row 75 booked as an external leak

- **Status:** accepted (ratified by the user 2026-09-17)
- **Date:** 2026-09-17
- **Corrects:** the ADR-0012 reading of the −2.6e-2 canary gap as "the raw
  table's own income-vs-expenditure reconciliation gap"; the gap is the
  unbooked third component of the purchaser-price intermediate bill
- **Related:** ADR-0010 (N−1+CPI clearing and the asserted external-account
  canary, review finding 2.1), ADR-0012 (A-bill calibration),
  `registry/preregistration.toml` (`designs.calibration_abill`),
  `tests/test_calibration.jl`, `tests/test_promoted_closures.jl`

## Context

ADR-0012 replaced the purchaser-price intermediate bill with the **domestic**
bill in the demand system. The raw table's row identity is

```
(1 − fs_u)·λ_u ≡ A_u + Imp_u + Tx_u      (rows 73 + 74 + 75)
```

where `A_u` is the domestic intermediate use (row 73), `Imp_u` the imported
intermediates (row 74) and `Tx_u` the product taxes on intermediate use
(row 75). ADR-0012 moved `Imp` out of the demand system into the external
account as `M_int`, but **`Tx` was booked nowhere**: it was charged neither to
domestic demand (correct — taxes demand no goods) nor to the external account.

The consequence is visible in the omitted-market canary. At the exact baseline
root of the full-71 A-bill calibration the omitted N-th market residual is
machine zero (7.8e-18), while the canary's `S − (I+X−M)` reads −2.574494e-2 —
so the review-2.1 identity, which ADR-0010 asserts at 1e-6, cannot hold.

Measured on the full-71 table (`GDP_P = 3 027 818` EUR m):

| term | value (GDP_P units) | source |
| --- | --- | --- |
| `Σ A_bill` (row 73) | 0.8627622267917028 | domestic intermediate bill |
| `Σ M_int` (row 74) | 0.2214030037472530 | imported intermediates |
| `Σ T_int` (row 75) | **0.0257449423974625** | product taxes on intermediate use |
| `Σ λ − 1` | 0.2471471729364182 | = the two non-domestic components |

`Σ T_int` equals the canary gap to 2e-16: **0.0257449423974625** vs the
measured −0.0257449423974625. The gap is therefore not a data discrepancy; it
is an omitted accounting term. It is also not the table's own
production-vs-expenditure residual (5.387 %, i.e. 0.053877 of GDP), which
remains a separate, unrelated open item.

Before the fix the gap appeared *near*-constant across financing closures
(F1 +2.5659e-2, F2 +2.5673e-2, F3 +2.5783e-2) — a drift of 1.2e-4 — because
the unbooked term scales with sectoral output and prices, exactly like
`M_int`. A constant "discrepancy import" could therefore never have closed the
identity without an unjustified tolerance.

## Decision

1. **`Tx` is booked as the third explicit leak of the intermediate bill.** The
   `Data` struct gains `T_int[u] = row75[u] / GDP_P` (per-user product taxes on
   intermediate use, model units), extracted by `_domestic_bills` alongside
   `A_bill` and `M_int` and carried through `assemble_data` and
   `recalibrate_open`.
2. **`external_balance_canary` adds the valued tax term**
   `T = dot(p .* (data.T_int ./ data.λ), y)` and returns it; the reported
   imbalance becomes `S − (I + X − M) + T`. The identity asserted for mobile
   η = 1 solutions is unchanged in form and tolerance.
3. **Nothing else changes.** `T_int` is not part of the demand system (it
   appears in neither `_mobile_market_demand` nor `problem_fixed`), so prices,
   quantities, employment, `s` and the household block are untouched. The
   saving rate stays `s = 0.1199` (full-71) / `0.1285` (70s).
4. **The government-side treatment is left open.** Product taxes are
   government revenue; the model's government is balanced by a lump-sum tax
   `τ0 = Σ gG`. Recycling `T_int` into that budget (reducing `τ0`, or funding
   `G`) is a modelling extension, not needed for the identity, and is recorded
   as an open item rather than decided here.

## Consequences

- **The canary identity is exact again**, at machine precision rather than at
  the 1e-6 gate, on both variants: full-71 baseline 7.7e-17 and solved root
  7.7e-17; 70s baseline 9.6e-17 and solved root 6.3e-16.
- **The 70s "microscopic clamp" (−2.4e-6) is the retained-economy version of
  the same term**, not a data defect: with `T_int` booked, the 70s baseline
  omitted-market residual (−2.3775e-6) and the canary (−2.3775e-6) agree to
  1e-16, so the disclosed residual is now visible on both sides of the
  identity instead of hidden in the assertion's slack.
- **The "discrepancy-import" decision is closed**: there was no discrepancy to
  book. ADR-0012's phrasing is corrected; the `designs/calibration_abill.toml`
  gate note is updated accordingly.
- **The reference continuation of the `matrix_5x3` batch no longer aborts.**
  The θ ladder always solved; the batch died in `assert_external_canary`
  afterwards. With `T_int` booked, `experiments/run.jl --design matrix_5x3`
  proceeds to the cells.
- **Tests:** `tests/test_calibration.jl` pins `Σ T_int` (0.0257 full-71 /
  0.0260 70s), the row identity `Σ A_bill + Σ M_int + Σ T_int = Σ λ − 1`, and
  the canary identity at the baseline and solved roots on both variants;
  `tests/test_promoted_closures.jl` checks that the tax term enters the canary
  linearly on the synthetic fixture (which runs without the real table).
- The raw table's 5.387 % production-vs-expenditure residual and the
  dataset-variant decision (`reduced` deferred, `70s` documented) stay as
  recorded open items.
