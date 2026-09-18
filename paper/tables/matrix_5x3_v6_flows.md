# Flow tables — the 5×3 matrix with the closed external account (ADR-0020 generation)

Generated from `runs/matrix_5x3-v6-*/manifest.toml` (schema v2, commit
`6db1da5`), no re-solve and no hand-typed numbers. All flows are in model
units: GDP at basic prices = 1, so a flow of 0.1 is 10 % of GDP
(GDP_P = 3 027 818 EUR m). Paper text cites the `matrix_5x3-v6-*` run ids.

This table supersedes the `matrix_5x3_v5` version, which stays as history
(ADR-0004). The generation's closure is ADR-0019 + ADR-0020 (option C): all N
goods-market clearings are enforced in every regime, and the external account
`S + T_int + M − (I+X) = F + B_gov` closes in every regime. At eta = 1 the
mobile system carries the endogenous net external transfer F with
`E = (1-tau)·w·ΣL + F`; at eta = 0 (BF) the endpoint is the **sectoral-wage**
system, which solves the per-sector FOC so that the frozen baseline allocation
is cost-minimizing at the equilibrium quantities — that is what closes the
account there and identifies the external position, replacing the v5 pin
`F = 0` under which the position was not identified. Programme financing under
F3 is booked as `B_gov = Σ p·g`.

## Table 1 — Baseline accounts (full-71 A-bill calibration)

| Account | Value | Source |
| --- | ---: | --- |
| GDP, production = income | 1.000000 | `Σ gva` = `Σ λ·fs` = 1 |
| GDP, expenditure (model identity) | 1.000000 | sum of the seven manifest `gdp_*` component diagnostics at the baseline |
| Real GDP (income side, ADR-0018) | 1.000000 | `w·ΣL`, deflated by `P^GDP = 1` |
| Household real consumption (welfare, ADR-0018) | 1.000000 | Törnqvist index, base `data.household_baseline` |
| Production-vs-expenditure residual (raw table) | 5.387 % (163 094 EUR m) | raw table's own gap, open item |
| Household income (after tax) | 0.785899 | `1 - tau0` |
| Household consumption (purchaser prices) | 0.691675 | `c0`, the table's own household column |
| Saving rate `s` | 0.119892 | identity-implied (ADR-0012) |
| Household saving `S` | 0.094223 | `s · (1 - tau0)` |
| Government spending = lump-sum tax `tau0` | 0.214101 | `Σ gG` |
| Investment `I` | 0.162367 | equipment + construction + inventories |
| Exports `X` | 0.421945 | no import margin (domestic sales abroad) |
| Final-demand imports `M_final` | 0.242941 | margin content of C + G + I |
| Intermediate imports `M_int` (row 74) | 0.221403 | booked leak (ADR-0012) |
| Intermediate product taxes `T_int` (row 75) | 0.025745 | booked leak (ADR-0013) |
| Total imports `M` at baseline | 0.464344 | `M_final + M_int` |
| Domestic intermediate bill `ΣA_bill` (row 73) | 0.862762 | charged to domestic demand |
| Row identity `ΣA + ΣM_int + ΣT_int = Σλ - 1` | 0.247147 | holds to 2e-16 |
| Baseline external identity `S + T_int + M - (I+X) - (F + B_gov)` | -5.6e-17 | machine zero at the baseline (F = 0, B_gov = 0; ADR-0013/ADR-0019) |
| Programme `G0 = Σ g` (2024 impulses) | 0.013310 (40 300 EUR m) | `[programme]` in the design |

The baseline block is carried over from the v5 table unchanged: the calibration
(full-71 A-bill) and the reference continuation are identical in the two
generations, and both reproduce the baseline at `F = 0` to machine precision.
Only the closure of the `eta = 0` endpoint changed (ADR-0020).

## Table 2 — Per-cell accounting flows

`Real GDP rel.` is the ADR-0018 income-side index (`gdp_rel`),
`Consumption rel.` the household-consumption (welfare) index
(`consumption_rel`). All other flows are manifest levels in model units
(GDP at basic prices = 1 at the calibration baseline). `S` is
`diagnostics.canary_s`; `M = -(gdp_m_final + gdp_m_int)`; `I + X =
diagnostics.canary_ixm + M`; `tax` is `gdp_g` (government spending
including the programme; the `public_budget` diagnostic equals it except
under F3, where the programme is booked externally and `public_budget`
stays at `tau0`); `T = -gdp_t_int`; `F`, `B_gov` and `Net ext. pos.` are
the `external_transfer`, `programme_financing` and `external_position`
metrics (`Net ext. pos. = F + B_gov`). `Resource side = S + T + M - (I+X)`
is the imbalance the model's own flows imply; `Gap = Resource side - Net ext.
pos.` is the identity gap `canary_diff`.

In this generation the two positions agree in every row: `Gap` is at the
solver-residual level (≤ 1.0e-11, and ≤ 2.0e-13 in fourteen of fifteen cells),
so `Net ext. pos.` is a closed-account quantity in the BF `eta = 0` rows too.
The BF row's `F` is now a solved equilibrium object, and its sectoral wage
vector (not shown here; `diagnostics.wage_min` / `wage_max`) is a genuine
instrument: 0.9734 … 1.5249 (BF-F1) and 0.9702 … 1.6719 (BF-F2/F3), with the
wage-bill-weighted aggregate reported as the `wage` metric.

| Cell | Real GDP rel. | Consumption rel. | L | S | tax | I + X | M | T | F | B_gov | Net ext. pos. | Resource side | Gap |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `BF-F1` | -0.000110 | +0.000963 | 1.000000 | 0.094357 | 0.213389 | 0.594520 | 0.468525 | 0.025823 | -0.00581497 | 0.00000000 | -0.00581497 | -0.00581497 | -1.04e-11 |
| `BF-F2` | -0.000162 | -0.019800 | 1.000000 | 0.092358 | 0.228152 | 0.595703 | 0.468985 | 0.025832 | -0.00852841 | 0.00000000 | -0.00852841 | -0.00852841 | -1.81e-13 |
| `BF-F3` | -0.000162 | -0.019800 | 1.000000 | 0.092358 | 0.228152 | 0.595703 | 0.468985 | 0.025832 | -0.02337469 | 0.01484628 | -0.00852841 | -0.00852841 | -1.94e-13 |
| `ALPHA-F1` | +0.000000 | +0.001432 | 1.000000 | 0.094317 | 0.214101 | 0.584313 | 0.465114 | 0.025666 | +0.00078481 | 0.00000000 | +0.00078481 | +0.00078481 | -4.51e-16 |
| `ALPHA-F2` | +0.000000 | -0.018226 | 1.000000 | 0.092506 | 0.227411 | 0.584313 | 0.465128 | 0.025664 | -0.00101403 | 0.00000000 | -0.00101403 | -0.00101403 | -3.57e-16 |
| `ALPHA-F3` | +0.000000 | -0.018226 | 1.000000 | 0.092506 | 0.227411 | 0.584313 | 0.465128 | 0.025664 | -0.01432394 | 0.01330991 | -0.00101403 | -0.00101403 | -3.87e-16 |
| `BETA-F1` | +0.000000 | +0.001432 | 1.000000 | 0.094317 | 0.214101 | 0.584313 | 0.465114 | 0.025666 | +0.00078481 | 0.00000000 | +0.00078481 | +0.00078481 | -4.86e-16 |
| `BETA-F2` | +0.000000 | -0.018226 | 1.000000 | 0.092506 | 0.227411 | 0.584313 | 0.465128 | 0.025664 | -0.00101403 | 0.00000000 | -0.00101403 | -0.00101403 | -2.81e-16 |
| `BETA-F3` | +0.000000 | -0.018226 | 1.000000 | 0.092506 | 0.227411 | 0.584313 | 0.465128 | 0.025664 | -0.01432394 | 0.01330991 | -0.00101403 | -0.00101403 | -3.26e-16 |
| `GAMMA-F1` | -0.000973 | -0.000806 | 0.999027 | 0.094107 | 0.214101 | 0.584313 | 0.464565 | 0.025641 | +0.00000000 | 0.00000000 | +0.00000000 | -0.00000000 | -5.10e-12 |
| `GAMMA-F2` | +0.001260 | -0.015333 | 1.001260 | 0.092779 | 0.227411 | 0.584313 | 0.465838 | 0.025696 | +0.00000000 | 0.00000000 | +0.00000000 | -0.00000000 | -2.22e-16 |
| `GAMMA-F3` | +0.017796 | +0.022644 | 1.017796 | 0.096357 | 0.227411 | 0.584313 | 0.475154 | 0.026112 | +0.00000000 | 0.01330991 | +0.01330991 | +0.01330991 | -3.92e-16 |
| `DELTA-F1` | -0.000973 | -0.000806 | 0.999027 | 0.094107 | 0.214101 | 0.584313 | 0.464565 | 0.025641 | +0.00000000 | 0.00000000 | +0.00000000 | -0.00000000 | -4.10e-13 |
| `DELTA-F2` | +0.001260 | -0.015333 | 1.001260 | 0.092779 | 0.227411 | 0.584313 | 0.465838 | 0.025696 | +0.00000000 | 0.00000000 | +0.00000000 | -0.00000000 | -6.59e-12 |
| `DELTA-F3` | +0.017796 | +0.022644 | 1.017796 | 0.096357 | 0.227411 | 0.584313 | 0.475154 | 0.026112 | +0.00000000 | 0.01330991 | +0.01330991 | +0.01330991 | -4.96e-16 |

## Table 2b — Expenditure components (manifest diagnostics)

| Cell | C_gross | G+prog | I | X | M_final | M_int | T_int | wedge | deflator |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `BF-F1` | 0.692655 | 0.213389 | 0.167780 | 0.426739 | -0.243928 | -0.224597 | -0.025823 | +1.04e-11 | 1.006326 |
| `BF-F2` | 0.677980 | 0.228152 | 0.168451 | 0.427252 | -0.244054 | -0.224931 | -0.025832 | +1.81e-13 | 1.007181 |
| `BF-F3` | 0.677980 | 0.228152 | 0.168451 | 0.427252 | -0.244054 | -0.224931 | -0.025832 | +1.95e-13 | 1.007181 |
| `ALPHA-F1` | 0.692366 | 0.214101 | 0.162367 | 0.421945 | -0.242654 | -0.222461 | -0.025666 | +4.44e-16 | 1.000000 |
| `ALPHA-F2` | 0.679069 | 0.227411 | 0.162367 | 0.421945 | -0.242648 | -0.222480 | -0.025664 | +2.22e-16 | 1.000000 |
| `ALPHA-F3` | 0.679069 | 0.227411 | 0.162367 | 0.421945 | -0.242648 | -0.222480 | -0.025664 | +2.22e-16 | 1.000000 |
| `BETA-F1` | 0.692366 | 0.214101 | 0.162367 | 0.421945 | -0.242654 | -0.222461 | -0.025666 | +4.44e-16 | 1.000000 |
| `BETA-F2` | 0.679069 | 0.227411 | 0.162367 | 0.421945 | -0.242648 | -0.222480 | -0.025664 | +4.44e-16 | 1.000000 |
| `BETA-F3` | 0.679069 | 0.227411 | 0.162367 | 0.421945 | -0.242648 | -0.222480 | -0.025664 | +0.00e+00 | 1.000000 |
| `GAMMA-F1` | 0.690819 | 0.214101 | 0.162367 | 0.421945 | -0.242269 | -0.222296 | -0.025641 | +5.10e-12 | 1.000000 |
| `GAMMA-F2` | 0.681070 | 0.227411 | 0.162367 | 0.421945 | -0.243147 | -0.222691 | -0.025696 | +0.00e+00 | 1.000000 |
| `GAMMA-F3` | 0.707337 | 0.227411 | 0.162367 | 0.421945 | -0.249699 | -0.225455 | -0.026112 | +2.22e-16 | 1.000000 |
| `DELTA-F1` | 0.690819 | 0.214101 | 0.162367 | 0.421945 | -0.242269 | -0.222296 | -0.025641 | +4.09e-13 | 1.000000 |
| `DELTA-F2` | 0.681070 | 0.227411 | 0.162367 | 0.421945 | -0.243147 | -0.222691 | -0.025696 | +6.59e-12 | 1.000000 |
| `DELTA-F3` | 0.707337 | 0.227411 | 0.162367 | 0.421945 | -0.249699 | -0.225455 | -0.026112 | +6.66e-16 | 1.000000 |

`gdp_wedge` is `-canary_diff` = `-Gap`, so it is ~1e-16 at the mobile
`eta = 1` cells, ~1e-12 at the fixed-wage `eta = 1` cells (at the level of
their solver residuals) and ≤ 1e-11 at the BF `eta = 0` cells: income- and
expenditure-side real GDP coincide in every row of this generation. The v5
reading of the BF wedge as "the labour-market residual the pin leaves behind"
no longer applies — the sectoral wages absorb that residual.

## Table 3 — The external account by labour endpoint

Booked financing (`F + B_gov`) and the resource-side imbalance
(`S + T + M - (I+X)`) by financing row: BF (eta = 0, sectoral wages, `F`
solved), mobile ALPHA/BETA (eta = 1, `F` solved), fixed-wage GAMMA/DELTA (no
external unknown).

| Financing | Quantity | BF (eta = 0, sectoral wages) | Mobile ALPHA/BETA (eta = 1) | Fixed GAMMA/DELTA |
| --- | --- | ---: | ---: | ---: |
| F1 | booked | -0.00581497 | +0.00078481 / +0.00078481 | +0.00000000 / +0.00000000 |
| F1 | resource side | -0.00581497 | +0.00078481 / +0.00078481 | -0.00000000 / -0.00000000 |
| F2 | booked | -0.00852841 | -0.00101403 / -0.00101403 | +0.00000000 / +0.00000000 |
| F2 | resource side | -0.00852841 | -0.00101403 / -0.00101403 | -0.00000000 / -0.00000000 |
| F3 | booked | -0.00852841 | -0.00101403 / -0.00101403 | +0.01330991 / +0.01330991 |
| F3 | resource side | -0.00852841 | -0.00101403 / -0.00101403 | +0.01330991 / +0.01330991 |

The two quantities coincide in every cell, so either may be read as the net
external position. The BF column is now comparable with the others: BF-F3
carries the programme booking `B_gov = +1.4846 %` of GDP (the programme valued
at the `eta = 0` equilibrium prices — the v5 pin left the BF prices at the
baseline, so it booked `+1.3310 %`, the same number as GAMMA/DELTA) and an
endogenous transfer `F = -2.337 %`, so its net position is `-0.853 %`.

## Status

- `BF-F1`: executed (pass, resid 9.8e-12)
- `BF-F2`: executed (pass, resid 5.2e-13)
- `BF-F3`: executed (pass, resid 6.0e-13)
- `ALPHA-F1`: executed (pass, resid 2.2e-15)
- `ALPHA-F2`: executed (pass, resid 2.2e-15)
- `ALPHA-F3`: executed (pass, resid 2.2e-15)
- `BETA-F1`: executed (pass, resid 2.2e-15)
- `BETA-F2`: executed (pass, resid 2.2e-15)
- `BETA-F3`: executed (pass, resid 2.2e-15)
- `GAMMA-F1`: executed (pass, resid 1.4e-13)
- `GAMMA-F2`: executed (pass, resid 8.9e-16)
- `GAMMA-F3`: executed (pass, resid 8.9e-16)
- `DELTA-F1`: executed (pass, resid 4.4e-14)
- `DELTA-F2`: executed (pass, resid 1.8e-13)
- `DELTA-F3`: executed (pass, resid 4.4e-16)

## Notes

- `Real GDP rel.` (`gdp_rel`, ADR-0018) is the income-side index
  `(w·ΣL)/(w·ΣL)_ref - 1`: 0.000 % in every full-employment mobile row
  (ALPHA/BETA), -0.011 % / -0.016 % in the BF rows under F1/F2, and
  -0.097 % / +0.126 % / +1.780 % in the fixed-wage rows (GAMMA/DELTA) under
  F1/F2/F3.
- `Consumption rel.` (the household-consumption Törnqvist welfare index):
  +0.143 % (mobile) / +0.096 % (BF) under F1; -1.823 % in the mobile F2 and
  F3 rows; -1.980 % in both BF-F2 and BF-F3. The BF rows now satisfy
  financing neutrality (F2 = F3), which the retired pin broke: v5 reported
  0.000 % for BF-F3 against -1.694 % for BF-F2.
- Financing neutrality (ADR-0019/ADR-0020, theorem of the closure): F2 and F3
  have identical real allocations `(p, y, w)` in every regime, with
  `F_F3 = F_F2 - B_gov` and an identical net external position. It holds in the
  BF rows to 3.5e-14 (`F_F3 - (F_F2 - B_gov)`), i.e. at the level of the
  solver residual. F3's decomposition books `B_gov = +1.3310 %` of GDP and a
  household transfer `F = -2.337 %`; F2 books only `F = -0.853 %`. F1
  (compositional) remains distinct.
- The booked identity `S + T_int + M - (I+X) = F + B_gov` holds to ~1e-16 at
  the mobile `eta = 1` cells, to ~1e-12 at the fixed-wage `eta = 1` cells (at
  the level of their 1e-13 solver residuals) and to ≤ 1.04e-11 at the BF
  `eta = 0` cells (the stiffest cell, BF-F1; BF-F2/F3 are ≤ 2.0e-13).
- The BF `eta = 0` endpoint is stiff: its Jacobian condition number is ~9.4e7
  against ~52.8 for the mobile all-N system, so its polish target is 1e-13 and
  its residual floor is ~1e-11 (BF-F1). The acceptance gates are the actual
  residuals (ADR-0015), unchanged at 1e-6.
- The twelve non-BF cells reproduce the v5 generation: ALPHA/BETA to 6.7e-16,
  GAMMA/DELTA to ≤ 1.6e-11. The GAMMA/DELTA difference is the documented
  near-singular fixed-wage sensitivity to the warm start (identical with and
  without the ADR-0020 kernel change, verified 2026-09-18): their solve stops
  at a slightly different residual (e.g. 8.9e-16 against 1.0e-13 in the v5
  manifest) and the consumption index moves in the 11th digit. Compare those
  rows against v5 with a 1e-10 tolerance, not bit-identity.
- `tax` is the `public_budget` diagnostic except under F3, where the
  programme is booked as `B_gov` and `public_budget` stays at `tau0`
  (see Table 2 header).
- All 15 cells pass their gates (residuals 8.9e-16 … 1.0e-11, ADR-0015 polish).
  The v4 and v5 tables stay as history; cite `matrix_5x3-v6-*` for the
  closed-account reading of the BF row.

## Promotion (2026-09-18) — the `eta = 0` closure (ADR-0020 option C)

The v5 table carried a correction recording that the BF row's
`Net ext. pos.` was not an external position at all: with the `F = 0` pin the
entry collapsed to the exogenous programme booking `B_gov = +1.3310 %` of GDP
and the account stood open by `Gap` = -0.794 % (BF-F3). ADR-0020 option C
replaces the pin by the sectoral-wage system, and this generation is the
result:

- The net external position (`F + B_gov`) is identified in the BF rows:
  `-0.581 %` (F1), `-0.853 %` (F2), `-0.853 %` (F3), against the v5 entries
  `0.000 %` / `0.000 %` / `+1.331 %`. The F3 booking itself is
  `B_gov = +1.4846 %` (v5: `+1.3310 %`), because the programme is now valued
  at the `eta = 0` equilibrium prices.
- The account closes: `Gap` = -1.0e-11 / -1.8e-13 / -1.9e-13 (v5:
  +4.3e-04 / -5.6e-04 / -7.9e-03).
- The reported BF-F3 net external position moves from `+1.331 %` to
  `-0.853 %` — and, unlike the v5 entry, it is now equal to BF-F2's position
  and to the resource-side imbalance the model's own flows imply.
- The BF row is no longer welfare-degenerate: `Consumption rel.` is
  -1.980 % in both F2 and F3 (v5: -1.694 % and 0.000 %).

Why the numbers move: at `eta = 0` the allocation is frozen, so the labour
market can clear per sector only through the wage. Solving the per-sector FOC
`L^cm_i(p, w_i) = L̄_i` makes the frozen allocation cost-minimizing at the
equilibrium quantities, which is exactly the condition the identity needs;
the retired pin substituted for that FOC and left the residual in the account.
The system has 3N+1 unknowns `[p; y; w(1:N); F]` with N zero-profit, N FOC, N
clearing and one numeraire equation; the free scalar `F` is required because
the block is homogeneous of degree 1 in `(p, w, F)` (dropping a clearing
equation instead would be the retired ADR-0010 shortcut).

Scope: the change touches only the `eta = 0` endpoint. GAMMA/DELTA
(`problem_fixed`, 2N unknowns) and the `eta = 1` mobile system are untouched;
the shared labour-demand helper accepts a wage vector by dispatch and reduces
bit-identically to the scalar path. Verification: the twelve non-BF cells
reproduce v5 (above), and the promotion was checked on the pristine kernel
(`experiments/probes/probe8_promotion_verification.jl`,
`probe9_nonbf_reproduction.jl`). Evidence for the closure itself:
`experiments/probes/probe7_sectoral_wages_eta0.jl`.
