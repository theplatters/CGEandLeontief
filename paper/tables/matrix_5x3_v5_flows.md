# Flow tables — the 5×3 matrix with the explicit external account (ADR-0019 generation)

Generated from `runs/matrix_5x3-v5-*/manifest.toml` (schema v2, commit
`f0af583`), no re-solve and no hand-typed numbers. All flows are in model
units: GDP at basic prices = 1, so a flow of 0.1 is 10 % of GDP
(GDP_P = 3 027 818 EUR m). Paper text cites the `matrix_5x3-v5-*` run ids.

This table supersedes the `matrix_5x3_v4` version. The v4 table stays as
history (ADR-0004). The generation's closure is ADR-0019: all N
goods-market clearings are enforced in every regime; the mobile eta = 1
system carries the endogenous net external transfer F with
E = (1-tau)w.SigmaL + F; the eta = 0 endpoint pins F = 0; programme
financing under F3 is booked as B_gov = Sigma p.g.

## Table 1 — Baseline accounts (full-71 A-bill calibration)

| Account | Value | Source |
| --- | ---: | --- |
| GDP, production = income | 1.000000 | `Sigma gva` = `Sigma lambda.fs` = 1 |
| GDP, expenditure (model identity) | 1.000000 | sum of the seven manifest `gdp_*` component diagnostics at the baseline |
| Real GDP (income side, ADR-0018) | 1.000000 | `w·SigmaL`, deflated by `P^GDP = 1` |
| Household real consumption (welfare, ADR-0018) | 1.000000 | Törnqvist index, base `data.household_baseline` |
| Production-vs-expenditure residual (raw table) | 5.387 % (163 094 EUR m) | raw table's own gap, open item |
| Household income (after tax) | 0.785899 | `1 - tau0` |
| Household consumption (purchaser prices) | 0.691675 | `c0`, the table's own household column |
| Saving rate `s` | 0.119892 | identity-implied (ADR-0012) |
| Household saving `S` | 0.094223 | `s . (1 - tau0)` |
| Government spending = lump-sum tax `tau0` | 0.214101 | `Sigma gG` |
| Investment `I` | 0.162367 | equipment + construction + inventories |
| Exports `X` | 0.421945 | no import margin (domestic sales abroad) |
| Final-demand imports `M_final` | 0.242941 | margin content of C + G + I |
| Intermediate imports `M_int` (row 74) | 0.221403 | booked leak (ADR-0012) |
| Intermediate product taxes `T_int` (row 75) | 0.025745 | booked leak (ADR-0013) |
| Total imports `M` at baseline | 0.464344 | `M_final + M_int` |
| Domestic intermediate bill `SigmaA_bill` (row 73) | 0.862762 | charged to domestic demand |
| Row identity `SigmaA + SigmaM_int + SigmaT_int = Sigmalambda - 1` | 0.247147 | holds to 2e-16 |
| Baseline external identity `S + T_int + M - (I+X) - (F + B_gov)` | -5.6e-17 | machine zero at the baseline (F = 0, B_gov = 0; ADR-0013/ADR-0019) |
| Programme `G0 = Sigma g` (2024 impulses) | 0.013310 (40 300 EUR m) | `[programme]` in the design |

The baseline block is unchanged from v4: the calibration (full-71 A-bill)
is identical, and ADR-0019 reproduces it at F = 0 to machine precision, so
only the external-identity line is reworded to the booked form.

## Table 2 -- Per-cell accounting flows

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
metrics (`Net ext. pos. = F + B_gov`).

| Cell | Real GDP rel. | Consumption rel. | L | S | tax | I + X | M | T | F | B_gov | Net ext. pos. |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `BF-F1` | +0.000000 | +0.000433 | 1.000000 | 0.094223 | 0.214101 | 0.584313 | 0.464869 | 0.025655 | +0.00000000 | 0.00000000 | +0.00000000 |
| `BF-F2` | +0.000000 | -0.016936 | 1.000000 | 0.092628 | 0.227411 | 0.584313 | 0.465445 | 0.025678 | +0.00000000 | 0.00000000 | +0.00000000 |
| `BF-F3` | +0.000000 | +0.000000 | 1.000000 | 0.094223 | 0.227411 | 0.584313 | 0.469599 | 0.025864 | +0.00000000 | 0.01330991 | +0.01330991 |
| `ALPHA-F1` | +0.000000 | +0.001432 | 1.000000 | 0.094317 | 0.214101 | 0.584313 | 0.465114 | 0.025666 | +0.00078481 | 0.00000000 | +0.00078481 |
| `ALPHA-F2` | +0.000000 | -0.018226 | 1.000000 | 0.092506 | 0.227411 | 0.584313 | 0.465128 | 0.025664 | -0.00101403 | 0.00000000 | -0.00101403 |
| `ALPHA-F3` | +0.000000 | -0.018226 | 1.000000 | 0.092506 | 0.227411 | 0.584313 | 0.465128 | 0.025664 | -0.01432394 | 0.01330991 | -0.00101403 |
| `BETA-F1` | +0.000000 | +0.001432 | 1.000000 | 0.094317 | 0.214101 | 0.584313 | 0.465114 | 0.025666 | +0.00078481 | 0.00000000 | +0.00078481 |
| `BETA-F2` | +0.000000 | -0.018226 | 1.000000 | 0.092506 | 0.227411 | 0.584313 | 0.465128 | 0.025664 | -0.00101403 | 0.00000000 | -0.00101403 |
| `BETA-F3` | +0.000000 | -0.018226 | 1.000000 | 0.092506 | 0.227411 | 0.584313 | 0.465128 | 0.025664 | -0.01432394 | 0.01330991 | -0.00101403 |
| `GAMMA-F1` | -0.000973 | -0.000806 | 0.999027 | 0.094107 | 0.214101 | 0.584313 | 0.464565 | 0.025641 | +0.00000000 | 0.00000000 | +0.00000000 |
| `GAMMA-F2` | +0.001260 | -0.015333 | 1.001260 | 0.092779 | 0.227411 | 0.584313 | 0.465838 | 0.025696 | +0.00000000 | 0.00000000 | +0.00000000 |
| `GAMMA-F3` | +0.017796 | +0.022644 | 1.017796 | 0.096357 | 0.227411 | 0.584313 | 0.475154 | 0.026112 | +0.00000000 | 0.01330991 | +0.01330991 |
| `DELTA-F1` | -0.000973 | -0.000806 | 0.999027 | 0.094107 | 0.214101 | 0.584313 | 0.464565 | 0.025641 | +0.00000000 | 0.00000000 | +0.00000000 |
| `DELTA-F2` | +0.001260 | -0.015333 | 1.001260 | 0.092779 | 0.227411 | 0.584313 | 0.465838 | 0.025696 | +0.00000000 | 0.00000000 | +0.00000000 |
| `DELTA-F3` | +0.017796 | +0.022644 | 1.017796 | 0.096357 | 0.227411 | 0.584313 | 0.475154 | 0.026112 | +0.00000000 | 0.01330991 | +0.01330991 |

## Table 2b -- Expenditure components (manifest diagnostics)

| Cell | C_gross | G+prog | I | X | M_final | M_int | T_int | wedge | deflator |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `BF-F1` | 0.691675 | 0.214101 | 0.162367 | 0.421945 | -0.242482 | -0.222387 | -0.025655 | -4.34e-04 | 1.000000 |
| `BF-F2` | 0.679961 | 0.227411 | 0.162367 | 0.421945 | -0.242871 | -0.222574 | -0.025678 | 5.62e-04 | 1.000000 |
| `BF-F3` | 0.691675 | 0.227411 | 0.162367 | 0.421945 | -0.245792 | -0.223807 | -0.025864 | 7.94e-03 | 1.000000 |
| `ALPHA-F1` | 0.692366 | 0.214101 | 0.162367 | 0.421945 | -0.242654 | -0.222461 | -0.025666 | -2.22e-16 | 1.000000 |
| `ALPHA-F2` | 0.679069 | 0.227411 | 0.162367 | 0.421945 | -0.242648 | -0.222480 | -0.025664 | 0.00e+00 | 1.000000 |
| `ALPHA-F3` | 0.679069 | 0.227411 | 0.162367 | 0.421945 | -0.242648 | -0.222480 | -0.025664 | 0.00e+00 | 1.000000 |
| `BETA-F1` | 0.692366 | 0.214101 | 0.162367 | 0.421945 | -0.242654 | -0.222461 | -0.025666 | -2.22e-16 | 1.000000 |
| `BETA-F2` | 0.679069 | 0.227411 | 0.162367 | 0.421945 | -0.242648 | -0.222480 | -0.025664 | 0.00e+00 | 1.000000 |
| `BETA-F3` | 0.679069 | 0.227411 | 0.162367 | 0.421945 | -0.242648 | -0.222480 | -0.025664 | 0.00e+00 | 1.000000 |
| `GAMMA-F1` | 0.690819 | 0.214101 | 0.162367 | 0.421945 | -0.242269 | -0.222296 | -0.025641 | 1.14e-12 | 1.000000 |
| `GAMMA-F2` | 0.681070 | 0.227411 | 0.162367 | 0.421945 | -0.243147 | -0.222691 | -0.025696 | 3.74e-12 | 1.000000 |
| `GAMMA-F3` | 0.707337 | 0.227411 | 0.162367 | 0.421945 | -0.249699 | -0.225455 | -0.026112 | 7.07e-12 | 1.000000 |
| `DELTA-F1` | 0.690819 | 0.214101 | 0.162367 | 0.421945 | -0.242269 | -0.222296 | -0.025641 | 4.75e-12 | 1.000000 |
| `DELTA-F2` | 0.681070 | 0.227411 | 0.162367 | 0.421945 | -0.243147 | -0.222691 | -0.025696 | 6.19e-12 | 1.000000 |
| `DELTA-F3` | 0.707337 | 0.227411 | 0.162367 | 0.421945 | -0.249699 | -0.225455 | -0.026112 | 6.66e-16 | 1.000000 |

`gdp_wedge` is `-canary_diff`: ~1e-16 at the mobile eta = 1 cells and ~1e-12
(at residual level) at the fixed-wage eta = 1 cells, so income- and
expenditure-side real GDP coincide there; at the BF eta = 0 cells it is the
documented factor-market gap (see Notes).

## Table 3 -- The external account by labour endpoint

Net external position by financing row: BF (eta = 0, position = B_gov,
F pinned to 0), mobile ALPHA/BETA (eta = 1, position = F + B_gov),
fixed-wage GAMMA/DELTA (position = B_gov, no external unknown).

| Financing | BF (eta = 0) | Mobile ALPHA/BETA (eta = 1) | Fixed GAMMA/DELTA |
| --- | ---: | ---: | ---: |
| F1 | +0.00000000 | +0.00078481 / +0.00078481 | +0.00000000 / +0.00000000 |
| F2 | +0.00000000 | -0.00101403 / -0.00101403 | +0.00000000 / +0.00000000 |
| F3 | +0.01330991 | -0.00101403 / -0.00101403 | +0.01330991 / +0.01330991 |

## Status

- `BF-F1`: executed (pass, resid 6.7e-16)
- `BF-F2`: executed (pass, resid 6.7e-16)
- `BF-F3`: executed (pass, resid 6.7e-16)
- `ALPHA-F1`: executed (pass, resid 6.7e-16)
- `ALPHA-F2`: executed (pass, resid 6.7e-16)
- `ALPHA-F3`: executed (pass, resid 6.7e-16)
- `BETA-F1`: executed (pass, resid 6.7e-16)
- `BETA-F2`: executed (pass, resid 6.7e-16)
- `BETA-F3`: executed (pass, resid 6.7e-16)
- `GAMMA-F1`: executed (pass, resid 6.7e-14)
- `GAMMA-F2`: executed (pass, resid 1.0e-13)
- `GAMMA-F3`: executed (pass, resid 2.0e-13)
- `DELTA-F1`: executed (pass, resid 1.3e-13)
- `DELTA-F2`: executed (pass, resid 1.7e-13)
- `DELTA-F3`: executed (pass, resid 4.4e-16)

## Notes

- `Real GDP rel.` (`gdp_rel`, ADR-0018) is the income-side index
  `(w.SigmaL)/(w.SigmaL)_ref - 1`: 0.000 % in every full-employment mobile row
  (BF/ALPHA/BETA), -0.097 % / +0.126 % / +1.780 % in the fixed-wage rows
  (GAMMA/DELTA) under F1/F2/F3. `Consumption rel.` is the household-consumption
  (welfare) Törnqvist index: +0.143 % (mobile) / +0.043 % (BF) under F1;
  -1.823 % in the mobile F2 and F3 rows, -1.694 % in BF-F2 and 0.000 % in
  BF-F3 (the eta = 0 endpoint pins F = 0, so its F2/F3 rows differ); +2.264 %
  in the fixed-wage F3 rows.
- The booked identity `S + T_int + M - (I+X) = F + B_gov` holds to
  ~1e-16 at the mobile eta = 1 cells and to ~1e-12 at the fixed-wage eta = 1
  cells (at the level of their 1e-13 solver residuals; worst over all eta = 1
  cells is 7.1e-12).
- Financing neutrality (ADR-0019, theorem of the closure): F2 and F3 have
  identical real allocations `(p, y, w)` in the mobile regime, with
  `F_F3 = F_F2 - B_gov` and an identical net external position. F3's
  decomposition books `B_gov = +1.331 %` of GDP and a household transfer
  `F = -1.432 %`; F2 books only `F = -0.101 %`. The F2/F3 differences
  reported in v1-v4 were artifacts of the omitted-market shortcut, not
  economics. F1 (compositional) remains distinct.
- The BF eta = 0 cells carry the documented factor-market gap (reported,
  not gated): `canary_diff` = +4.3438e-04 (F1), -5.6182e-04 (F2),
  -7.9361e-03 (F3). Zero-profit prices the cost-minimizing labour
  demand, not the frozen baseline allocation; all N clearings are enforced
  with the F = 0 pin.
- The external position is now the booked `F + B_gov`, not the v1-v4
  residual canary. The fixed-wage rows close through employment (L moves)
  with no external transfer.
- `tax` is the `public_budget` diagnostic except under F3, where the
  programme is booked as `B_gov` and `public_budget` stays at `tau0`
  (see Table 2 header).
- All 15 cells pass their gates (residuals 4.4e-16..2.0e-13, ADR-0015 polish);
  v5 reproduces the v4 equilibrium values on the fixed rows and BF rows,
  while the mobile F1/F2/F3 rows move with the explicit external account
  (see the v4 table for the superseded decomposition).
