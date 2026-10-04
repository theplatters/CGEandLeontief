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

The `Resource side` and `Gap` columns were added on 2026-09-18 (correction,
see Notes). `Resource side = S + T + M - (I+X)` is the imbalance the model's
own flows imply; `Gap = Resource side - Net ext. pos.` is the identity gap
`canary_diff`. At the twelve `eta = 1` cells the two positions agree to
machine precision, so `Net ext. pos.` is a closed-account quantity there. At
the three BF `eta = 0` cells they do not: `Net ext. pos.` is the *booked*
financing entry `F + B_gov`, which the `F = 0` pin (ADR-0019 D2) reduces to
the programme cost `B_gov`, and the account is open by `Gap`. Read the BF
rows from `Resource side` plus the Notes, never from `Net ext. pos.`.

| Cell | Real GDP rel. | Consumption rel. | L | S | tax | I + X | M | T | F | B_gov | Net ext. pos. | Resource side | Gap |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `BF-F1` | +0.000000 | +0.000433 | 1.000000 | 0.094223 | 0.214101 | 0.584313 | 0.464869 | 0.025655 | +0.00000000 | 0.00000000 | +0.00000000 | +0.00043438 | 4.34e-04 |
| `BF-F2` | +0.000000 | -0.016936 | 1.000000 | 0.092628 | 0.227411 | 0.584313 | 0.465445 | 0.025678 | +0.00000000 | 0.00000000 | +0.00000000 | -0.00056182 | -5.62e-04 |
| `BF-F3` | +0.000000 | +0.000000 | 1.000000 | 0.094223 | 0.227411 | 0.584313 | 0.469599 | 0.025864 | +0.00000000 | 0.01330991 | +0.01330991 | +0.00537381 | -7.94e-03 |
| `ALPHA-F1` | +0.000000 | +0.001432 | 1.000000 | 0.094317 | 0.214101 | 0.584313 | 0.465114 | 0.025666 | +0.00078481 | 0.00000000 | +0.00078481 | +0.00078481 | -2.00e-16 |
| `ALPHA-F2` | +0.000000 | -0.018226 | 1.000000 | 0.092506 | 0.227411 | 0.584313 | 0.465128 | 0.025664 | -0.00101403 | 0.00000000 | -0.00101403 | -0.00101403 | -4.09e-16 |
| `ALPHA-F3` | +0.000000 | -0.018226 | 1.000000 | 0.092506 | 0.227411 | 0.584313 | 0.465128 | 0.025664 | -0.01432394 | 0.01330991 | -0.00101403 | -0.00101403 | -4.35e-16 |
| `BETA-F1` | +0.000000 | +0.001432 | 1.000000 | 0.094317 | 0.214101 | 0.584313 | 0.465114 | 0.025666 | +0.00078481 | 0.00000000 | +0.00078481 | +0.00078481 | -2.05e-16 |
| `BETA-F2` | +0.000000 | -0.018226 | 1.000000 | 0.092506 | 0.227411 | 0.584313 | 0.465128 | 0.025664 | -0.00101403 | 0.00000000 | -0.00101403 | -0.00101403 | -4.13e-16 |
| `BETA-F3` | +0.000000 | -0.018226 | 1.000000 | 0.092506 | 0.227411 | 0.584313 | 0.465128 | 0.025664 | -0.01432394 | 0.01330991 | -0.00101403 | -0.00101403 | -4.37e-16 |
| `GAMMA-F1` | -0.000973 | -0.000806 | 0.999027 | 0.094107 | 0.214101 | 0.584313 | 0.464565 | 0.025641 | +0.00000000 | 0.00000000 | +0.00000000 | -0.00000000 | -1.14e-12 |
| `GAMMA-F2` | +0.001260 | -0.015333 | 1.001260 | 0.092779 | 0.227411 | 0.584313 | 0.465838 | 0.025696 | +0.00000000 | 0.00000000 | +0.00000000 | -0.00000000 | -3.74e-12 |
| `GAMMA-F3` | +0.017796 | +0.022644 | 1.017796 | 0.096357 | 0.227411 | 0.584313 | 0.475154 | 0.026112 | +0.00000000 | 0.01330991 | +0.01330991 | +0.01330991 | -7.07e-12 |
| `DELTA-F1` | -0.000973 | -0.000806 | 0.999027 | 0.094107 | 0.214101 | 0.584313 | 0.464565 | 0.025641 | +0.00000000 | 0.00000000 | +0.00000000 | -0.00000000 | -4.75e-12 |
| `DELTA-F2` | +0.001260 | -0.015333 | 1.001260 | 0.092779 | 0.227411 | 0.584313 | 0.465838 | 0.025696 | +0.00000000 | 0.00000000 | +0.00000000 | -0.00000000 | -6.19e-12 |
| `DELTA-F3` | +0.017796 | +0.022644 | 1.017796 | 0.096357 | 0.227411 | 0.584313 | 0.475154 | 0.026112 | +0.00000000 | 0.01330991 | +0.01330991 | +0.01330991 | -5.03e-16 |

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

`gdp_wedge` is `-canary_diff` = `-Gap`: ~1e-16 at the mobile eta = 1 cells and
~1e-12 (at residual level) at the fixed-wage eta = 1 cells, so income- and
expenditure-side real GDP coincide there; at the BF eta = 0 cells it is the
labour-market residual the `F = 0` pin leaves behind,
`-w * (sum L^cm - sum L)`, where `sum L^cm` is the cost-minimizing aggregate
labour demand at the pinned equilibrium and `sum L` the labour supplied (the
frozen bar at eta = 0, the demand-determined allocation in every eta = 1
regime; see Notes).

## Table 3 -- The external account by labour endpoint

Booked financing (`F + B_gov`) and the resource-side imbalance
(`S + T + M - (I+X)`) by financing row: BF (eta = 0, `F` pinned to 0), mobile
ALPHA/BETA (eta = 1, `F` solved), fixed-wage GAMMA/DELTA (no external unknown).

| Financing | Quantity | BF (eta = 0) | Mobile ALPHA/BETA (eta = 1) | Fixed GAMMA/DELTA |
| --- | --- | ---: | ---: | ---: |
| F1 | booked | +0.00000000 | +0.00078481 / +0.00078481 | +0.00000000 / +0.00000000 |
| F1 | resource side | +0.00043438 | +0.00078481 / +0.00078481 | -0.00000000 / -0.00000000 |
| F2 | booked | +0.00000000 | -0.00101403 / -0.00101403 | +0.00000000 / +0.00000000 |
| F2 | resource side | -0.00056182 | -0.00101403 / -0.00101403 | -0.00000000 / -0.00000000 |
| F3 | booked | +0.01330991 | -0.00101403 / -0.00101403 | +0.01330991 / +0.01330991 |
| F3 | resource side | +0.00537381 | -0.00101403 / -0.00101403 | +0.01330991 / +0.01330991 |

Correction (2026-09-18): in the `eta = 1` columns the two quantities agree to
machine precision, so either may be read as the net external position. In the
BF column they differ by the open account, and the `booked` entry is not an
external position at all: it is the exogenous F3 programme booking
`B_gov = sum p . g` = 1.3310 percent of GDP, identical to the last digit in
BF-F3, GAMMA-F3 and DELTA-F3 because it is a design constant, not a result.
At the eta = 0 endpoint the external position is not identified (the closure
pins `F = 0`; see the Notes).

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
  identical real allocations `(p, y, w)` in the MOBILE regime, with
  `F_F3 = F_F2 - B_gov` and an identical net external position. F3's
  decomposition books `B_gov = +1.331 %` of GDP and a household transfer
  `F = -1.432 %`; F2 books only `F = -0.101 %`. The F2/F3 differences
  reported in v1-v4 were artifacts of the omitted-market shortcut, not
  economics. F1 (compositional) remains distinct. The eta = 0 rows are
  excluded from this reading: the theorem needs a free `F`, and the BF
  endpoint pins it (ADR-0019 D2, decision item 7).
- The BF eta = 0 cells carry the labour-market residual that the `F = 0` pin
  leaves behind (reported, not gated): `Gap` = `canary_diff` = +4.3438e-04
  (F1), -5.6182e-04 (F2), -7.9361e-03 (F3). Zero-profit prices the
  cost-minimizing labour demand, and at the pinned point that demand stands
  at `sum L^cm` = 0.9995656 / 1.0005618 / 1.0079361 against the frozen bar
  `sum L` = 1, so `Gap = -w * (sum L^cm - sum L)`. The frozen sectoral
  allocation is not what creates the gap: the gap closes when the pin is set
  to the labour-clearing value (see the correction below).
- The external position is the booked `F + B_gov`, not the v1-v4 residual
  canary. The fixed-wage rows close through employment (L moves) with no
  external transfer, and their account closes to ~1e-12, so `Net ext. pos.`
  is a genuine position there. At the BF eta = 0 rows the booked entry is
  the exogenous financing entry alone and the account is open by `Gap`: see
  the correction below before reading any BF entry as an external position.
- `tax` is the `public_budget` diagnostic except under F3, where the
  programme is booked as `B_gov` and `public_budget` stays at `tau0`
  (see Table 2 header).
- All 15 cells pass their gates (residuals 4.4e-16..2.0e-13, ADR-0015 polish);
  v5 reproduces the v4 equilibrium values on the fixed rows and BF rows,
  while the mobile F1/F2/F3 rows move with the explicit external account
  (see the v4 table for the superseded decomposition).

## Correction (2026-09-18) -- the BF eta = 0 external position

Audit of the BF-F3 entry (+0.01330991), prompted by its move from -0.007947
in v4. The manifests are reproduced exactly (independently re-solved through
the same harness: every metric of all 15 cells agrees to <= 1e-11, mostly
1e-16), so no arithmetic error is involved. The finding is about what the
entry is.

It is the booked financing, not an equilibrium object. `external_position =
F + B_gov`, and the eta = 0 endpoint pins `F = 0`, so the entry collapses to
`B_gov = sum p . g` = the programme cost = 1.3310 percent of GDP. It is
identical to the last digit in BF-F3, GAMMA-F3 and DELTA-F3 -- a quantity no
labour closure can move is not a labour-market result. The three rows are
not alike either: GAMMA/DELTA-F3 close their account to ~1e-12 (employment
rises to 1.017796 and the resources genuinely arrive from abroad), while
BF-F3 is open by -0.7936 percent of GDP.

Its own resource side says something else. `S + T + M - (I+X)` = +0.5374
percent of GDP at BF-F3, against the booked +1.3310 percent; the difference
is `Gap` = -0.7936 percent, i.e. the labour-market residual of the pin. So
the entry overstates the imbalance the model's flows imply by 0.79 percentage
points -- 60 percent of its own size. The same holds at BF-F1 (booked
+0.0000 vs resource +0.0434) and BF-F2 (booked +0.0000 vs resource -0.0562).

The position is not identified at eta = 0. Replacing the closure's pin by
`F = c` leaves every `c` an exact root (residual 2.2e-15) and moves the
reported entry one-for-one: `c` in {-0.02, -0.01, 0, +0.01, +0.02} gives
-0.669 / +0.331 / +1.331 / +2.331 / +3.331 percent of GDP, and the resource
side moves with it (-0.355 / +0.091 / +0.537 / +0.983 / +1.429 percent). The
`F = 0` pin (ADR-0019 D2) is a normalisation of the immobile benchmark, not
an equilibrium condition on the external account; it substitutes for the
labour-supply row, which is identically zero at eta = 0.

What the pin does and does not carry. Solving the eta = 0 system with the pin
set to ALPHA's solved `F` reproduces the ALPHA cell to 2.8e-17 in the
canonical vector, in household expenditure and in the consumption block;
replacing the pin by the labour equation `sum L^cm = Lbar` closes the account
to ~4e-16 and returns BF-F1/F2/F3 = ALPHA-F1/F2/F3 exactly (+0.00078481 /
-0.00101403 / -0.00101403 percent of GDP), with financing neutrality then
holding at eta = 0 as well (`F_F3 = F_F2 - B_gov` to 4e-17). In this
demand-only design the frozen sectoral allocation enters nothing aggregate
(only `sum L_i`, which equals the bar either way), so the entire BF/ALPHA
difference -- including BF-F3's "welfare unchanged (0.000 percent)" entry,
which becomes -1.8226 percent at the labour-clearing `F` -- is carried by the
pin. The allocation margin needs a design in which it reaches prices or
aggregate demand (sector-specific wages at eta = 0, a supply-side shock, or a
second factor).

Consequences for the paper. Read the BF rows from `Resource side` and `Gap`;
do not read `Net ext. pos.` at the BF rows as an external position, and do
not compare BF-F3's entry with the mobile or fixed-wage entries. ADR-0020
(accepted, option C; executed as the `matrix_5x3_v6` generation, whose flow
table `paper/tables/matrix_5x3_v6_flows.md` carries the current BF row)
records the closure options for eta = 0 and their scope for a v6
generation. Evidence: `experiments/probes/probe5_bf_f3_external_position.jl`,
`probe5b_bf_f3_followups.jl`, `probe5c_bf_pin_equivalence.jl`.
