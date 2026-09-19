# Flow tables — the 5x3 matrix and the sectoral family (ADR-0022 generation)

Generated from `runs/matrix_5x3-v9-*/manifest.toml` (schema v2, commit `048c8fa`),
no re-solve and no hand-typed numbers: every figure below is a manifest field or a
sum of manifest fields, with the mapping stated in the table notes. All flows are in
model units: GDP at basic prices = 1, so a flow of 0.1 is 10 % of GDP
(GDP_P = 3 027 818 EUR m). Paper text cites the `matrix_5x3-v9-*` run ids.

This generation supersedes `matrix_5x3_v7` and `matrix_5x3_v8` (both aborted
mid-batch by promotion-plumbing defects, retained as records) and, through them,
`matrix_5x3_v6` (ADR-0004). Its closure is ADR-0019 + ADR-0020 option C as in v6,
plus ADR-0022: a BETA cell carrying an elasticity vector solves the N sectoral
labour markets as the 3N+1 system. With `eta_s,i = 0` for every sector that system
*is* the v6 `BF` endpoint (proved in `paper/equivalence.tex` v3, Proposition 3), so
no separate rigid sectoral cell is reported: the corner is read from the `BF` cells.

## Table 1 — Baseline accounts

Unchanged from `paper/tables/matrix_5x3_v6_flows.md`, Table 1: the calibration
(full-71 A-bill, ADR-0012/ADR-0013) and the reference continuation are identical in
the two generations, and the fifteen matrix cells of this generation reproduce the
v6 cells bit-for-bit (Table 2). The baseline block is therefore not duplicated here.

## Table 2 — The fifteen matrix cells: bit-identity with v6

The `Delta` column is `max |v9 - v6|` over the six headline metrics of that cell
(income-side real GDP, household-consumption welfare index, employment, GDP
deflator, `max abs(p-1)`, external transfer) plus the `wage_min`/`wage_max` and
`canary_diff` diagnostics. It is zero in every cell: the ADR-0022 extension leaves
the executed matrix untouched, which is why the v6 flow decomposition (Tables 2,
2b and 3 of `matrix_5x3_v6_flows.md`) remains the citation for these rows.

| Cell | Real GDP rel. | Consumption rel. | L | deflator | `max abs(p-1)` | Delta |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `BF-F1` | -0.000110 | +0.000963 | 1.000000 | 1.006326 | 0.221442 | +0.00e+00 |
| `BF-F2` | -0.000162 | -0.019800 | 1.000000 | 1.007181 | 0.279009 | +0.00e+00 |
| `BF-F3` | -0.000162 | -0.019800 | 1.000000 | 1.007181 | 0.279009 | +0.00e+00 |
| `ALPHA-F1` | +0.000000 | +0.001432 | 1.000000 | 1.000000 | 0.000000 | +0.00e+00 |
| `ALPHA-F2` | +0.000000 | -0.018226 | 1.000000 | 1.000000 | 0.000000 | +0.00e+00 |
| `ALPHA-F3` | +0.000000 | -0.018226 | 1.000000 | 1.000000 | 0.000000 | +0.00e+00 |
| `BETA-F1` | +0.000000 | +0.001432 | 1.000000 | 1.000000 | 0.000000 | +0.00e+00 |
| `BETA-F2` | +0.000000 | -0.018226 | 1.000000 | 1.000000 | 0.000000 | +0.00e+00 |
| `BETA-F3` | +0.000000 | -0.018226 | 1.000000 | 1.000000 | 0.000000 | +0.00e+00 |
| `GAMMA-F1` | -0.000973 | -0.000806 | 0.999027 | 1.000000 | 0.000000 | +0.00e+00 |
| `GAMMA-F2` | +0.001260 | -0.015333 | 1.001260 | 1.000000 | 0.000000 | +0.00e+00 |
| `GAMMA-F3` | +0.017796 | +0.022644 | 1.017796 | 1.000000 | 0.000000 | +0.00e+00 |
| `DELTA-F1` | -0.000973 | -0.000806 | 0.999027 | 1.000000 | 0.000000 | +0.00e+00 |
| `DELTA-F2` | +0.001260 | -0.015333 | 1.001260 | 1.000000 | 0.000000 | +0.00e+00 |
| `DELTA-F3` | +0.017796 | +0.022644 | 1.017796 | 1.000000 | 0.000000 | +0.00e+00 |

## Table 3 — Headline metrics, the eighteen sectoral cells

`wage max/min` is `diagnostics.wage_max / wage_min`; the reported `wage` metric is
the wage-bill-weighted aggregate (`sum_i w_i L_i / sum_i L_i`). F2 and F3 agree in
every sectoral cell to all printed decimals (financing neutrality, ADR-0019).

| Cell | Real GDP rel. | Consumption rel. | L | deflator | `max abs(p-1)` | wage max/min |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `BETA-F1-etas025` | -0.007086 | +0.003252 | 1.000809 | 1.003776 | 0.133610 | 1.331995 |
| `BETA-F2-etas025` | -0.007671 | -0.016888 | 1.000864 | 1.004096 | 0.152897 | 1.381559 |
| `BETA-F3-etas025` | -0.007671 | -0.016888 | 1.000864 | 1.004096 | 0.152897 | 1.381559 |
| `BETA-F1-etas05` | -0.004414 | +0.004227 | 1.001191 | 1.002701 | 0.095839 | 1.234983 |
| `BETA-F2-etas05` | -0.004690 | -0.015712 | 1.001255 | 1.002872 | 0.105666 | 1.259605 |
| `BETA-F3-etas05` | -0.004690 | -0.015712 | 1.001255 | 1.002872 | 0.105666 | 1.259605 |
| `BETA-F1-etas1` | -0.001981 | +0.005119 | 1.001564 | 1.001726 | 0.061283 | 1.148364 |
| `BETA-F2-etas1` | -0.002069 | -0.014669 | 1.001627 | 1.001802 | 0.065422 | 1.158491 |
| `BETA-F3-etas1` | -0.002069 | -0.014669 | 1.001627 | 1.001802 | 0.065422 | 1.158491 |
| `BETA-F1-etas2` | -0.000184 | +0.005784 | 1.001859 | 1.001003 | 0.035629 | 1.085418 |
| `BETA-F2-etas2` | -0.000189 | -0.013912 | 1.001914 | 1.001033 | 0.037171 | 1.089127 |
| `BETA-F3-etas2` | -0.000189 | -0.013912 | 1.001914 | 1.001033 | 0.037171 | 1.089127 |
| `BETA-F1-rigidprog` | -0.012839 | -0.005416 | 0.996999 | 1.005456 | 0.215061 | 1.544801 |
| `BETA-F2-rigidprog` | -0.014525 | -0.026914 | 0.996621 | 1.006206 | 0.272033 | 1.697924 |
| `BETA-F3-rigidprog` | -0.014525 | -0.026914 | 0.996621 | 1.006206 | 0.272033 | 1.697924 |
| `BETA-F1-rigidhalf` | -0.007748 | +0.001496 | 1.000188 | 1.005646 | 0.109583 | 1.247681 |
| `BETA-F2-rigidhalf` | -0.008652 | -0.019053 | 1.000196 | 1.006294 | 0.123362 | 1.275892 |
| `BETA-F3-rigidhalf` | -0.008652 | -0.019053 | 1.000196 | 1.006294 | 0.123362 | 1.275892 |

## Table 4 — Per-cell accounting flows, the eighteen sectoral cells

`S` is `diagnostics.canary_s`; `M = -(gdp_m_final + gdp_m_int)`; `I + X =`
`diagnostics.canary_ixm + M`; `tax` is `gdp_g`; `T = -gdp_t_int`; `F`, `B_gov` and
`Net ext. pos.` are the `external_transfer`, `programme_financing` and
`external_position` metrics (`Net ext. pos. = F + B_gov`). `Resource side =`
`S + T + M - (I+X)` is the imbalance the model's own flows imply; `Gap =`
`Resource side - Net ext. pos.` is the identity gap `canary_diff`.

| Cell | S | tax | I + X | M | T | F | B_gov | Net ext. pos. | Resource side | Gap |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `BETA-F1-etas025` | 0.094541 | 0.213727 | 0.588550 | 0.465915 | 0.025784 | -0.00230972 | +0.00000000 | -0.00230972 | -0.00230972 | -8.51e-13 |
| `BETA-F2-etas025` | 0.092632 | 0.227904 | 0.588840 | 0.465997 | 0.025787 | -0.00442385 | +0.00000000 | -0.00442385 | -0.00442385 | -4.49e-13 |
| `BETA-F3-etas025` | 0.092632 | 0.227904 | 0.588840 | 0.465997 | 0.025787 | -0.01863008 | +0.01420623 | -0.00442385 | -0.00442385 | -5.09e-13 |
| `BETA-F1-etas05` | 0.094619 | 0.213850 | 0.587313 | 0.466054 | 0.025767 | -0.00087292 | +0.00000000 | -0.00087292 | -0.00087292 | -4.49e-13 |
| `BETA-F2-etas05` | 0.092743 | 0.227780 | 0.587468 | 0.466128 | 0.025768 | -0.00282844 | +0.00000000 | -0.00282844 | -0.00282844 | -8.78e-13 |
| `BETA-F3-etas05` | 0.092743 | 0.227780 | 0.587468 | 0.466128 | 0.025768 | -0.01677266 | +0.01394422 | -0.00282844 | -0.00282844 | -9.92e-13 |
| `BETA-F1-etas1` | 0.094689 | 0.213951 | 0.586212 | 0.466184 | 0.025751 | +0.00041235 | +0.00000000 | +0.00041235 | +0.00041235 | +2.80e-13 |
| `BETA-F2-etas1` | 0.092841 | 0.227656 | 0.586282 | 0.466245 | 0.025751 | -0.00144371 | +0.00000000 | -0.00144371 | -0.00144371 | +9.23e-13 |
| `BETA-F3-etas1` | 0.092841 | 0.227656 | 0.586282 | 0.466245 | 0.025751 | -0.01515439 | +0.01371069 | -0.00144371 | -0.00144371 | -1.52e-12 |
| `BETA-F1-etas2` | 0.094742 | 0.214019 | 0.585410 | 0.466280 | 0.025740 | +0.00135125 | +0.00000000 | +0.00135125 | +0.00135125 | +2.48e-13 |
| `BETA-F2-etas2` | 0.092913 | 0.227557 | 0.585438 | 0.466329 | 0.025739 | -0.00045731 | +0.00000000 | -0.00045731 | -0.00045731 | +5.39e-13 |
| `BETA-F3-etas2` | 0.092913 | 0.227557 | 0.585438 | 0.466329 | 0.025739 | -0.01399821 | +0.01354089 | -0.00045731 | -0.00045731 | +5.34e-13 |
| `BETA-F1-rigidprog` | 0.093750 | 0.213676 | 0.590320 | 0.464062 | 0.025718 | -0.00678953 | +0.00000000 | -0.00678953 | -0.00678953 | -9.14e-12 |
| `BETA-F2-rigidprog` | 0.091687 | 0.228378 | 0.590983 | 0.463954 | 0.025713 | -0.00962905 | +0.00000000 | -0.00962905 | -0.00962905 | -1.11e-13 |
| `BETA-F3-rigidprog` | 0.091687 | 0.228378 | 0.590983 | 0.463954 | 0.025713 | -0.02438239 | +0.01475334 | -0.00962905 | -0.00962905 | -1.11e-13 |
| `BETA-F1-rigidhalf` | 0.094397 | 0.213467 | 0.591769 | 0.466578 | 0.025820 | -0.00497369 | +0.00000000 | -0.00497369 | -0.00497369 | +5.58e-13 |
| `BETA-F2-rigidhalf` | 0.092428 | 0.228018 | 0.592484 | 0.466765 | 0.025828 | -0.00746270 | +0.00000000 | -0.00746270 | -0.00746270 | -8.13e-13 |
| `BETA-F3-rigidhalf` | 0.092428 | 0.228018 | 0.592484 | 0.466765 | 0.025828 | -0.02207582 | +0.01461312 | -0.00746270 | -0.00746270 | -1.04e-12 |

## Notes

- **The account closes in every sectoral cell.** `Gap` is at the solver-residual
  level (|Gap| <= 1.1e-11 across all thirty-three cells, and <= 2e-13 in the
  mobile matrix cells), so `Net ext. pos.` is a closed-account quantity in the
  sectoral rows too: `F` is a solved equilibrium object there, not a pin.
- **Prices move with the demand composition.** In every sectoral cell
  `max abs(p-1)` differs between F1 and F2/F3, and the deflator is different from
  one; all twelve `eta = 1` matrix cells read exactly `0.000000` on both. This is
  the demand-sensitivity signature no single-wage closure produces.
- **The elasticity splits the adjustment.** As `eta_s` rises 0.25 -> 2 the price
  response falls (0.1529 -> 0.0372 at F2) while employment rises
  (1.000864 -> 1.001914) and the welfare cost falls (-1.69 % -> -1.39 %).
- **Incidence dominates the level.** Making the programme sectors rigid
  (`eta_s,i = 0` there) more than doubles the price response (0.2720 against
  0.1057) and reverses the employment effect (-0.34 % against +0.13 %); rigidity
  on the largest half by baseline employment gives 0.1234 with a near-neutral
  +0.02 %.
- **Wage dispersion is large.** `max_i w_i / min_i w_i` reads 1.26 at
  `eta_s = 0.5` and 1.70 in the rigid-programme variant, against 1.72 at the
  `eta = 0` endpoint, for a programme worth 1.33 % of GDP.
- **F1 is the welfare exception.** Its consumption index is positive
  (+0.42 % at `eta_s = 0.5`), because its programme is financed by a preference
  tilt rather than a tax.

## Status

All thirty-three cells of `matrix_5x3_v9` are `executed` and pass every gate
(residual, budget, sectoral labour-market gap; ADR-0019 acceptance). The
generation is the citable one for the ADR-0022 sectoral family; the fifteen
matrix rows may be cited from either v6 (the flows above) or v9 (bit-identical).

## Open items

- The **level** of `eta_s` is unidentified under demand-only shocks: the ladder is
  a sensitivity band, not an estimate. Closing it needs the supply-shock arm.
- The **grouping rule** (which sectors are rigid) is a modelling choice and flips
  the sign of the employment effect; report it as a band unless a principled rule
  is found.
- The **welfare metric** (Tornqvist on household consumption) is dispersion-blind.
- A **finer ladder** (`eta_s` in 0.1, 5) and a robustness probe on the sign under a
  different shock vector remain to be run.
