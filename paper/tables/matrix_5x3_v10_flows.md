# Flow tables — the 5x3 matrix and the sectoral family (ADR-0022 corrected generation)

Generated from `runs/matrix_5x3-v10-*/manifest.toml` (schema v2, design commit `979a58e`, run commit `c2c5b7d`),
no re-solve and no hand-typed numbers: every figure below is a manifest field or a
sum of manifest fields, with the mapping stated in the table notes. All flows are in
model units: GDP at basic prices = 1, so a flow of 0.1 is 10 % of GDP
(GDP_P = 3 027 818 EUR m). Paper text cites the `matrix_5x3-v10-*` run ids.

This generation re-mints all 33 v9 cells (15 matrix + 18 sectoral) on the C1-fixed kernel,
which uses the sectoral wage vector in `gdp_components` on the ADR-0022 cells
(η = 1 with `eta_s_vec`). It supersedes `matrix_5x3_v9` for the 18 sectoral cells'
`Real GDP rel.`/`deflator` columns (Table 3) and their `M`/`I + X` columns (Table 4);
the v9 runs stay as issued (ADR-0004). Its closure is ADR-0019 + ADR-0020 option C as in v6,
plus ADR-0022: a BETA cell carrying an elasticity vector solves the N sectoral
labour markets as the 3N+1 system. With `eta_s,i = 0` for every sector that system
*is* the v6 `BF` endpoint (proved in `paper/equivalence.tex` v3, Proposition 3), so
no separate rigid sectoral cell is reported: the corner is read from the `BF` cells.

## Table 1 — Baseline accounts

Unchanged from `paper/tables/matrix_5x3_v6_flows.md`, Table 1: the calibration
(full-71 A-bill, ADR-0012/ADR-0013) and the reference continuation are identical in
the two generations, and the fifteen matrix cells of this generation reproduce the
v9 cells to solver precision (Table 2). The baseline block is therefore not duplicated here.

## Table 2 — The fifteen matrix cells: solver-precision agreement with v9

The `Delta` column is `max |v10 - v9|` over that cell's `[metrics]` and `[diagnostics]`
fields. v9's Table 2 records a zero v9-vs-v6 delta in every cell, so v6 and v9 are
bit-identical to each other; v10 ran in this working copy, so the fifteen cells reproduce
v9 to solver precision — max 9.34e-09 (at `BF-F1`, the stiff corner's `diagnostics.wage_max` diagnostic),
1.63e-11 in the worst ALPHA/BETA/GAMMA/DELTA row (`GAMMA-F3`, `metrics.consumption`),
4.55e-15 where the equilibrium is the exact baseline — and no printed figure is affected
at the table's precision (15/15 rows identical to v9 to all printed decimals).

| Cell | Real GDP rel. | Consumption rel. | L | deflator | `max abs(p-1)` | Delta |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `BF-F1` | -0.000110 | +0.000963 | 1.000000 | 1.006326 | 0.221442 | +9.34e-09 |
| `BF-F2` | -0.000162 | -0.019800 | 1.000000 | 1.007181 | 0.279009 | +1.22e-13 |
| `BF-F3` | -0.000162 | -0.019800 | 1.000000 | 1.007181 | 0.279009 | +8.04e-14 |
| `ALPHA-F1` | +0.000000 | +0.001432 | 1.000000 | 1.000000 | 0.000000 | +4.22e-15 |
| `ALPHA-F2` | +0.000000 | -0.018226 | 1.000000 | 1.000000 | 0.000000 | +4.44e-15 |
| `ALPHA-F3` | +0.000000 | -0.018226 | 1.000000 | 1.000000 | 0.000000 | +4.55e-15 |
| `BETA-F1` | +0.000000 | +0.001432 | 1.000000 | 1.000000 | 0.000000 | +4.00e-15 |
| `BETA-F2` | +0.000000 | -0.018226 | 1.000000 | 1.000000 | 0.000000 | +4.44e-15 |
| `BETA-F3` | +0.000000 | -0.018226 | 1.000000 | 1.000000 | 0.000000 | +4.55e-15 |
| `GAMMA-F1` | -0.000973 | -0.000806 | 0.999027 | 1.000000 | 0.000000 | +9.11e-12 |
| `GAMMA-F2` | +0.001260 | -0.015333 | 1.001260 | 1.000000 | 0.000000 | +8.62e-12 |
| `GAMMA-F3` | +0.017796 | +0.022644 | 1.017796 | 1.000000 | 0.000000 | +1.63e-11 |
| `DELTA-F1` | -0.000973 | -0.000806 | 0.999027 | 1.000000 | 0.000000 | +1.00e-11 |
| `DELTA-F2` | +0.001260 | -0.015333 | 1.001260 | 1.000000 | 0.000000 | +9.28e-13 |
| `DELTA-F3` | +0.017796 | +0.022644 | 1.017796 | 1.000000 | 0.000000 | +5.33e-15 |

## Table 3 — Headline metrics, the eighteen sectoral cells (corrected)

`wage max/min` is `diagnostics.wage_max / wage_min`; the reported `wage` metric is
the wage-bill-weighted aggregate (`sum_i w_i L_i / sum_i L_i`). F2 and F3 agree in
every sectoral cell to all printed decimals (financing neutrality, ADR-0019; 6/6 variant pairs identical).

| Cell | Real GDP rel. | Consumption rel. | L | deflator | `max abs(p-1)` | wage max/min |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `BETA-F1-etas025` | +0.000811 | +0.003252 | 1.000809 | 1.003771 | 0.133610 | 1.331995 |
| `BETA-F2-etas025` | +0.000864 | -0.016888 | 1.000864 | 1.004088 | 0.152897 | 1.381559 |
| `BETA-F3-etas025` | +0.000864 | -0.016888 | 1.000864 | 1.004088 | 0.152897 | 1.381559 |
| `BETA-F1-etas05` | +0.001217 | +0.004227 | 1.001191 | 1.002699 | 0.095839 | 1.234983 |
| `BETA-F2-etas05` | +0.001287 | -0.015712 | 1.001255 | 1.002868 | 0.105666 | 1.259605 |
| `BETA-F3-etas05` | +0.001287 | -0.015712 | 1.001255 | 1.002868 | 0.105666 | 1.259605 |
| `BETA-F1-etas1` | +0.001597 | +0.005119 | 1.001564 | 1.001725 | 0.061283 | 1.148364 |
| `BETA-F2-etas1` | +0.001667 | -0.014669 | 1.001627 | 1.001800 | 0.065422 | 1.158491 |
| `BETA-F3-etas1` | +0.001667 | -0.014669 | 1.001627 | 1.001800 | 0.065422 | 1.158491 |
| `BETA-F1-etas2` | +0.001886 | +0.005784 | 1.001859 | 1.001003 | 0.035629 | 1.085418 |
| `BETA-F2-etas2` | +0.001945 | -0.013912 | 1.001914 | 1.001033 | 0.037171 | 1.089127 |
| `BETA-F3-etas2` | +0.001945 | -0.013912 | 1.001914 | 1.001033 | 0.037171 | 1.089127 |
| `BETA-F1-rigidprog` | -0.003012 | -0.005416 | 0.996999 | 1.005448 | 0.215061 | 1.544801 |
| `BETA-F2-rigidprog` | -0.003414 | -0.026914 | 0.996621 | 1.006190 | 0.272033 | 1.697924 |
| `BETA-F3-rigidprog` | -0.003414 | -0.026914 | 0.996621 | 1.006190 | 0.272033 | 1.697924 |
| `BETA-F1-rigidhalf` | +0.000144 | +0.001496 | 1.000188 | 1.005639 | 0.109583 | 1.247681 |
| `BETA-F2-rigidhalf` | +0.000124 | -0.019053 | 1.000196 | 1.006281 | 0.123362 | 1.275892 |
| `BETA-F3-rigidhalf` | +0.000124 | -0.019053 | 1.000196 | 1.006281 | 0.123362 | 1.275892 |

### Correction record (C1)

`gdp_components` collapsed the sectoral wage vector to its first entry on the eighteen
ADR-0022 cells, so the v9 `gdp_rel` below is wrong in every sectoral row (the expenditure
side and the deflator inherit the error). The corrected income-side `gdp_rel` flips sign in
15/18 cells; the v9 `gdp_wedge` (up to 3.8e-03) collapses to 1.5e-12 or better in every
sectoral cell (`gdp_wedge = -canary_diff` is now hard-asserted; max violation 4.4e-16). All other
metrics/diagnostics fields agree with v9 to cross-copy float noise (max 7.856e-09 at
`BETA-F1-rigidprog`, `diagnostics.wage_max`).

| Cell | v9 Real GDP rel. (wrong) | v10 Real GDP rel. | v9 deflator | v10 deflator |
| --- | ---: | ---: | ---: | ---: |
| `BETA-F1-etas025` | -0.007086 | +0.000811 | 1.003776 | 1.003771 |
| `BETA-F2-etas025` | -0.007671 | +0.000864 | 1.004096 | 1.004088 |
| `BETA-F3-etas025` | -0.007671 | +0.000864 | 1.004096 | 1.004088 |
| `BETA-F1-etas05` | -0.004414 | +0.001217 | 1.002701 | 1.002699 |
| `BETA-F2-etas05` | -0.004690 | +0.001287 | 1.002872 | 1.002868 |
| `BETA-F3-etas05` | -0.004690 | +0.001287 | 1.002872 | 1.002868 |
| `BETA-F1-etas1` | -0.001981 | +0.001597 | 1.001726 | 1.001725 |
| `BETA-F2-etas1` | -0.002069 | +0.001667 | 1.001802 | 1.001800 |
| `BETA-F3-etas1` | -0.002069 | +0.001667 | 1.001802 | 1.001800 |
| `BETA-F1-etas2` | -0.000184 | +0.001886 | 1.001003 | 1.001003 |
| `BETA-F2-etas2` | -0.000189 | +0.001945 | 1.001033 | 1.001033 |
| `BETA-F3-etas2` | -0.000189 | +0.001945 | 1.001033 | 1.001033 |
| `BETA-F1-rigidprog` | -0.012839 | -0.003012 | 1.005456 | 1.005448 |
| `BETA-F2-rigidprog` | -0.014525 | -0.003414 | 1.006206 | 1.006190 |
| `BETA-F3-rigidprog` | -0.014525 | -0.003414 | 1.006206 | 1.006190 |
| `BETA-F1-rigidhalf` | -0.007748 | +0.000144 | 1.005646 | 1.005639 |
| `BETA-F2-rigidhalf` | -0.008652 | +0.000124 | 1.006294 | 1.006281 |
| `BETA-F3-rigidhalf` | -0.008652 | +0.000124 | 1.006294 | 1.006281 |

## Table 4 — Per-cell accounting flows, the eighteen sectoral cells

`S` is `diagnostics.canary_s`; `M = -(gdp_m_final + gdp_m_int)`; `I + X =`
`diagnostics.canary_ixm + M`; `tax` is `gdp_g`; `T = -gdp_t_int`; `F`, `B_gov` and
`Net ext. pos.` are the `external_transfer`, `programme_financing` and
`external_position` metrics (`Net ext. pos. = F + B_gov`). `Resource side =`
`S + T + M - (I+X)` is the imbalance the model's own flows imply; `Gap =`
`Resource side - Net ext. pos.` is the identity gap `canary_diff`.

| Cell | S | tax | I + X | M | T | F | B_gov | Net ext. pos. | Resource side | Gap |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `BETA-F1-etas025` | 0.094541 | 0.213727 | 0.590284 | 0.467649 | 0.025784 | -0.00230972 | +0.00000000 | -0.00230972 | -0.00230972 | -8.93e-13 |
| `BETA-F2-etas025` | 0.092632 | 0.227904 | 0.590719 | 0.467877 | 0.025787 | -0.00442385 | +0.00000000 | -0.00442385 | -0.00442385 | -4.49e-13 |
| `BETA-F3-etas025` | 0.092632 | 0.227904 | 0.590719 | 0.467877 | 0.025787 | -0.01863008 | +0.01420623 | -0.00442385 | -0.00442385 | -4.77e-13 |
| `BETA-F1-etas05` | 0.094619 | 0.213850 | 0.588549 | 0.467290 | 0.025767 | -0.00087292 | +0.00000000 | -0.00087292 | -0.00087292 | -4.85e-13 |
| `BETA-F2-etas05` | 0.092743 | 0.227780 | 0.588783 | 0.467443 | 0.025768 | -0.00282844 | +0.00000000 | -0.00282844 | -0.00282844 | -8.66e-13 |
| `BETA-F3-etas05` | 0.092743 | 0.227780 | 0.588783 | 0.467443 | 0.025768 | -0.01677266 | +0.01394422 | -0.00282844 | -0.00282844 | -9.46e-13 |
| `BETA-F1-etas1` | 0.094689 | 0.213951 | 0.586997 | 0.466968 | 0.025751 | +0.00041235 | +0.00000000 | +0.00041235 | +0.00041235 | 4.11e-13 |
| `BETA-F2-etas1` | 0.092841 | 0.227656 | 0.587103 | 0.467067 | 0.025751 | -0.00144371 | +0.00000000 | -0.00144371 | -0.00144371 | 1.05e-12 |
| `BETA-F3-etas1` | 0.092841 | 0.227656 | 0.587103 | 0.467067 | 0.025751 | -0.01515439 | +0.01371069 | -0.00144371 | -0.00144371 | -1.49e-12 |
| `BETA-F1-etas2` | 0.094742 | 0.214019 | 0.585864 | 0.466734 | 0.025740 | +0.00135125 | +0.00000000 | +0.00135125 | +0.00135125 | 2.65e-13 |
| `BETA-F2-etas2` | 0.092913 | 0.227557 | 0.585907 | 0.466798 | 0.025739 | -0.00045731 | +0.00000000 | -0.00045731 | -0.00045731 | 5.49e-13 |
| `BETA-F3-etas2` | 0.092913 | 0.227557 | 0.585907 | 0.466798 | 0.025739 | -0.01399821 | +0.01354089 | -0.00045731 | -0.00045731 | 4.66e-13 |
| `BETA-F1-rigidprog` | 0.093750 | 0.213676 | 0.592481 | 0.466223 | 0.025718 | -0.00678953 | +0.00000000 | -0.00678953 | -0.00678953 | -4.30e-13 |
| `BETA-F2-rigidprog` | 0.091687 | 0.228378 | 0.593434 | 0.466404 | 0.025713 | -0.00962905 | +0.00000000 | -0.00962905 | -0.00962905 | -1.31e-13 |
| `BETA-F3-rigidprog` | 0.091687 | 0.228378 | 0.593434 | 0.466404 | 0.025713 | -0.02438239 | +0.01475334 | -0.00962905 | -0.00962905 | -1.24e-13 |
| `BETA-F1-rigidhalf` | 0.094397 | 0.213467 | 0.593504 | 0.468314 | 0.025820 | -0.00497369 | +0.00000000 | -0.00497369 | -0.00497369 | 5.85e-13 |
| `BETA-F2-rigidhalf` | 0.092428 | 0.228018 | 0.594420 | 0.468701 | 0.025828 | -0.00746270 | +0.00000000 | -0.00746270 | -0.00746270 | -8.08e-13 |
| `BETA-F3-rigidhalf` | 0.092428 | 0.228018 | 0.594420 | 0.468701 | 0.025828 | -0.02207582 | +0.01461312 | -0.00746270 | -0.00746270 | -1.04e-12 |

## Notes

- **The account closes in every sectoral cell.** `Gap` is at the solver-residual
  level (|Gap| <= 7.1e-12 across all 33 cells — worst at `GAMMA-F3` — and <= 5.41e-13 in the
  non-fixed-wage matrix cells), so `Net ext. pos.` is a closed-account quantity in the
  sectoral rows too: `F` is a solved equilibrium object there, not a pin.
- **Prices move with the demand composition.** In every sectoral cell
  `max abs(p-1)` differs between F1 and F2/F3 (6/6 variants), and the deflator differs from
  one in all 18 sectoral rows; all 12 `eta = 1` matrix cells read exactly `0.000000` on both. This is
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

All 33 cells of `matrix_5x3_v10` are `executed` (33/33) and pass every gate (33/33)
(residual, budget, sectoral labour-market gap; ADR-0019 acceptance). The
generation is the citable one for the ADR-0022 sectoral family and for the matrix; the fifteen
matrix rows may be cited from either v6/v9 (solver-precision agreement, Table 2) or v10.

## Open items

- The **level** of `eta_s` is unidentified under demand-only shocks: the ladder is
  a sensitivity band, not an estimate. Closing it needs the supply-shock arm.
- The **grouping rule** (which sectors are rigid) is a modelling choice and flips
  the sign of the employment effect; report it as a band unless a principled rule
  is found.
- The **welfare metric** (Tornqvist on household consumption) is dispersion-blind.
- A **finer ladder** (`eta_s` in 0.1, 5) and a robustness probe on the sign under a
  different shock vector remain to be run.
