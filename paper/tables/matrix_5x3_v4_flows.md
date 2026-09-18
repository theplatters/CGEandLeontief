# Flow tables — the 5×3 matrix on the A-bill calibration (ADR-0018 generation)

Generated from `runs/matrix_5x3-v4-*/manifest.toml` (schema v2, commit
`f360d0c`), no re-solve and no hand-typed numbers. All flows are in model
units: GDP at basic prices = 1, so a flow of 0.1 is 10 % of GDP
(GDP_P = 3 027 818 EUR m). Paper text cites the `matrix_5x3-v4-*` run ids.

Supersedes the `matrix_5x3_v3` version of this table. The v3 table's
"Real GDP rel." column is the household-consumption (welfare) Törnqvist index;
under ADR-0018 that column is reported here as `Consumption rel.`, and real GDP
is the income-side index `gdp_rel = (w·ΣL)/(w·ΣL)_ref − 1` (the GDP deflator is
exactly one under demand-only shocks, so it drops out of this matrix). The
v1–v3 runs stay as history (ADR-0004).

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
| Baseline external identity `S - (I+X-M) + T` | -5.6e-17 | machine zero (ADR-0013) |
| Programme `G0 = Sigma g` (2024 impulses) | 0.013310 (40 300 EUR m) | `[programme]` in the design |

## Table 2 -- Per-cell accounting flows

`Real GDP rel.` is the ADR-0018 income-side index (`gdp_rel`),
`Consumption rel.` the household-consumption (welfare) index
(`consumption_rel`). All other flows are manifest levels in model units
(GDP at basic prices = 1 at the calibration baseline).

| Cell | Real GDP rel. | Consumption rel. | L | S | tax | I + X | M | T | Net ext. position |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `BF-F1` | +0.000000 | +0.000433 | 1.000000 | 0.094223 | 0.214101 | 0.584313 | 0.464884 | 0.025657 | +0.000452 |
| `BF-F2` | +0.000000 | -0.016936 | 1.000000 | 0.092628 | 0.227411 | 0.584313 | 0.465457 | 0.025680 | -0.000548 |
| `BF-F3` | +0.000000 | +0.000000 | 1.000000 | 0.094223 | 0.214101 | 0.584313 | 0.469590 | 0.025862 | -0.007947 |
| `ALPHA-F1` | +0.000000 | +0.000433 | 1.000000 | 0.094223 | 0.214101 | 0.584313 | 0.464895 | 0.025659 | +0.000465 |
| `ALPHA-F2` | +0.000000 | -0.016936 | 1.000000 | 0.092628 | 0.227411 | 0.584313 | 0.465411 | 0.025673 | -0.000601 |
| `ALPHA-F3` | +0.000000 | +0.000000 | 1.000000 | 0.094223 | 0.214101 | 0.584313 | 0.469120 | 0.025783 | -0.008495 |
| `BETA-F1` | +0.000000 | +0.000433 | 1.000000 | 0.094223 | 0.214101 | 0.584313 | 0.464895 | 0.025659 | +0.000465 |
| `BETA-F2` | +0.000000 | -0.016936 | 1.000000 | 0.092628 | 0.227411 | 0.584313 | 0.465411 | 0.025673 | -0.000601 |
| `BETA-F3` | +0.000000 | +0.000000 | 1.000000 | 0.094223 | 0.214101 | 0.584313 | 0.469120 | 0.025783 | -0.008495 |
| `GAMMA-F1` | -0.000973 | -0.000806 | 0.999027 | 0.094107 | 0.214101 | 0.584313 | 0.464565 | 0.025641 | -0.000000 |
| `GAMMA-F2` | +0.001260 | -0.015333 | 1.001260 | 0.092779 | 0.227411 | 0.584313 | 0.465838 | 0.025696 | -0.000000 |
| `GAMMA-F3` | +0.017796 | +0.022644 | 1.017796 | 0.096357 | 0.214101 | 0.584313 | 0.475154 | 0.026112 | -0.000000 |
| `DELTA-F1` | -0.000973 | -0.000806 | 0.999027 | 0.094107 | 0.214101 | 0.584313 | 0.464565 | 0.025641 | -0.000000 |
| `DELTA-F2` | +0.001260 | -0.015333 | 1.001260 | 0.092779 | 0.227411 | 0.584313 | 0.465838 | 0.025696 | -0.000000 |
| `DELTA-F3` | +0.017796 | +0.022644 | 1.017796 | 0.096357 | 0.214101 | 0.584313 | 0.475154 | 0.026112 | -0.000000 |

## Table 2b -- Expenditure components (manifest diagnostics)

| Cell | C_gross | G+prog | I | X | M_final | M_int | T_int | wedge | deflator |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `BF-F1` | 0.691675 | 0.214101 | 0.162367 | 0.421945 | -0.242482 | -0.222402 | -0.025657 | -4.52e-04 | 1.000000 |
| `BF-F2` | 0.679961 | 0.227411 | 0.162367 | 0.421945 | -0.242871 | -0.222586 | -0.025680 | +5.48e-04 | 1.000000 |
| `BF-F3` | 0.691675 | 0.227411 | 0.162367 | 0.421945 | -0.245792 | -0.223798 | -0.025862 | +7.95e-03 | 1.000000 |
| `ALPHA-F1` | 0.691675 | 0.214101 | 0.162367 | 0.421945 | -0.242482 | -0.222413 | -0.025659 | -4.65e-04 | 1.000000 |
| `ALPHA-F2` | 0.679961 | 0.227411 | 0.162367 | 0.421945 | -0.242871 | -0.222540 | -0.025673 | +6.01e-04 | 1.000000 |
| `ALPHA-F3` | 0.691675 | 0.227411 | 0.162367 | 0.421945 | -0.245792 | -0.223328 | -0.025783 | +8.50e-03 | 1.000000 |
| `BETA-F1` | 0.691675 | 0.214101 | 0.162367 | 0.421945 | -0.242482 | -0.222413 | -0.025659 | -4.65e-04 | 1.000000 |
| `BETA-F2` | 0.679961 | 0.227411 | 0.162367 | 0.421945 | -0.242871 | -0.222540 | -0.025673 | +6.01e-04 | 1.000000 |
| `BETA-F3` | 0.691675 | 0.227411 | 0.162367 | 0.421945 | -0.245792 | -0.223328 | -0.025783 | +8.50e-03 | 1.000000 |
| `GAMMA-F1` | 0.690819 | 0.214101 | 0.162367 | 0.421945 | -0.242269 | -0.222296 | -0.025641 | -2.22e-16 | 1.000000 |
| `GAMMA-F2` | 0.681070 | 0.227411 | 0.162367 | 0.421945 | -0.243147 | -0.222691 | -0.025696 | +0.00e+00 | 1.000000 |
| `GAMMA-F3` | 0.707337 | 0.227411 | 0.162367 | 0.421945 | -0.249699 | -0.225455 | -0.026112 | +0.00e+00 | 1.000000 |
| `DELTA-F1` | 0.690819 | 0.214101 | 0.162367 | 0.421945 | -0.242269 | -0.222296 | -0.025641 | +5.04e-12 | 1.000000 |
| `DELTA-F2` | 0.681070 | 0.227411 | 0.162367 | 0.421945 | -0.243147 | -0.222691 | -0.025696 | +1.54e-12 | 1.000000 |
| `DELTA-F3` | 0.707337 | 0.227411 | 0.162367 | 0.421945 | -0.249699 | -0.225455 | -0.026112 | +5.35e-12 | 1.000000 |

## Table 3 -- The external account by financing closure

| Financing | Net ext. position (mobile rows) | Net ext. position (fixed-wage rows) |
| --- | ---: | ---: |
| F1 | +0.000465 / +0.000465 | -0.000000 / -0.000000 |
| F2 | -0.000601 / -0.000601 | -0.000000 / -0.000000 |
| F3 | -0.008495 / -0.008495 | -0.000000 / -0.000000 |

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
- `GAMMA-F1`: executed (pass, resid 8.9e-16)
- `GAMMA-F2`: executed (pass, resid 8.9e-16)
- `GAMMA-F3`: executed (pass, resid 8.9e-16)
- `DELTA-F1`: executed (pass, resid 1.3e-13)
- `DELTA-F2`: executed (pass, resid 4.5e-14)
- `DELTA-F3`: executed (pass, resid 1.4e-13)

## Notes

- `Real GDP rel.` (`gdp_rel`, ADR-0018) is the income-side index
  `(w·SigmaL)/(w·SigmaL)_ref - 1`: 0.000 % in every full-employment mobile row
  (BF/ALPHA/BETA), -0.097 % / +0.126 % / +1.780 % in the fixed-wage rows
  (GAMMA/DELTA) under F1/F2/F3. `Consumption rel.` is the household-consumption
  (welfare) Törnqvist index: +0.043 % under F1, -1.694 % (mobile) / -1.533 %
  (fixed) under F2, 0.000 % (mobile) / +2.264 % (fixed) under F3.
- The F2 column is the sharpest separation: the tax-financed programme leaves
  real GDP flat at full employment (the pre-registered "aggregate ~ 0"
  signature holds for GDP) while the household welfare index falls by the tax.
  The v3 table reported the welfare fall as "real GDP"; ADR-0018 separates the
  two and v4 is the citable generation for either.
- The **wedge** in Table 2b is `SigmaV - w·SigmaL` (the manifest `gdp_wedge`),
  i.e. minus the external position `p·rho` carried by the omitted N-th market.
  It is zero to machine precision at the baseline and at the fixed-wage cells
  (all N markets clear) and equals the net external position in the mobile
  eta = 1 rows (+0.06 % of GDP for F2, +0.85 % for F3). The expenditure-side
  Divisia index counts it; the income-side real GDP does not.
- `M = M_final + M_int` (the `gdp_m_final`/`gdp_m_int` manifest diagnostics,
  reported negative there); `T` is `gdp_t_int`. `S + M + T = I + X` holds to
  machine precision at every cell except where the net external position is
  nonzero (F1/F2 mobile and F3 mobile), which is exactly the `Net ext.
  position` column.
- `tax` is the `public_budget` diagnostic: government spending including the
  programme under F2 (`tau0 + G0`), and `tau0` otherwise.
- All 15 cells pass their gates (residual <= 1.4e-13, ADR-0015 polish);
  v4 reproduces the v3 equilibrium values to <= 1.2e-11 (the two generations
  differ only in the measurement layer).
