# Flow tables — the 5×3 matrix on the A-bill calibration (ADR-0015 generation)

Generated from `runs/matrix_5x3-v3-*/manifest.toml` and a re-solve of the same 15
design cells (commit `89cfde5`). All flows are in model units: GDP at basic
prices = 1, so a flow of 0.1 is 10 % of GDP (GDP_P = 3 027 818 EUR m). Paper text
cites the `matrix_5x3-v3-*` run ids; never retype these numbers. Supersedes the
`matrix_5x3-v2` version of this table (four cells moved at the 1e-5 relative
level under the ADR-0015 monotone polish; the other eleven are unchanged to
machine precision).

## Table 1 — Baseline accounts (full-71 A-bill calibration)

| Account | Value | Source |
| --- | ---: | --- |
| GDP, production = income | 1.000000 | `Σ gva` = `Σ λ·fs` = 1 |
| GDP, expenditure | 0.946136 | domestic final demand at basic prices |
| Production-vs-expenditure residual | 5.387 % (163 094 EUR m) | raw table's own gap, open item |
| Household income (after tax) | 0.785899 | `1 − τ0` |
| Household consumption (purchaser prices) | 0.691675 | `c0`, the table's own household column |
| Saving rate `s` | 0.119892 | identity-implied (ADR-0012) |
| Household saving `S` | 0.094223 | `s · (1 − τ0)` |
| Government spending = lump-sum tax `τ0` | 0.214101 | `Σ gG` |
| Investment `I` | 0.162367 | equipment + construction + inventories |
| Exports `X` | 0.421945 | no import margin (domestic sales abroad) |
| Final-demand imports `M_final` | 0.242941 | margin content of C + G + I |
| Intermediate imports `M_int` (row 74) | 0.221403 | booked leak (ADR-0012) |
| Intermediate product taxes `T_int` (row 75) | 0.025745 | booked leak (ADR-0013) |
| Total imports `M` at baseline | 0.464344 | `M_final + M_int` |
| Domestic intermediate bill `ΣA_bill` (row 73) | 0.862762 | charged to domestic demand |
| Row identity `ΣA + ΣM_int + ΣT_int = Σλ − 1` | 0.247147 | holds to 2e-16 |
| Baseline external identity `S − (I+X−M) + T` | −5.6e-17 | machine zero (ADR-0013) |
| Programme `G0 = Σ g` (2024 impulses) | 0.013310 (40 300 EUR m) | `[programme]` in the design |

## Table 2 — Per-cell accounting flows

`Real GDP rel.` is the Törnqvist index against the baseline (baseline = 1); `L` is
total employment (`L̄` = 1); `Net ext. position` is the model's own
identity-consistent external balance, `S − (I+X−M) + T`, i.e. the omitted-market
canary. Negative values are net inflows (external deficits). Every cell's
residual is at or below 1e-12 except where noted in its manifest (ADR-0015
monotone polish).

| Cell | Real GDP rel. | L | S | τ (tax) | I + X | M | T | Net ext. position |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `BF-F1` | +0.000433 | 1.000000 | 0.094223 | 0.214101 | 0.584313 | 0.464884 | 0.025657 | +0.000452 |
| `BF-F2` | -0.016936 | 1.000000 | 0.092628 | 0.227411 | 0.584313 | 0.465459 | 0.025680 | -0.000546 |
| `BF-F3` | -0.000000 | 1.000000 | 0.094223 | 0.214101 | 0.584313 | 0.456280 | 0.025862 | -0.007947 |
| `ALPHA-F1` | +0.000433 | 1.000000 | 0.094223 | 0.214101 | 0.584313 | 0.464895 | 0.025659 | +0.000465 |
| `ALPHA-F2` | -0.016936 | 1.000000 | 0.092628 | 0.227411 | 0.584313 | 0.465411 | 0.025673 | -0.000601 |
| `ALPHA-F3` | +0.000000 | 1.000000 | 0.094223 | 0.214101 | 0.584313 | 0.455811 | 0.025783 | -0.008495 |
| `BETA-F1` | +0.000433 | 1.000000 | 0.094223 | 0.214101 | 0.584313 | 0.464895 | 0.025659 | +0.000465 |
| `BETA-F2` | -0.016936 | 1.000000 | 0.092628 | 0.227411 | 0.584313 | 0.465411 | 0.025673 | -0.000601 |
| `BETA-F3` | +0.000000 | 1.000000 | 0.094223 | 0.214101 | 0.584313 | 0.455811 | 0.025783 | -0.008495 |
| `GAMMA-F1` | -0.000806 | 0.999027 | 0.094107 | 0.214101 | 0.584313 | 0.464565 | 0.025641 | +0.000000 |
| `GAMMA-F2` | -0.015333 | 1.001260 | 0.092779 | 0.227411 | 0.584313 | 0.465838 | 0.025696 | -0.000000 |
| `GAMMA-F3` | +0.022644 | 1.017796 | 0.096357 | 0.214101 | 0.584313 | 0.461844 | 0.026112 | -0.000000 |
| `DELTA-F1` | -0.000806 | 0.999027 | 0.094107 | 0.214101 | 0.584313 | 0.464565 | 0.025641 | -0.000000 |
| `DELTA-F2` | -0.015333 | 1.001260 | 0.092779 | 0.227411 | 0.584313 | 0.465838 | 0.025696 | -0.000000 |
| `DELTA-F3` | +0.022644 | 1.017796 | 0.096357 | 0.214101 | 0.584313 | 0.461844 | 0.026112 | -0.000000 |

## Table 3 — The external account by financing closure

`F = dot(p, g) = 0.013310` is the programme's value under F3 (externally
financed); the programme's own import content is `0.002851`. The net external
position is the canary of Table 2.

| Financing | What the closure does | Net ext. position (mobile rows) | Net ext. position (fixed-wage rows) |
| --- | --- | ---: | ---: |
| F1 | pure composition shift, budget neutral | +0.000465 | +0.000008 / −0.000007 |
| F2 | programme financed by a lump-sum tax (`τ0` 0.2141 → 0.2274) | −0.000601 | −0.000002 |
| F3 | programme financed externally (inflow `F` = 0.013310) | −0.008495 | −2.5e-16 |

Reading F3: the inflow is `F` = 0.013310, of which the programme's own import
content is 0.002851; the *net* external position is −0.008495 for the mobile η = 1
rows, so the endogenous response offsets 0.004815 (36 %) of the inflow. The
fixed-wage rows have all N markets clearing, so the inflow is absorbed entirely
and the net position is zero to machine precision. The pre-registered signature
"external deficit = F" therefore describes the *gross inflow*, not the net
position; the paper should report both.

## Notes

- The three leaks are separated: `M` is the import content of final demand plus
  `M_int`; `T` is the product taxes on intermediate use (row 75); `S` is
  household saving. `S + M + T = I + X` holds to machine precision at every
  baseline; at a cell the difference is the net external position.
- `τ (tax)` is the `public_budget` diagnostic: government spending including the
  programme under F2 (`τ0 + G0`), and `τ0` otherwise.
- Employment is pinned at `L̄ = 1` in the mobile rows (ALPHA/BETA) and endogenous
  in the fixed-wage rows (GAMMA/DELTA), which is why the F3 column separates them.
