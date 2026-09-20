---
title: "Concise summary: the 5x3 matrix, the sectoral labour door, and the GAMMA menu"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-09-19"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [summary, matrix, labour-closure, prices, gamma, handoff]
last-updated: September 2026
---

**Version 2** \textcolor{revisionV1}{(September 2026)}
**Version 1** (September 2026)

One page of results, no derivation: the executed 5x3 matrix, the sectoral
labour-market family that gives the flexible-wage rows a demand channel, and the
seven fixed-wage (GAMMA) variants measured side by side. Everything below is a
manifest field or a probe measurement; the sources are listed at the end, and
the derivations live in `paper/equivalence.tex` (v3), `paper/framing_gamma.tex`
and `docs/WORKPLAN_SENSITIVE_PRICES.md` (v7).

# The 5x3 matrix, current results

Generation `matrix_5x3_v9` (33 cells, all gates pass). The fifteen matrix cells
reproduce the `matrix_5x3_v6` values **bit-for-bit** (`max|Delta| = 0` on every
metric and diagnostic), so the equivalence results of `paper/equivalence.tex`
are invariant to the extension. Deviations from the no-programme baseline; the
programme is 1.33 % of GDP.

| Cell | Real GDP rel. | Consumption rel. | Employment | Deflator | `max abs(p-1)` |
| --- | ---: | ---: | ---: | ---: | ---: |
| `BF-F1` | -0.000110 | +0.000963 | 1.000000 | 1.006326 | 0.221442 |
| `BF-F2` | -0.000162 | -0.019800 | 1.000000 | 1.007181 | 0.279009 |
| `BF-F3` | -0.000162 | -0.019800 | 1.000000 | 1.007181 | 0.279009 |
| `ALPHA-F1` | +0.000000 | +0.001432 | 1.000000 | 1.000000 | 0.000000 |
| `ALPHA-F2` | +0.000000 | -0.018226 | 1.000000 | 1.000000 | 0.000000 |
| `ALPHA-F3` | +0.000000 | -0.018226 | 1.000000 | 1.000000 | 0.000000 |
| `BETA-F1` | +0.000000 | +0.001432 | 1.000000 | 1.000000 | 0.000000 |
| `BETA-F2` | +0.000000 | -0.018226 | 1.000000 | 1.000000 | 0.000000 |
| `BETA-F3` | +0.000000 | -0.018226 | 1.000000 | 1.000000 | 0.000000 |
| `GAMMA-F1` | -0.000973 | -0.000806 | 0.999027 | 1.000000 | 0.000000 |
| `GAMMA-F2` | +0.001260 | -0.015333 | 1.001260 | 1.000000 | 0.000000 |
| `GAMMA-F3` | +0.017796 | +0.022644 | 1.017796 | 1.000000 | 0.000000 |
| `DELTA-F1` | -0.000973 | -0.000806 | 0.999027 | 1.000000 | 0.000000 |
| `DELTA-F2` | +0.001260 | -0.015333 | 1.001260 | 1.000000 | 0.000000 |
| `DELTA-F3` | +0.017796 | +0.022644 | 1.017796 | 1.000000 | 0.000000 |

Three readings.

- **Only the BF row moves prices.** `max abs(p-1)` and the deflator are exactly
  zero/one in the twelve `eta = 1` cells of ALPHA, BETA, GAMMA and DELTA, and
  non-zero in BF (0.221-0.279; deflator 1.0063-1.0072). With one factor,
  constant returns and a demand-free price block, demand can reach prices only
  through an endogenous wage - which is why the repair belongs in the wage
  block.
- **The two equivalences are exact.** ALPHA = BETA (the elastic-supply row
  reproduces the full-employment row) and GAMMA = DELTA (the Leontief corner
  reproduces the fixed-wage row), both to all printed decimals in every
  financing column. Both are now *proved* ex ante, not merely observed
  (Propositions 1 and 2 of `paper/equivalence.tex`).
- **The financing spread is the live dimension in the fixed-wage row.** F3 is
  the Keynesian case (+1.78 % employment, +2.26 % consumption), F2 the
  tax-financed loss (-1.53 %), F1 near-neutral. In the mobile rows F2 and F3 are
  real-identical (`F_F3 = F_F2 - B_gov`).

# The labour door: the sectoral family

`N` sectoral labour markets with supply elasticities `eta_s,i` (ADR-0022), whose
**rigid corner is exactly the `eta = 0` endpoint** - the same equations, not a
limit (Proposition 3 of `paper/equivalence.tex`; measured `max|Delta| = 0.00e+00`
warm-started, `< 1e-10` cold, and bit-exact under a supply shock too). Executed
in `matrix_5x3_v9`:

| Sectoral cell | `max abs(p-1)` F1 / F2 / F3 | Employment (F2) | Consumption (F2) | Deflator (F2) | Wage max/min |
| --- | --- | ---: | ---: | ---: | ---: |
| uniform `eta_s = 0.25` | 0.133610 / 0.152897 / 0.152897 | 1.000864 | -1.689% | 1.004096 | 1.38 |
| uniform `eta_s = 0.5` | 0.095839 / 0.105666 / 0.105666 | 1.001255 | -1.571% | 1.002872 | 1.26 |
| uniform `eta_s = 1` | 0.061283 / 0.065422 / 0.065422 | 1.001627 | -1.467% | 1.001802 | 1.16 |
| uniform `eta_s = 2` | 0.035629 / 0.037171 / 0.037171 | 1.001914 | -1.391% | 1.001033 | 1.09 |
| rigid programme sectors, `eta_s = 0.5` | 0.215061 / 0.272033 / 0.272033 | 0.996621 | -2.691% | 1.006206 | 1.70 |
| rigid largest half, `eta_s = 0.5` | 0.109583 / 0.123362 / 0.123362 | 1.000196 | -1.905% | 1.006294 | 1.28 |
| BF (`eta = 0` endpoint) | 0.221442 / 0.279009 / 0.279009 | 1.000000 | -1.980% | 1.007181 | 1.72 |
| GAMMA (fixed wage, contrast) | 0.000000 / 0.000000 / 0.000000 | 1.001260 | -1.533% | 1.000000 | 1.00 |

- **Demand sensitivity, defined strictly.** In every sectoral cell `max abs(p-1)`
  differs between F1 and F2/F3 (F2 = F3 throughout), while all twelve `eta = 1`
  matrix cells read exactly zero. That difference - not the mere fact that prices
  move - is what the door delivers.
- **The elasticity splits the adjustment.** As `eta_s` rises 0.25 to 2 the price
  response falls (0.1529 to 0.0372 at F2), employment rises (1.000864 to
  1.001914) and the welfare cost falls (-1.69 % to -1.39 %).
- **Incidence dominates the level.** Making the seven programme sectors rigid
  more than doubles the price response (0.2720 against 0.1057) and *reverses* the
  employment effect (-0.34 % against +0.13 %).
- **The grouping rule is the open choice**, and it is now motivated rather than
  arbitrary: the programme's incidence is 63.4 % specialised construction works
  plus six material sectors, i.e. exactly the sectors where capacity pressure
  binds, so `eta_s,i = 0` there is a testable statement (`framing_gamma.tex`
  Section 4). The sweep gives the band: price response `[0.097, 0.279]`,
  employment `[-0.31 %, +0.16 %]`.
- **The level of `eta_s` is unidentified under demand-only shocks** - the ladder
  is a sensitivity band. A supply shock separates the rungs (employment
  1.000000 / 1.002187 / 1.008218 at `eta_s` = 0 / 0.5 / 2), so the supply arm is
  the route to a number.

# The GAMMA menu: seven fixed-wage variants

The second route to prices, and the key outcome of the latest round. The
discriminator is the same: "sees demand" means `max abs(p-1)` differs across the
financing columns.

| Variant | Financing | `max abs(p-1)` | Deflator | Employment | Consumption | Sees demand |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| 1. baseline (`w = 1`) | F1 | 0.000000 | 1.000000 | 0.999027 | -0.081% | no |
| | F2 | 0.000000 | 1.000000 | 1.001260 | -1.533% | no |
| | F3 | 0.000000 | 1.000000 | 1.017796 | +2.264% | no |
| 2. DELTA (= baseline) | all | 0.000000 | 1.000000 | = baseline | = baseline | no |
| 3. pinned wage structure (+10 % tilt) | all | 0.057820 | --- | 1.006865 | --- | no |
| 4. capacity channel `delta = +0.5` | F1 | 0.098457 | 1.003027 | 1.004167 | +0.459% | yes |
| | F2 | 0.110761 | 1.005379 | 1.008489 | -1.082% | yes |
| | F3 | 0.128463 | 1.021461 | 1.038814 | +2.199% | yes |
| 5. Kaldor-Verdoorn `delta = -0.5` | F1 | 0.123462 | 1.003454 | 0.997837 | -1.211% | yes |
| | F2 | 0.115629 | 1.001283 | 0.998527 | -2.464% | yes |
| | F3 | 0.128014 | 0.985336 | 1.003370 | +2.266% | yes |
| 6. BF allocation rule (`eta = 0`) | F1 | 0.000000 | 1.000000 | 1.000000 | +0.043% | no |
| | F2 | 0.000000 | 1.000000 | 1.000000 | -1.694% | no |
| | F3 | 0.000000 | 1.000000 | 1.000000 | +0.000% | no |
| 7. dual labour market (insiders rigid) | F1 | 0.024701 | 1.000438 | 0.989569 | -1.338% | yes |
| | F2 | 0.027768 | 1.000492 | 0.991179 | -2.881% | yes |
| | F3 | 0.030189 | 1.000648 | 1.005996 | +0.689% | yes |

- **Two variants see demand, through different margins.** The capacity channel
  through the *cost* (unit cost rises with own output); the dual labour market
  through *rationing* (insider prices 0.0278 against outsider 0.0015 at F2 - a
  19x dualism signature, with identical wages in both segments).
- **The dual labour market is the only fixed-wage variant where the programme
  destroys employment** (-1.04 % / -0.88 % / +0.60 % across F1/F2/F3), because
  the bottleneck is transmitted economy-wide through intermediate costs.
  Quantity rationing is harsher than price-based rigidity (which reads -0.34 %).
- **The capacity channel's sign is the economics.** `delta > 0` is capacity
  pressure (F3 deflator 1.0215), `delta < 0` the Kaldor-Verdoorn case (F3
  deflator 0.9853, prices falling with demand): a 3.6 pp deflator range from one
  parameter's sign.
- **The BF allocation rule is the double-rigidity corner.** Employment is a
  datum, prices are pinned, and under F3 household consumption is literally
  unchanged (+0.000 %). It is the benchmark the other variants are read against.
- **The pinned wage structure buys distribution, not demand** (prices move
  0.058 identically in F1/F2/F3). It is the instrument for *who works where*.

# The exact equivalences

| Equivalence | Status | Source |
| --- | --- | --- |
| ALPHA = BETA | proved ex ante; demand-only scope | Prop. 1, `equivalence.tex` |
| GAMMA = DELTA | definitional (Leontief corner) | Prop. 2, `equivalence.tex` |
| mobile F2 = F3 | proved (`F_F3 = F_F2 - B_gov`) | `equivalence.tex` Section 5 |
| sectoral corner = `eta = 0` endpoint | proved, and shock-independent | Prop. 3, `equivalence.tex` v3 |
| BF vs ALPHA | **not** a coincidence (ADR-0020 option C) | `equivalence.tex` classification |

# Open items

1. **The level of `eta_s`** (the supply-shock design; ADR-0018's schema item) -
   the one remaining batch that converts a band into a number.
2. **The grouping rule** - now motivated by capacity pressure, still to be
   preregistered with its testable evidence (vacancies, overtime, backlogs).
3. \textcolor{revisionV1}{**The programme total** - price basis resolved
   (2026-09-20): the raw impulse is in current (2023) prices (horizon mean
   58.89 bn) and the design's G0 = 40.3 bn is that mean at 2019 prices
   (deflator 1.46, construction-cost based; 58.89/1.46 = 40.34 bn). The
   remaining question is the composition (the raw 2024 current-price vector)
   and the manuscript wording; ADR-0024 records the evidence.}
4. **The welfare metric** - dispersion-blind, and the kernel has no leisure term
   or income effect, so employment and wage dispersion are allocation facts, not
   welfare claims.
5. **The manuscript pass** and the `docs/DOCS_ASSESSMENT.md` block (still citing
   v3/v4 numbers).

# Where the numbers come from

- Runs: `runs/matrix_5x3-v9-*` (33 cells; the matrix rows are bit-identical to
  `matrix_5x3-v6-*`). Flow tables: `paper/tables/matrix_5x3_v9_flows.md` and
  `matrix_5x3_v6_flows.md`.
- Probes: `probe13` / `probe14` (the sectoral closure and the S1-S5 slots),
  `probe15` (the supply arm), `probe16` (the rigidity band), `probe17` (the
  capacity channel in the fixed-wage closure), `probe18` (the BF allocation rule
  in GAMMA), `probe19` (the dual labour market).
- Decisions: ADR-0019 (external account), ADR-0020 (the `eta = 0` endpoint),
  ADR-0021 (wage structure, proposed), ADR-0022 (sectoral labour markets,
  accepted). Derivations: `paper/equivalence.tex` v3, `paper/framing_gamma.tex`.
  Steps and decision points: `docs/WORKPLAN_SENSITIVE_PRICES.md` v7.

# Revision Log

- **Version 2** \textcolor{revisionV1}{(September 2026)} --- Open item 3
  updated: the programme total's price basis is resolved (the raw impulse is
  in current prices; G0 = 40.3 bn is the horizon mean at 2019 prices,
  deflator 1.46; 58.89/1.46 = 40.34 bn), per ADR-0024. No numbers change.
- **Version 1** (September 2026)
