---
title: "Concise summary: the 5x3 matrix, the sectoral labour door, and the GAMMA menu"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-09-19"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [summary, matrix, labour-closure, prices, gamma, handoff]
last-updated: September 2026
---

**Version 3** \textcolor{revisionV2}{(September 2026)}
**Version 2** \textcolor{revisionV1}{(September 2026)}
**Version 1** (September 2026)

One page of results, no derivation: the executed 5x3 matrix, the sectoral
labour-market family that gives the flexible-wage rows a demand channel, and the
seven fixed-wage (GAMMA) variants measured side by side. Everything below is a
manifest field or a probe measurement; the sources are listed at the end, and
the derivations live in `paper/equivalence.tex` (v3), `paper/framing_gamma.tex`
and the wage-structure door in `docs/VariationinGamma.md` (ADR-0021).

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
in `matrix_5x3_v9` (deflators corrected from `matrix_5x3_v10`):

| Sectoral cell | `max abs(p-1)` F1 / F2 / F3 | Employment (F2) | Consumption (F2) | Deflator (F2) | Wage max/min |
| --- | --- | ---: | ---: | ---: | ---: |
| uniform `eta_s = 0.25` | 0.133610 / 0.152897 / 0.152897 | 1.000864 | -1.689% | 1.004088 | 1.38 |
| uniform `eta_s = 0.5` | 0.095839 / 0.105666 / 0.105666 | 1.001255 | -1.571% | 1.002868 | 1.26 |
| uniform `eta_s = 1` | 0.061283 / 0.065422 / 0.065422 | 1.001627 | -1.467% | 1.001800 | 1.16 |
| uniform `eta_s = 2` | 0.035629 / 0.037171 / 0.037171 | 1.001914 | -1.391% | 1.001033 | 1.09 |
| rigid programme sectors, `eta_s = 0.5` | 0.215061 / 0.272033 / 0.272033 | 0.996621 | -2.691% | 1.006190 | 1.70 |
| rigid largest half, `eta_s = 0.5` | 0.109583 / 0.123362 / 0.123362 | 1.000196 | -1.905% | 1.006281 | 1.28 |
| BF (`eta = 0` endpoint) | 0.221442 / 0.279009 / 0.279009 | 1.000000 | -1.980% | 1.007181 | 1.72 |
| GAMMA (fixed wage, contrast) | 0.000000 / 0.000000 / 0.000000 | 1.001260 | -1.533% | 1.000000 | 1.00 |

Note (2026-09-20): the deflators above are corrected from the
`matrix_5x3_v10` manifests — the previous values were measured with the C1
defect (`gdp_components` collapsed the sectoral wage vector; ADR-0022
amendment; `paper/tables/matrix_5x3_v10_flows.md`). All other columns
(prices, employment, consumption, wage dispersion) are solutions and are
unchanged between v9 and v10.

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
3. \textcolor{revisionV2}{**The programme total** - settled (2026-09-20):
   the principal shock is the deflated horizon mean, G0 = 40.3 bn (40.34 bn
   to 0.1 %; the raw impulse is in current (2023) prices, deflator 1.46);
   ADR-0024 is accepted. The remaining question is the composition (the raw
   2024 current-price vector; a constant-2019-price incidence vector needs
   sectoral deflators) and the manuscript wording.}
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
  accepted), \textcolor{revisionV1}{ADR-0024 (programme total; price basis: 2023
  prices per Hornykewycz et al. 2025; G0 = 40,300 at 2019 prices)}.
  Derivations: `paper/equivalence.tex` v3, `paper/framing_gamma.tex`.
  Steps and decision points: `docs/VariationinGamma.md` (ADR-0021) and `docs/ETAs.md`.


# Supply-shock exercise (merged from `CONCISE_SUMMARY_MERITOFSUPPLYSHOCKs.md`, ADR-0027)
# The problem the matrix left on the table

The published 5x3 matrix is demand-only: a green-investment programme
shifted final demand, nothing else. Under demand-only shocks the model's
price block is pinned by technology alone (the zero-profit conditions are
homogeneous in prices and the wage, contain no demand term; Lemma 1 of
`paper/equivalence.tex`), so relative prices stay at one and the real wage
sits at its anchor. Two consequences look like defects:

- **BETA reproduces ALPHA bit-for-bit.** The elastic-supply row returns
  exactly the full-employment row, because with the real wage pinned the
  supply curve returns baseline employment whatever its elasticity.
- **`eta_s` (the labour-supply elasticity) is unidentified.** Every value
  0.5, 2, 5 gives the same solution to machine precision. Reviewer
  point **R2.1** (labour-supply elasticity as bridging parameter) is
  precisely this demand.

The supply-shock exercise is the experiment that resolves the demand: the
paper always knew the only channel that moves the real wage is technology
(DE-0011 documented that the numeraire is *not* the lever). The arm makes
technology shocks a first-class scenario and measures the wage-employment
response.

# What was measured

Three shock shapes, on the same full-71-sector calibration and the same
40.3 bn programme (G0 = 40,300 EUR m at 2019 prices, ADR-0024), 63 cells
each, F1/F2/F3 financing, identification matrix ALPHA + scalar BETA
(`eta_s` = 0.5/1/2/5) + sectoral uniform BETA:

| Design | Shock | Wage (F2, mid ladder) | max abs(p-1) | Implied `eta_s` |
| --- | --- | ---: | ---: | --- |
| `supply_etas_s1` | +10/20/30 % productivity in sector 1 | 1.0049 | 0.196 | exact (0.5000/1.0000/2.0000/5.0000) |
| `supply_etas_prog` | +5/10/20 % on the programme's own sectors | 1.0046 | 0.0625 | exact |
| `supply_etas_unif` | +1/2/3 % in every sector | 1.0408 | 0.0446 | exact |

`ln L / ln w` recovers the input elasticity to four decimals in every
scalar BETA cell and is invariant to the shock magnitude (spread < 1e-6
for s1, < 1e-3 for prog/unif). The ALPHA control pins employment at
exactly 1 with the wage moving. All three arms satisfy the model's own
accounting: F2 = F3 financing neutrality and the closed external-account
identity hold under the technology shock on all 189 cells (the identity
had previously only ever been asserted at prices equal to one).

# What it adds

**1. The identification result --- `eta_s` is a real parameter, and the
null was the experiment, not the model.** The map `eta_s` $\to$
(wage, employment, real GDP) is strictly monotone and invertible under a
supply shock, and the model recovers the input elasticity *exactly* from
its own response. The demand-only degeneracy is now explained as a
property of the demand-only *design* (Lemma 1), not a hidden failure of
the closure. That is the direct, measured answer to R2.1.

**2. The sensitivity of the published results to the elasticity.** The
band that the matrix could not resolve is now a quantitative statement:
under one representative technology shock the employment effect spans
+0.24 % (`eta_s` = 0.5) to +2.48 % (`eta_s` = 5) --- a 2.2 percentage-point
spread attributable to the labour-supply assumption alone, and real GDP
spans a similar range. The paper can state how much of its headline rests
on the assumed elasticity per financing column, instead of leaving the
parameter as an undisclosed dial.

**3. BETA becomes an informative row.** Under a technology shock BETA no
longer equals ALPHA: employment rises with `eta_s` while full employment
stays flat, and prices move. The closure family shows variety exactly
where the theory says it should. The "BETA is a free rider / the model
cannot distinguish closures" objection loses its force.

**4. The distribution of a productivity dividend between wages and
jobs.** The uniform arm is the cleanest reading: a 2 % economy-wide
productivity gain raises the real wage 4 % and, depending on the
labour-supply regime, yields real GDP growth of about +8 % (`eta_s` = 1)
to +27 % (`eta_s` = 5) as employment expands from +4 % to +22 %. In the
paper's two-paradigm framing: a full-employment labour market (ALPHA)
takes the whole dividend in wages; an elastic supply shares it with jobs.
The exercise quantifies the elasticity of that distribution --- a
distributional statement the demand-only matrix could not touch.

**5. A green-productivity scenario.** The programme-own-sector arm
(`prog`) asks: what if the decarbonisation programme itself raises
productivity in construction and materials (learning, prefabrication)
by up to 20 %? The answer includes a supply-side offset to the cost side
--- prices respond (max abs(p-1) = 0.0625) and real GDP rises above the
demand-only picture. That ties the exercise to the paper's motivation
(investment in foundational infrastructures) and to the climate-
productivity idea recorded as IDEA-0001.

**6. Architectural robustness, demonstrated, not assumed.** Two
properties the paper already leans on --- financing neutrality (F2 = F3)
and the closed external account --- are verified under a second family of
shocks on 189 cells. The F2/F3 finding is not an artifact of the
demand-only design; the accounting is not fragile to $A \neq 1$.

**7. The recombination, realized.** The matrix deliberately left out the
($\eta$, `eta_s`) corner (immobile allocation with elastic supply) as a
demand-only null result. Under supply shocks that corner is informative,
and the arm separates the two margins empirically: the allocation
friction moves composition, the supply elasticity moves the level ---
exactly as the note anticipated.

# What it could add, if drafted into the manuscript

Concrete, ready-draftable paragraphs (decide which earn their place):

- **The R2.1 referee paragraph**: one page stating the identification
  result, the Lemma-1 explanation of the demand-only equality, and the
  measured recovery table --- with the flow tables cited by run id.
- **The disclosure paragraph**: the `eta_s`-sensitivity band, reported
  per financing column, with the sentence "the employment spread across
  the elasticity grid [X pp] is the part of the result that rests on an
  assumption the demand-only data cannot pin down".
- **The distributional reading** (item 4 above) --- the strongest
  narrative hook for the heterodox framing: labour-market institutions
  determine who gets the productivity dividend.
- **The green-productivity scenario** (item 5) as a supplementary
  scenario alongside IDEA-0001.

# What it does not add (so the claims stay honest)

- It does **not** repair the demand-only matrix: prices stay at one there
  by Lemma 1, and demand-sensitive prices in *that* matrix remain the
  job of the labour door and the capacity channel --- the supply arm is a
  different experiment, clearly labelled as such.
- It does **not** estimate a real-world `eta_s` from data. It builds the
  identification machinery and measures the band; an empirical elasticity
  (from an observed wage-employment response, e.g. a construction-sector
  shock) enters the model as a scenario, as it should.
- The `unif` magnitudes are partly architecture: the constant-returns
  open structure lets elastic supply absorb productivity gains
  aggressively. Report them as levels, do not over-read the multiplier.
- None of the demand-side matrix numbers changed anywhere (the generation
  re-ran bit-identical).

# Where the numbers live

Flow tables `paper/tables/supply_etas_{s1,prog,unif}_flows.md` (generated
from the manifests, no hand-typed numbers); runs `runs/supply_etas_s1-*`,
`runs/supply_etas_prog-*`, `runs/supply_etas_unif-*`; the workplan and
signatures in `docs/ETAs.md` v5; the assessment note in
`docs/DOCS_ASSESSMENT.md` v8 (finding 2 upgrade); the schema decision in
ADR-0023, the programme-size settlement in ADR-0024.

# Revision Log

- **Version 3** \textcolor{revisionV2}{(September 2026)} --- Open item 3
  updated: the operator keeps the deflated horizon mean (G0 = 40.3 bn) as the
  principal programme; ADR-0024 accepted. The composition remains the open
  piece. The six sectoral `Deflator (F2)` values are corrected from the
  `matrix_5x3_v10` manifests (the v9 values carried the C1 defect); no other
  numbers change.
- **Version 2** \textcolor{revisionV1}{(September 2026)} --- Open item 3
  updated: the programme total's price basis is resolved (the impulse table is
  in current (2023) prices per Hornykewycz et al. 2025; G0 = 40.3 bn is its
  horizon mean at 2019 prices, deflator 1.46; 58.89 / 1.46 = 40.338 bn =
  40.337 bn), per ADR-0024 (accepted). No numbers change.
- **Version 1** (September 2026)

