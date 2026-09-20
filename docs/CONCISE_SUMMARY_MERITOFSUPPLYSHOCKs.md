---
title: "What the supply-shock exercise adds to the paper"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-09-20"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [summary, supply-arm, identification, eta_s, merit, paper]
last-updated: September 2026
---

**Version 1** (September 2026)

One question, answered plainly: what did running 189 supply-shock cells
(`supply_etas_s1/prog/unif`, all executed, all gates pass) buy the paper?
Short answer: it turns a *defect-looking degeneracy* into a *measured
identification result*, it prices the elasticity that the matrix cannot
see, and it restores the BETA row as an informative part of the closure
family. It does not change a single demand-side number.

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

- **Version 1** (September 2026)