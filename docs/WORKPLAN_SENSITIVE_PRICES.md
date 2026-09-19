---
title: "Workplan: demand-sensitive prices from the wage block (segmented wage setting)"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-09-18"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [workplan, labour-closure, segmented-wages, prices, adr-annexe, handoff]
last-updated: September 2026
---

**Version 3** \textcolor{revisionV2}{(September 2026)}
**Version 2** \textcolor{revisionV1}{(September 2026)}
**Version 1** (September 2026)

This is the implementation annexe to ADR-0021 and to the segmented-wage ADR it
calls for. It lives in `docs/` at the operator's request, as the executable
companion to those decisions: the ADRs carry the choice, this file carries the
steps, the tests and the decision points for the session that executes it.
Nothing here is a decision, and nothing here has been run.

# The goal, and what it excludes

Goal: **demand-only scenarios in which real prices respond to demand.** The
paper's shocks stay where they are (the programme impulse, the F1 tilt, the F2
tax, the F3 external financing); no supply-side shock is introduced. The price
response is to be generated inside the wage block instead.

Excluded, with reasons: supply-side shocks (awkward in a demand-shock framing,
recorded thematically as `docs/ideas/IDEA-0001-climate-productivity-shocks.md`),
imperfect competition and markups (a different price block), a second production
factor (a different model), and the retired allocation wedge (never derived as a
coefficient of the model).

# The measurement that defines the problem

From the `matrix_5x3_v6` manifests, the price response of every labour row:

| Row | `max abs(p - 1)` (F1 / F2 / F3) | GDP deflator | Wage |
| --- | --- | --- | --- |
| BF (eta = 0, sectoral wages) | 0.221 / 0.279 / 0.279 | 1.006326 / 1.007181 / 1.007181 | solved, 0.9734 to 1.5249 |
| ALPHA, BETA (eta = 1) | 2.7e-15 / 2.7e-15 / 2.7e-15 | 1.000000 | 1.000000 |
| GAMMA, DELTA (fixed wage) | 1.9e-13 / 2.2e-15 / 2.2e-15 | 1.000000 | 1.000000 |

Two conclusions. First, the `eta = 0` row is already a demand-sensitive-price
row: its sectoral wages are solved given the demand composition, so prices and
the deflator move with the financing cell. The paper can say something about
prices in the demand-only design today, from that row, with no new machinery.
Second, in every other row prices are pinned to the baseline to machine
precision, and the reason is structural rather than numerical: with one factor,
constant returns and a demand-free price block, $p_i = cost_i(p, w, ip)$ depends
on the wage and on intermediate prices only. Demand can reach prices **only
through an endogenous wage**. That is why the fix belongs in the wage block, and
it also explains why BETA and DELTA coincide with their neighbours: with prices
fixed, the labour-supply elasticity has nothing to bite on and the CES-versus-
Leontief contrast disappears.

One clarification of the earlier option list, because it matters for scoping: a
pinned wage **vector** (Door 1, ADR-0021) produces price *heterogeneity* --- a
tilted pin gives `max abs(p - 1)` of order 0.1 --- but not price *sensitivity*:
the pin is exogenous, so the price block still never sees demand. Door 1 is the
complement (it makes the allocation margin observable); the route to
demand-sensitive prices is a **free** wage, which is Door 3.

# The specification to implement

Partition the sectors into a sticky set $S$ (wages pinned) and a flexible set
$F$ (one common wage $w_f$, solved). Unknowns
$[p(1:N); y(1:N); w_f; F_{ext}]$, that is $2N + 2$; equations:

- $N$ zero-profit conditions, with the sector wage $w_i = \bar w_i$ for
  $i \in S$ and $w_i = w_f$ for $i \in F$;
- $N$ goods-market clearings, all of them, as in every regime since ADR-0019;
- one labour-market condition for the flexible segment,
  $\sum_{i \in F} L^{cm}_i = \bar L_f \cdot ((w_f / P) / (\bar w_f / \bar P))^{\eta_s}$,
  so the sticky segment's employment stays demand-determined and uncapped, as in
  the fixed-wage regime;
- one numeraire, $P^{CPI} = 1$.

Properties to rely on, and to check rather than assume:

- **Nesting.** $S = \emptyset$ reproduces ALPHA (and BETA at $\eta_s = 0.5$);
  $S$ = all sectors reproduces GAMMA, where the free wage disappears and the
  external unknown is dropped. Three existing closures become endpoints of one
  family, which is the presentational prize: GAMMA stops being a row and becomes
  the "all pinned" corner.
- **Exogenous parameter and endogenous outcome.** Only the sticky *set* is
  exogenous; the *level* of the pinned wages is a numeraire (a uniform rescale of
  $(\bar w_S, w_f)$ leaves the real allocation invariant, as measured for the
  pure pin in `probe11`). The **segment wage gap** $w_f / \bar w_S$ is therefore
  an equilibrium outcome, not a knob --- which is the result worth reporting.
- **Demand-sensitive prices.** $w_f$ appears in the clearing block through
  $L^{cm}$, so it responds to the demand composition, and with it every price
  through the intermediate block.
- **Determinacy.** The ADR-0017 round-gain criterion applies to the fully pinned
  system and does **not** carry over unchanged once a wage is free; the free wage
  plausibly supplies the feedback the fully pinned system lacks, so the
  admission question should relax --- but this must be measured, not asserted.

# Scenarios to evaluate \textcolor{revisionV2}{\normalsize [added v3]}

\textcolor{revisionV2}{The specification above fixes a system; this section fixes what is run and what is reported, so that the next session can preregister a design without re-deciding the experiment. The slots cover both specifications, because the choice between them (next section) is a choice about what wage rigidity does to employment, and it should be made against evidence rather than in advance.}

| Slot | Closure and parameters | Cells | What it demonstrates | What to report |
| --- | --- | --- | --- | --- |
| `S0` | the existing fifteen cells, re-run for provenance | `matrix_5x3-v6-*` | the baseline, and the rows in which prices do not move | the price table as it stands: `max abs(p - 1)` and the deflator |
| `S1` | `eta = 0` (BF): all-sectoral rigid wages at the frozen allocation | `BF-F1`, `BF-F2`, `BF-F3` (already executed and gated) | the largest demand-driven price response in the family, with **no elasticity parameter at all**: the wages solve the first-order conditions at the frozen allocation | `max abs(p - 1)` of 0.221 to 0.279, the deflator, the wage dispersion, the external position |
| `S2` | general closure, uniform elasticity ladder: `eta_s,i = eta_s` for every sector | `eta_s` in 0.1, 0.5, 1, 2, each across F1, F2, F3 | the central exhibit of the general formulation: the price response as a function of the sectoral supply elasticity | the pass-through measure (below), plus employment and welfare at each rung |
| `S3` | general closure, two-group vector: a sticky group at `eta_s,i = 0`, a flexible group at `eta_s = 0.5` | sticky group in {the programme sectors, the largest half by baseline employment}, each across F1, F2, F3 | the rigidity share as a continuous dimension, with the segment wage gap as an outcome | the segment wage gap, employment by group, and the price response against the sticky share |
| `S4` | pinned wage vector (Door 1, ADR-0021) | the BF-transplant pin plus two institutional pins, each across F1, F2, F3 | price *heterogeneity* without price *sensitivity*, and the allocation-margin comparison against `S1` | prices, and the allocation difference against the BF cells that share the same wage structure |
| `S5` | optional companion: the utilization externality, `delta` in {+0.5, -0.5} | the programme cells | the opposite price sign (the Kaldor--Verdoorn case at `delta < 0`) | the deflator under both signs, as a bracket on the price response |

## Controls, not scenarios

\textcolor{revisionV2}{These belong in the test suite as much as in the design, and they are the implementation's canaries:}

- \textcolor{revisionV2}{**The nesting check.** The general closure at `eta_s,i = 0` for every sector must reproduce the `S1` manifests, because those are the same equations. This is the canary for the translation: if `S2` at `eta_s = 0` does not equal `S1`, the kernel change is wrong. It replaces the ALPHA and GAMMA nesting of the restricted specification, which the general one does not have.}
- \textcolor{revisionV2}{**The neutrality control.** `F2` and `F3` must have identical prices, wages and allocations, as in every regime since ADR-0019; a difference beyond solver noise means the new closure broke the accounting.}
- \textcolor{revisionV2}{**The level control.** A uniform rescale of the anchor wages must leave the real allocation, employment and the price ratios unchanged (Corollary 1 of `paper/equivalence.tex`).}

## What the general formulation captures that the restricted one does not

\textcolor{revisionV2}{This is the evidence the specification choice should rest on, and it is why the slots above are written for the general closure:}

- \textcolor{revisionV2}{**The elasticity ladder.** `S2` traces the price response from its maximum at `eta_s,i = 0` to zero as sectoral supply becomes perfectly elastic. That is a statement about the model, not about one parameter choice; the restricted variant has no such ladder, its only parameter being which sectors are sticky.}
- \textcolor{revisionV2}{**What rigidity does to employment.** In the general closure a sector at `eta_s,i = 0` has its employment *fixed* at the baseline, because the labour-market condition pins `L^{cm}_i = L-bar_i`. In the restricted variant a sticky sector's employment is *demand-determined and uncapped*, as in the fixed-wage regime. Both are defensible closures, but they are different economies, and that difference is why the choice between the two specifications is substantive rather than cosmetic.}
- \textcolor{revisionV2}{**A free endpoint that already exists.** The general closure's all-rigid corner is exactly the `eta = 0` row, so the family is anchored on an executed, gated and published set of cells rather than on a new construction.}

## The reporting unit

\textcolor{revisionV2}{The paper needs one number per cell, not a price vector. The natural measure is a **pass-through**: the change in the GDP deflator (or in `max abs(p - 1)`) per percentage point of programme spending relative to GDP. It is comparable across labour rows, it is defined in every regime, and for the existing cells it is computable from the manifests --- zero to machine precision in ALPHA, BETA, GAMMA and DELTA, and of order 0.5 in the BF row, where the deflator moves 0.63 to 0.72 per cent against a programme of 1.33 per cent of GDP. Every slot above should be reported with it, so that `S2` and `S3` read as a curve in one parameter rather than as a table of vectors.}

# Reconciliation with the five-door menu of `docs/ETAs.md` \textcolor{revisionV1}{\normalsize [added v2]}

`docs/ETAs.md` derives the same obstruction this workplan starts from --- Lemma 1
of `paper/equivalence.tex` (the price block is homogeneous of degree one in
$(p, w)$ and demand-free, so $p = w = \mathrm{CPI} = 1$ under $A = 1$) and
Corollary 1 (the real wage is pinned at its anchor) --- and ranks the doors that
lift it. Four consequences for this workplan, none of which changes its priority.

- **The labour door is candidate 1 there too.** The sequencing note in
  `docs/ETAs.md` keeps sector-specific wages as candidate 1 for the published
  matrix, ahead of the utilization/Verdoorn arm and the terms-of-trade arm, with
  the supply-side arm first only because it is cheap and earns the `eta_s`
  identification claim. This workplan pursues the labour door, which is the only
  door that repairs the demand-only matrix; the supply arm is parked for framing
  reasons (`docs/ideas/IDEA-0001-climate-productivity-shocks.md`). The deviation
  from ETAs' ordering is deliberate and is recorded here rather than taken
  silently.
- **Specification: general or restricted.** `docs/ETAs.md` specifies the labour
  door as $N$ sectoral labour markets with a vector of sectoral supply
  elasticities; the specification above is the restricted variant (a sticky set
  plus one free wage), chosen because it nests ALPHA and GAMMA. The two are
  different systems: the general one makes every sectoral wage demand-sensitive
  and its all-rigid corner ($\eta_{s,i} = 0$ everywhere, so $L^{cm}_i = \bar L_i$)
  reproduces the `eta = 0` row, but it does not nest ALPHA or GAMMA; the
  restricted one nests those two rows but makes only one wage demand-sensitive. \textcolor{revisionV2}{They also differ in what rigidity does to employment, which is the substantive part: at `eta_s,i = 0` the general closure pins `L^{cm}_i = L-bar_i`, so a rigid sector's employment is fixed, while in the restricted variant it is demand-determined and uncapped. See the scenario section above for what each captures.}
  **Recommended: implement the general closure**, with the scenarios reported at
  a low-dimensional parametrisation (a sticky group at $\eta_{s,i} = 0$ and a
  flexible group at a common elasticity), because $N$ elasticities are weakly
  identified under demand-only shocks. The restricted variant then becomes a
  presentation choice rather than a separate implementation.
- **The Verdoorn arm is the companion, not a rival.** The capacity door's reduced
  form --- an endogenous utilization externality in the unit-cost hook, with
  $\delta < 0$ giving the Kaldor--Verdoorn case --- is the cheapest door and the
  only one with a *negative* price sign, so running it alongside the labour door
  brackets the price response instead of predicting one direction. It is recorded,
  with the markup, capacity and external doors, as
  `docs/ideas/IDEA-0002-five-doors-demand-sensitive-prices.md`.
- **One schema change serves both routes.** The per-cell parameter source that
  `docs/ETAs.md` records as the blocking item of its own workplan (the harness
  hard-codes a null shock and exposes no per-cell parameter) is the same change
  this workplan needs for the utilization arm and for the parked supply route. It
  should be done once --- and the identifier should be checked first: the id
  `ADR-0018` cited in `docs/ETAs.md` is already taken by the real-GDP measurement
  decision, so a supply-shock schema ADR must take the next free number.

# Steps

- **Step 1, prototype (no `src/` change).** Write
  `experiments/probes/probe12_segmented_wages.jl`, mirroring the fixed-wage
  residual for the sticky sectors and the mobile residual for the flexible ones,
  reusing `_mobile_market_demand` with a full per-sector wage vector (it accepts
  one since ADR-0020). Measure, on the full-71 calibration: the endpoint nesting
  against the committed `ALPHA` and `GAMMA` manifests; the demand sensitivity of
  prices across the three financing cells; the segment wage gap and employment
  as the sticky share rises; the identity gap through the price-weighted clearing
  residual (the proxy validated in `probe11`, ratio 1.0000 against the kernel
  canary along quantity perturbations); the Jacobian conditioning and the
  residual floor. Record the falsifiable predictions before running: prices move
  with the financing cell; the wage gap widens with the sticky share; employment
  absorbs more as the sticky share rises.
- **Step 2, ADR.** Draft the segmented-wage ADR (proposed) with the
  specification above, the nesting table, the measured demand sensitivity and the
  admission check. Registry: a new labour-closure entry, with `labor.GAMMA`
  pointing at it as the all-pinned endpoint if the operator decides to present
  the family instead of the row.
- **Step 3, promotion.** Dispatch-only kernel change: a segmented closure type
  carrying the sticky set and the pinned vector, its residual function, the
  `solve` plumbing, and the vector-wage handling in `equilibrium_residuals`,
  `market_clearing_residuals`, `external_balance_canary` and `gdp_components`.
  The harness needs a third gate for the segmented regime (the flexible
  segment's labour-market residual) and, as at eta = 0, the scalar `wage` metric
  should become the wage-bill-weighted aggregate with `wage_min` and `wage_max`
  carrying the structure.
- **Step 4, tests.** The endpoint nesting on the fixtures and on one full-71
  cell; the level invariance under a uniform rescale; the demand-sensitivity
  assertion (the price deviation must move with the financing cell, which no
  existing closure satisfies); the identity gate; the determinacy check.
- **Step 5, design and generation.** A design carrying the existing fifteen cells
  plus the segmented cells (sticky share at several levels, each across F1, F2,
  F3), preregistered before execution. A `src/` change supersedes the current
  generation's provenance, so the existing cells are re-run in the same batch
  rather than mixed with it.
- **Step 6, records.** Flow table from the new manifests, the assessment block on
  demand-sensitive prices, the registry and log entries, the board, the gate.

# Decision points for the operator

- **The sticky set.** Which sectors are pinned, and by what rule. This is the
  research choice, and it should be preregistered; a rule ("the programme
  sectors", "the largest half by baseline employment") is preferable to a
  hand-picked list, because a hand-picked list invites the tuning objection.
- **The flexible segment's supply elasticity.** $\eta_s = 0$ makes the family
  nest ALPHA at one end; $\eta_s = 0.5$ nests BETA. Either is defensible; the
  choice determines what the paper can say about the supply side.
- **GAMMA's role.** Keep it as a matrix row for continuity, or present the family
  with GAMMA as its all-pinned endpoint. The second is cleaner and answers the
  redundancy question, but it changes the paper's table structure.
- **Scope.** Whether the rigidity share becomes a headline dimension of the
  revision or a follow-up paper. Nothing in steps 1 to 3 commits either way; the
  probe is worth doing regardless.

# Risks and how they surface

- **The mirrored prototype.** The eta = 0 work showed that a correct prototype
  can still be mistranslated into the kernel, where the defect is structural and
  silent. The endpoint nesting is the self-diagnosing check: if $S = \emptyset$
  does not reproduce ALPHA and $S$ = all does not reproduce GAMMA, the
  translation is wrong.
- **Determinacy.** Measured in step 1, not assumed; if the free wage does not
  stabilise the system, the family is admitted only on a restricted sticky set.
- **Provenance.** A new generation is required; the existing cells must be
  re-run in the same batch and reproduce their manifests.
- **Interpretation.** The sticky set is exogenous and the wage gap is an outcome.
  A welfare difference across sticky shares is not a statement about wage policy,
  and the rigidity share does not identify a labour-supply elasticity.

# Evidence base

`experiments/probes/probe11_gamma_wage_structure.jl` (level invariance, the
identity proxy and its validation, conditioning); `probe7_sectoral_wages_eta0.jl`
(the eta = 0 closure the vector-wage apparatus comes from);
`probe8_promotion_verification.jl` and `probe9_nonbf_reproduction.jl` (the
promotion's verification pattern); `runs/matrix_5x3-v6-*` (the price table
above); `docs/VariationinGamma.md` (the four doors and the claim-by-claim
provenance); ADR-0020 and ADR-0021; the `eta_s` identifiability item in
`docs/DOCS_ASSESSMENT.md`.

# Revision Log

- **Version 1** (September 2026)
- **Version 2** \textcolor{revisionV1}{(September 2026)} --- Added the reconciliation with the five-door menu of `docs/ETAs.md`: the labour door is that document's candidate 1 as well, so the priority stands; the specification choice between the general sectoral-labour-market closure and the restricted sticky-set variant is flagged with a recommendation (general, reported at a low-dimensional parametrisation, because N sectoral elasticities are weakly identified under demand-only shocks); the Verdoorn arm is recorded as the cheap companion with the opposite price sign; and the shared per-cell parameter source is identified as one schema change serving both routes, with the note that the `ADR-0018` id cited in `docs/ETAs.md` is already taken.
- **Version 3** \textcolor{revisionV2}{(September 2026)} --- Added the scenario section (section 4, with its three subsections), which the operator flagged as missing: six slots (`S0` the re-run baseline, `S1` the already-executed `eta = 0` row, `S2` the uniform elasticity ladder, `S3` the two-group rigidity share, `S4` the pinned-vector complement, `S5` the optional Verdoorn companion), three controls that act as implementation canaries (the nesting check at `eta_s,i = 0`, the F2/F3 neutrality control, the level control), an explicit statement of what the general formulation captures that the restricted one does not, and a single reporting unit (a pass-through per percentage point of programme spending) so that the slots read as a curve rather than a table of vectors. The specification bullet of the reconciliation section gained the employment-behaviour difference between the two variants.