---
title: "Workplan: demand-sensitive prices from the wage block (sectoral labour markets)"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-09-18"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [workplan, labour-closure, sectoral-wages, prices, adr-annexe, handoff]
last-updated: September 2026
---

**Version 5** \textcolor{revisionV4}{(September 2026)}
**Version 4** \textcolor{revisionV3}{(September 2026)}
**Version 3** \textcolor{revisionV2}{(September 2026)}
**Version 2** \textcolor{revisionV1}{(September 2026)}
**Version 1** (September 2026)

This is the implementation annexe to ADR-0021 and to the sectoral-labour-market
ADR it calls for. It lives in `docs/` at the operator's request, as the executable
companion to those decisions: the ADRs carry the choice, this file carries the
steps, the tests and the decision points for the session that executes it.
\textcolor{revisionV3}{Version 3 stated "nothing here has been run"; Version 4
records that Step 1 has now been run (probe 12) and that its result replaced the
restricted sticky-set specification with the general sectoral-labour-market
closure. The general closure is now the specification; the restricted variant is
kept below as a measured dead end.}

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

\textcolor{revisionV3}{The $eta = 0$ row already realizes the general
construction proposed below, at its all-rigid corner: there the per-sector
first-order condition ties each wage to that sector's quantity, which is exactly
what lets demand move prices. What Version 4 adds is the elasticity that makes
that tie flexible, $\eta_{s,i} > 0$, so that the same machinery traces a price
response from its maximum (the $eta = 0$ row) to zero.}

One clarification of the earlier option list, because it matters for scoping: a
pinned wage **vector** (Door 1, ADR-0021) produces price *heterogeneity* --- a
tilted pin gives `max abs(p - 1)` of order 0.1 --- but not price *sensitivity*:
the pin is exogenous, so the price block still never sees demand. Door 1 is the
complement (it makes the allocation margin observable); the route to
demand-sensitive prices is a **free** wage, which is Door 3.

# The specification to implement \textcolor{revisionV3}{\normalsize [revised v4: the general sectoral-labour-market closure]}

\textcolor{revisionV3}{Version 3 specified a \emph{restricted} sticky-set variant (a
sticky set $S$ with pinned wages and the remainder sharing one solved wage
$w_f$) and named it "the specification to implement". Step 1 of Version 3 was
run as \texttt{experiments/probes/probe12\_segmented\_wages.jl} and
\textbf{refutes that variant}: with a single free wage, zero profit ($N$
equations) together with
the CPI numeraire (one equation) already determine the $N$ prices and $w_f$, so
the flexible labour-market condition never reaches the wage --- it only fixes
the flexible segment's employment. Measured on the full-71 calibration: the
restricted system nests ALPHA/BETA bit-exactly at $S = \emptyset$ and GAMMA
bit-exactly at $S$ = all, its external account closes to $\sim 10^{-16}$, and
employment responds to the sticky share --- but $w_f \equiv 1$ and
$\max|p - 1| \sim 10^{-14}$ in every financing cell, and a tilted pin (+10\% on
the programme sectors) moves prices by $5.5 \times 10^{-2}$ \textbf{identically} across
F1, F2 and F3. That is Door-1 heterogeneity, not demand-sensitivity. The
restricted variant therefore stays only as a measured dead end; the
specification is the general closure.}

Partition is not used. Every sector has its own labour market with its own
sectoral real-wage supply curve:

- Unknowns $[\,p(1{:}N);\; y(1{:}N);\; w(1{:}N);\; F\,]$, that is $3N + 1$
  \textcolor{revisionV3}{(the same unknown count as the $eta = 0$ endpoint, ADR-0020 option C)}.
- $N$ zero-profit conditions, $p_i = cost_i(p, w_i)$ (a sectoral wage per
  sector);
- $N$ goods-market clearings, all of them, with the cost-minimizing allocation
  $L_i = L^{cm}_i(p, y_i, w_i)$ and household wage income $\sum_i w_i L_i$;
- $N$ sectoral labour-market conditions,
  $L^{cm}_i(p, y_i, w_i) = \bar L_i\,\big[ (w_i / \Pi(p)) / (\bar w_i / \bar \Pi) \big]^{\eta_{s,i}}$,
  with a **vector** of sectoral supply elasticities $\eta_{s,i} \geq 0$ and the
  anchors $\bar w_i / \bar \Pi = 1$ at the calibration;
- one numeraire, $\Pi(p) = 1$;
- the scalar $F$ is free, exactly as at $eta = 0$: replacing the mobile
  system's single aggregate labour equation by $N$ sectoral conditions adds
  $N - 1$ equations to a block that is homogeneous of degree one in
  $(p, w, F)$, so one free scalar is required for consistency with the demand
  block (dropping a clearing equation instead is the retired ADR-0010 shortcut).

Properties to rely on, and to check rather than assume:

- **Nesting onto an executed row.** $\eta_{s,i} = 0$ for every sector pins
  $L^{cm}_i = \bar L_i$, so the allocation is the frozen baseline and each $w_i$
  clears that sector's frozen allocation --- which is exactly the
  `problem_sectoral` system of ADR-0020 option C, executed as the `BF` row of
  `matrix_5x3_v6` (\texttt{S1}). The all-rigid corner is therefore anchored on an
  executed, gated and published set of cells, and it is the implementation
  canary (below). \textcolor{revisionV3}{The restricted variant has no such
  anchor, which is the second reason it is not the specification.}
- **The elasticity ladder.** Raising a uniform $\eta_s$ traces the price
  response from its maximum at $\eta_s = 0$ to zero in the limit: as sectoral
  supply becomes perfectly elastic, each real wage is pinned at its anchor and
  prices stop moving. That is a statement about the model, not about one
  parameter choice.
- **Demand-sensitive prices.** \textcolor{revisionV3}{Each supply condition ties
  $w_i$ to quantities through $L^{cm}_i(p, y_i, w_i)$, so the zero-profit block
  acquires an $y$-dependence it does not have in any single-wage closure; the
  free wage is not pinned by technology and the numeraire alone, which was the
  restricted variant's defect.}
- **Determinacy.** The ADR-0017 round-gain criterion applies to the fully pinned
  fixed-wage system and does **not** carry over. The sectoral-wage systems are
  stiff (the $eta = 0$ endpoint's Jacobian condition number is $\sim 9.4\times10^7$
  against $\sim 52.8$ mobile), so the admission question and the residual floor
  must be measured, not asserted.

# Scenarios to evaluate

The specification above fixes a system; this section fixes what is run and what
is reported. The slots cover the general closure.

| Slot | Closure and parameters | Cells | What it demonstrates | What to report |
| --- | --- | --- | --- | --- |
| `S0` | the existing fifteen cells, re-run for provenance | `matrix_5x3-v6-*` | the baseline, and the rows in which prices do not move | the price table as it stands: `max abs(p - 1)` and the deflator |
| `S1` | `eta = 0` (BF): all-sectoral rigid wages at the frozen allocation | `BF-F1`, `BF-F2`, `BF-F3` (already executed and gated) | the largest demand-driven price response in the family, with **no elasticity parameter at all**: the wages solve the first-order conditions at the frozen allocation | `max abs(p - 1)`, the deflator, the wage dispersion, the external position |
| `S2` | general closure, uniform elasticity ladder: `eta_s,i = eta_s` for every sector | `eta_s` in 0.1, 0.5, 1, 2, each across F1, F2, F3 | the central exhibit: the price response as a function of the sectoral supply elasticity, from its `S1` maximum toward zero | the pass-through measure (below), plus employment and welfare at each rung |
| `S3` | general closure, two-group vector: a rigid group at `eta_s,i = 0`, a flexible group at `eta_s = 0.5` | rigid group in {the programme sectors, the largest half by baseline employment}, each across F1, F2, F3 | the rigidity share as a continuous dimension | employment by group, and the price response against the rigid share |
| `S4` | pinned wage vector (Door 1, ADR-0021) | the BF-transplant pin plus two institutional pins, each across F1, F2, F3 | price *heterogeneity* without price *sensitivity*, and the allocation-margin comparison against `S1` | prices, and the allocation difference against the BF cells that share the same wage structure |
| `S5` | optional companion: the utilization externality, `delta` in {+0.5, -0.5} | the programme cells | the opposite price sign (the Kaldor--Verdoorn case at `delta < 0`) | the deflator under both signs, as a bracket on the price response |

## The five slots, measured \textcolor{revisionV4}{\normalsize [added v5]}

\textcolor{revisionV4}{The operator asked for the slots compared. \texttt{experiments/probes/probe14\_s1\_s5\_comparison.jl} measures all five on the full-71 calibration (no \texttt{src/} change). The discriminator is the one the workplan's goal turns on: \textbf{demand sensitivity} means \texttt{max abs(p-1)} differs across the financing columns F1/F2/F3; \textbf{exogenous heterogeneity} means prices move but by the same amount in every column.}

| Slot (as measured) | `max abs(p-1)` F1 / F2 / F3 | employment (F2) | wage max/min | verdict |
| --- | --- | ---: | ---: | --- |
| \textcolor{revisionV4}{\texttt{S1} rigid ($eta = 0$, the executed v6 \texttt{BF} row)} | 0.2214 / 0.2790 / 0.2790 | 1.00000000 | 1.72 | demand-sensitive |
| \textcolor{revisionV4}{\texttt{S2} uniform $\eta_s = 0.5$} | 0.0958 / 0.1057 / 0.1057 | 1.00125507 | 1.26 | demand-sensitive |
| \textcolor{revisionV4}{\texttt{S2} uniform $\eta_s = 2$} | 0.0356 / 0.0372 / 0.0372 | 1.00191380 | 1.09 | demand-sensitive |
| \textcolor{revisionV4}{\texttt{S3} two-group (rigid = programme sectors)} | 0.2151 / 0.2720 / 0.2720 | 0.99662119 | 1.70 | demand-sensitive |
| \textcolor{revisionV4}{\texttt{S4} pinned vector (programme sectors +10\%)} | 0.05782 / 0.05782 / 0.05782 | 1.00686452 | 1.10 | heterogeneity only |
| \textcolor{revisionV4}{\texttt{S5} utilization $\delta = +0.5$} | 0.0948 / 0.1045 / 0.1045 | 1.00000000 | 1.00 | demand-sensitive |
| \textcolor{revisionV4}{\texttt{S5} utilization $\delta = -0.5$} | 0.1361 / 0.1239 / 0.1239 | 1.00000000 | 1.00 | demand-sensitive |

\textcolor{revisionV4}{Three readings, all of which bear on the design choice. (i) \texttt{S4} is the only slot whose price movement is exogenous --- identical to five digits across the financing columns --- so the Door-1/Door-3 distinction is visible in the numbers, not only in the derivation. (ii) The uniform-$\eta_s$ slot \texttt{S2} needs \textbf{no sectoral heterogeneity assumption}: one elasticity for every sector, the same parameter count as BETA, and it already moves prices with demand. The heterogeneity that produces the variety sits in the quantities (N separate labour markets), not in the preferences. The two-group \texttt{S3} is therefore an optional refinement rather than the headline. (iii) The capacity door's reduced form \texttt{S5} reaches the \texttt{S2} magnitude with a single parameter, no wage dispersion and no employment movement at all, and with $\delta < 0$ it reverses the ordering of the financing columns (the Kaldor--Verdoorn sign) --- so it is a genuine alternative route to the same goal, not merely a robustness arm.}

## Controls, not scenarios

These belong in the test suite as much as in the design, and they are the
implementation's canaries:

- **The nesting check.** The general closure at `eta_s,i = 0` for every sector
  must reproduce the `S1` manifests, because those are the same equations. This
  is the canary for the translation: if `S2` at `eta_s = 0` does not equal `S1`,
  the kernel change is wrong. It replaces the ALPHA and GAMMA nesting of the
  restricted specification, which the general one does not have.
- **The neutrality control.** `F2` and `F3` must have identical prices, wages
  and allocations, as in every regime since ADR-0019; a difference beyond solver
  noise means the new closure broke the accounting.
- \textcolor{revisionV3}{\textbf{The level control, sharpened.} A uniform rescale of
  the anchor wages is a change of units only when the whole real wage path is
  rescaled with it; with the CPI numeraire fixed, rescaling the anchors alone is
  a real shock. Probe 12 measured both: the pure pin (no numeraire beyond the
  wage) is invariant to $2.7\times10^{-12}$, while a segmented pin with a
  fixed CPI numeraire moves the allocation by 70\%. The general closure's own
  invariance --- a uniform rescale of $(\bar w, \bar\Pi)$ --- should be stated
  and measured explicitly, not assumed from Corollary 1 of
  \texttt{paper/equivalence.tex}, which holds on the single-wage demand-only
  slice. This replaces the Version-3 level control, which asserted the
  invariance outright.}

## The reporting unit

The paper needs one number per cell, not a price vector. The natural measure is a
**pass-through**: the change in the GDP deflator (or in `max abs(p - 1)`) per
percentage point of programme spending relative to GDP. It is comparable across
labour rows, it is defined in every regime, and for the existing cells it is
computable from the manifests --- zero to machine precision in ALPHA, BETA,
GAMMA and DELTA, and of order 0.5 in the BF row, where the deflator moves 0.63 to
0.72 per cent against a programme of 1.33 per cent of GDP. Every slot above
should be reported with it, so that `S2` and `S3` read as a curve in one
parameter rather than as a table of vectors.

# Reconciliation with the five-door menu of `docs/ETAs.md`

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
  elasticities; Version 3 of this workplan specified the restricted sticky-set
  variant instead, chosen because it nests ALPHA and GAMMA. \textcolor{revisionV3}{That choice is now
  settled by measurement: the restricted variant does not create a demand
  channel to prices (probe 12), whereas the general closure does and nests onto
  the executed $eta = 0$ row, so it is not a choice between two working routes
  but between a working route and a decorative one. \textbf{Decision (v4): implement
  the general closure}, with the scenarios reported at a low-dimensional
  parametrisation (a rigid group at $\eta_{s,i} = 0$ and a flexible group at a
  common elasticity), because $N$ elasticities are weakly identified under
  demand-only shocks.}
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

- **Step 1, the two prototypes (no `src/` change).**
  \textcolor{revisionV3}{\texttt{experiments/probes/probe12\_segmented\_wages.jl}
  (the restricted variant) is done and refutes the sticky-set specification:
  nesting ALPHA/BETA/GAMMA bit-exact, account closed, employment responsive, but
  prices and the free wage demand-free.
  \texttt{experiments/probes/probe13\_sectoral\_labour.jl} is the successor: the
  general closure on the full-71 calibration, measuring the endpoint nesting onto
  the committed \texttt{BF} (v6) manifests, the uniform $\eta_s$ ladder and the
  two-group variant, the demand sensitivity of prices across the three financing
  cells, the identity gap through the price-weighted clearing residual (the proxy
  validated in probe 11, ratio 1.0000), and the Jacobian conditioning and residual
  floor. Falsifiable predictions, recorded before running: prices move with the
  financing cell; the price response falls monotonically as $\eta_s$ rises;
  employment absorbs more as $\eta_s$ rises; the $\eta_s = 0$ corner reproduces
  \texttt{S1} bit-for-bit.}
- **Step 2, ADR.** Draft the sectoral-labour-market ADR (proposed) with the
  specification above, the nesting table, the measured demand sensitivity and the
  admission check. Registry: a new labour-closure entry (working id `BETA-S`, the
  sectoral generalisation of BETA).
- **Step 3, promotion.** Dispatch-only kernel change: a sectoral closure type
  carrying the elasticity vector (and, for the two-group reporting, a rigid
  group), its residual function, the `solve` plumbing, and the vector-wage
  handling in `equilibrium_residuals`, `market_clearing_residuals`,
  `external_balance_canary` and `gdp_components` (all of which already accept a
  vector since ADR-0020). The harness needs the third gate for the sectoral
  regime (the vector labour-market residual) and, as at eta = 0, the scalar
  `wage` metric should stay the wage-bill-weighted aggregate with `wage_min` and
  `wage_max` carrying the structure.
- **Step 4, tests.** The endpoint nesting at `eta_s,i = 0` on the fixtures and on
  one full-71 cell (`S1` reproduction); the level invariance under a uniform
  rescale of the anchors; the demand-sensitivity assertion (the price deviation
  must move with the financing cell, which no existing closure satisfies); the
  identity gate; the determinacy check.
  \textcolor{revisionV4}{Status (same session): Step 2 done --- ADR-0022
  (accepted 2026-09-19). Step 3 done --- the dispatch-only kernel change:
  \texttt{MobileLaborCESElasticities} carries the optional
  \texttt{eta\_s\_vec}, \texttt{SectoralElasticLaborClosure} describes the
  N-market form, \texttt{problem\_sectoral} gains the real-wage supply term with
  the $\eta_{s,i} = 0$ branch kept arithmetically separate for bit-identity, and
  \texttt{solve} and the residual dispatchers select the 3N+1 system. Step 4 done
  --- \texttt{tests/test\_sectoral\_labour.jl}, wired into the suite: the
  $\eta_{s,i} = 0$ nesting canary against the executed \texttt{BF} endpoint, the
  demand-sensitivity assertion, the scalar-BETA contrast, and the identity.
  Step 5 done --- the per-cell \texttt{eta\_s\_rigid\_group} schema, the
  preregistered design and the \texttt{matrix\_5x3\_v9} generation: 33 cells, all
  executed, all gates pass. Two aborted generations are retained as the record of
  the promotion plumbing rather than of the model --- \texttt{v7} (the new
  closure had no diagnostics method, so the eighteen sectoral cells failed at
  gate evaluation after a successful solve) and \texttt{v8} (the ADR-0019
  acceptance gate collapsed the sectoral wage vector to a scalar, so the
  goods-market clearing check used a wrong common wage). Both were caught by
  existing gates, and the corrected batch was opened only after a full-path
  pre-flight over every cell type. The executed numbers reproduce probe14 to six
  digits.}
- **Step 5, design and generation.** A design carrying the existing fifteen cells
  plus the sectoral cells (the `eta_s` ladder and the two-group variant, each
  across F1, F2, F3), preregistered before execution. A `src/` change supersedes
  the current generation's provenance, so the existing cells are re-run in the
  same batch rather than mixed with it.
- **Step 6, records.** Flow table from the new manifests, the assessment block on
  demand-sensitive prices, the registry and log entries, the board, the gate.

# Decision points for the operator

- \textcolor{revisionV3}{\textbf{The elasticity vector.} How the
  $N$ elasticities are pinned. The recommendation is a low-dimensional
  parametrisation --- a rigid group at $\eta_{s,i} = 0$ and a flexible group at a
  common $\eta_s$ --- because $N$ elasticities are weakly identified under
  demand-only shocks. The rule that selects the groups (``the programme
  sectors'', ``the largest half by baseline employment'') should be
  preregistered.}
  \textcolor{revisionV4}{The measured comparison (above) settles the default: the \textbf{uniform}
  $\eta_s$ slot \texttt{S2} is the headline, because it needs one elasticity and no
  statement that any sector is more wage- or leisure-responsive than another;
  the two-group \texttt{S3} becomes the robustness slot. If the operator prefers a
  single-parameter route with no labour-market heterogeneity at all, the
  capacity door's reduced form (\texttt{S5}) reaches the same price magnitude with one
  parameter $\delta$ and no wage dispersion.}
- **The flexible segment's supply elasticity.** `eta_s = 0` makes the family
  nest `S1` at one end; larger values trace the ladder. \textcolor{revisionV3}{The
  calibration choice determines how far the family departs from the executed
  rigid corner and is the single number the ladder reports against.}
- **The closure's id and role.** Whether the sectoral closure is presented as a
  family with `S1` as its all-rigid endpoint (which answers the redundancy
  question) or as a new row. The first is cleaner; the second changes the paper's
  table structure. \textcolor{revisionV3}{(This replaces the Version-3 GAMMA-role
  question, which the restricted variant raised and which the general one does
  not: the general closure's rigid corner is the BF row, not GAMMA.)}
- **Scope.** Whether the sectoral price response becomes a headline dimension of
  the revision or a follow-up paper. Nothing in steps 1 to 3 commits either way;
  the probe is worth doing regardless.

# Risks and how they surface

- **The mirrored prototype.** The eta = 0 work showed that a correct prototype
  can still be mistranslated into the kernel, where the defect is structural and
  silent. The endpoint nesting is the self-diagnosing check: if `S2` at
  `eta_s,i = 0` does not reproduce the `S1` manifests, the translation is wrong.
- **Determinacy.** Measured in step 1, not assumed; the sectoral systems are
  stiff, so the acceptance gate and the polish target must be set from the
  measured residual floor.
- **Provenance.** A new generation is required; the existing cells must be
  re-run in the same batch and reproduce their manifests.
- **Interpretation.** \textcolor{revisionV3}{The elasticity vector is exogenous
  and the sectoral wage dispersion is an outcome. A welfare difference across
  $\eta_s$ is not a statement about wage policy, and the elasticity does not
  identify a labour-supply elasticity from the demand-only design alone --- it
  is a sensitivity ladder over an assumed value.} \textcolor{revisionV3}{The
  restricted variant's measured failure is itself the demonstration that the
  distinction between the two doors is substantive, and it should be cited when
  the four closures are compared.}

# Evidence base

`experiments/probes/probe11_gamma_wage_structure.jl` (level invariance, the
identity proxy and its validation, conditioning);
`probe7_sectoral_wages_eta0.jl`
(the eta = 0 closure the vector-wage apparatus comes from);
\textcolor{revisionV3}{\texttt{probe12\_segmented\_wages.jl} (the restricted
variant: non-nesting of the free wage, acute demand-freeness, the tilted-pin
case, the pin-level behaviour); \texttt{probe13\_sectoral\_labour.jl} (the
general closure: nesting onto \texttt{S1}, the $\eta_s$ ladder, the two-group
variant, demand sensitivity);}
`probe8_promotion_verification.jl` and `probe9_nonbf_reproduction.jl` (the
promotion's verification pattern); `runs/matrix_5x3-v6-*` (the price table
above, and the `S1` cells the general closure nests onto);
`docs/VariationinGamma.md` (the four doors and the claim-by-claim
provenance); ADR-0020 and ADR-0021; the `eta_s` identifiability item in
`docs/DOCS_ASSESSMENT.md`.

# Revision Log

- **Version 1** (September 2026)
- **Version 2** \textcolor{revisionV1}{(September 2026)} --- Added the reconciliation with the five-door menu of `docs/ETAs.md`: the labour door is that document's candidate 1 as well, so the priority stands; the specification choice between the general sectoral-labour-market closure and the restricted sticky-set variant is flagged with a recommendation (general, reported at a low-dimensional parametrisation, because N sectoral elasticities are weakly identified under demand-only shocks); the Verdoorn arm is recorded as the cheap companion with the opposite price sign; and the shared per-cell parameter source is identified as one schema change serving both routes, with the note that the `ADR-0018` id cited in `docs/ETAs.md` is already taken.
- **Version 3** \textcolor{revisionV2}{(September 2026)} --- Added the scenario section (section 4, with its three subsections), which the operator flagged as missing: six slots (`S0` the re-run baseline, `S1` the already-executed `eta = 0` row, `S2` the uniform elasticity ladder, `S3` the two-group rigidity share, `S4` the pinned-vector complement, `S5` the optional Verdoorn companion), three controls that act as implementation canaries (the nesting check at `eta_s,i = 0`, the F2/F3 neutrality control, the level control), an explicit statement of what the general formulation captures that the restricted one does not, and a single reporting unit (a pass-through per percentage point of programme spending) so that the slots read as a curve rather than a table of vectors. The specification bullet of the reconciliation section gained the employment-behaviour difference between the two variants.
- **Version 4** \textcolor{revisionV3}{(September 2026)} --- Step 1 executed (`probe12_segmented_wages.jl`) and the restricted sticky-set specification **refuted**: with a single free wage, zero profit plus the CPI numeraire already determine the prices and that wage, so the flexible labour-market condition never reaches the wage (positions: S = empty nests ALPHA/BETA and S = all nests GAMMA bit-exactly; the account closes to 1e-16; employment responds; but `w_f = 1`, `max abs(p-1)` ~ 1e-14 in every financing cell, and a +10\% tilted pin moves prices 5.5e-2 identically across F1/F2/F3). The title, the specification section and the reconciliation now state the **general sectoral-labour-market closure** (N supply conditions, `eta_s` vector), whose all-rigid corner is the executed `BF` (v6) row, and Step 1 gains `probe13_sectoral_labour.jl` as its successor. Decision points updated (the elasticity-vector rule replaces the sticky-set rule; the closure-role question replaces the GAMMA-role question), the level control is sharpened (Corollary 1 holds on the single-wage slice, not in general), and the evidence base and risks record probe 12.
- **Version 5** \textcolor{revisionV4}{(September 2026)} --- Added the measured comparison of the five scenario slots (`S1`..`S5`) from `experiments/probes/probe14_s1_s5_comparison.jl`: `S4` (the pinned vector) is the only slot whose price movement is exogenous (identical across F1/F2/F3), while `S1`, `S2`, `S3` and `S5` are demand-sensitive; the uniform-`eta_s` slot reaches the goal with a single elasticity and no sectoral-heterogeneity assumption, and the capacity door's reduced form reaches the same magnitude with one parameter and no wage or employment movement. The decision points record the resulting default (uniform `eta_s` as the headline, the two-group as robustness).