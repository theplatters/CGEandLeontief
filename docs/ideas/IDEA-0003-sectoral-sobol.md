---
title: "IDEA-0003: Sobol decomposition of the sectoral labour response"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-09-20"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [ideas, sobol, sensitivity, sectoral, labour, variance-decomposition, closure-taxonomy]
last-updated: September 2026
---

**Version 2** (September 2026)
**Version 1** (September 2026)

Status: `scoped` (the factor space, the block structure and the cost are specified; nothing measured).

# What the original idea was

Version 1 proposed to variance-decompose the **sectoral** response of the ADR-0022 family --- the price response `max |p-1|`, wage dispersion, employment, the consumption/welfare index --- over the elasticity parametrisation, with Sobol first-order and total indices. The question: *which sectors' wage-responsiveness drives the headline, and how much of the outcome sits in interactions?*

It fixed three scope constraints and they all still stand.

- **Not 71-dimensional.** A full Sobol over `eta_s,i` (`i = 1..71`) is wrong: stiff 3N+1 solves are the cost unit and the individual elasticities are weakly identified.
- **A preregistered measure.** Sobol indices presuppose a probability measure on the inputs; the ranges come from `docs/grouping_rule_evidence.md` and are preregistered like any other design cell.
- **Monotone responses expected.** First-order indices should dominate; reporting ST close to Sf is then an honest near-separability result, not an embarrassment.

It named two routes: (a) a grouped Sobol over K ~ 5-8 sector blocks x financing, Saltelli sampling, first and total order (the publishable headline); (b) Morris screening first, for the cheap 71-sector map of which sectors matter. Both routes stay open; what changes below is the *factor space they run over*.

# What version 2 adds

Three findings from the 2026-09-20 review of the executed work change what a sectoral Sobol means here. None of them weakens the idea; they relocate it.

- **The aggregate Sobol is done and cited.** Factors `(eta, epsilon, theta, sigma)`, `eta` in {0, 1} only, output aggregate real GDP: `theta` 0.3951 > `sigma` 0.2734 > `epsilon` 0.1650 > `eta` 0.1571, interaction share `1 - sum(Sf)` = 0.0093, 32/32 grid points, 0 failures. The sectoral version was never run.
- **The sectoral dimension was promoted from garnish to mechanism** by ADR-0022. Under Lemma 1 of `paper/equivalence.tex` every single-wage closure is demand-free in prices (`max |p-1|` of order 1e-14); the only price-moving rows in the kernel are the wage-vector rows. "What happens sectorally" is therefore the channel through which prices respond *at all*, not a downstream detail --- which is what makes the sectoral decomposition load-bearing rather than decorative.
- **The closure axis is categorical, and its levels are nested.** Version 1's sketch silently assumed the old continuous factor set could be re-pointed at the new closures. It cannot, for four reasons set out below. The honest object is a *segmented* design, not one grand factorial.

# The aggregate Sobol as built (the reference point)

Everything below refers to this machinery, which already exists in the kernel.

| Property | As built |
| --- | --- |
| Factors | `eta`, `epsilon`, `theta`, `sigma` |
| `eta` levels | {0, 1} only --- the ADR-0010 endpoints; intermediates retired |
| Design | full factorial, balanced; ANOVA sum-of-squares |
| Output | aggregate real GDP (`nominal_gdp`, `:sectoral_q`, `:sectoral_p` exist in code, never run) |
| Indices | `Sf = SSf/SStot`, `STf = 1 - SS_-f/SStot`, interaction `STf - Sf`, `1 - sum(Sf)` |
| Failure handling | failed grid points retained as NaN and counted in `n_failed`; DE-0006 forbids renormalising over main effects |
| Result | `theta` 0.3951, `sigma` 0.2734, `epsilon` 0.1650, `eta` 0.1571; interaction 0.0093 |

# What a sectoral Sobol entails

## Why the version-1 factor set cannot be reused unchanged

- **`eta` is no longer a factor --- it is a closure level.** BF is the `eta = 0` endpoint, ALPHA the `eta = 1` endpoint. Carrying both "closure" and "`eta`" is collinear by construction; the old `eta` slot *becomes* the labour-closure factor.
- **The closure factor is categorical with five levels, and the levels are partly nested.** BETA in its scalar form is identical to ALPHA under demand-only shocks (the real-wage anchor binds, so the supply curve returns `Lbar` for any elasticity); DELTA is GAMMA plus the Leontief limit and reproduces GAMMA to six digits; GAMMA's F2 and F3 are real-neutral. A categorical factor is a perfectly legal ANOVA factor, but its main effect is a contrast whose levels partly coincide, and some apparent interactions are definitional nesting rather than economic interaction. That must be stated in the table, not smoothed.
- **The new inputs are vectors, not scalars.** `eta_s,i` (71 entries) and the GAMMA pin `wbar_i` (71 entries). A scalar Sobol cannot host them; they must be reduced by grouping (K ~ 5-8 blocks), by the rigid-group share, or by a tilt magnitude.
- **Structural zeros break the balanced factorial.** `eta_s` only exists inside the BETA sectoral row; the tilt only inside GAMMA/DELTA. A single grand design has impossible cells, and DE-0006 forbids silently completing an incomplete design.

## The segmented design

Three blocks, each balanced within itself. The closure map is a *designed contrast*; only blocks B and C are sensitivity analyses in the classic sense.

| Block | Factors | Levels | Output | Scope |
| --- | --- | --- | --- | --- |
| A. Closure map | labour (cat.) x financing (cat.) x `theta` x `epsilon` x `sigma` | 5 x 3 x 3 x 3 x 3 = 405 | aggregate + sectoral vectors | whole matrix |
| B. Sectoral elasticity | `eta_s` structure x `theta` x `epsilon` x `sigma` x financing | 6 x 3 x 3 x 3 x 3 = 486 | prices, employment, allocation | BETA sectoral row only |
| C. GAMMA structure | tilt pattern (cat.) x tilt magnitude x `theta` x `epsilon` x `sigma` x financing | 3 x 3 x 3 x 3 x 3 x 3 = 729 | allocation, employment, real income | GAMMA/DELTA rows only |

Those counts are the full-resolution version. A reduced core at two elasticity levels (`theta`, `epsilon`, `sigma` in a corner set, 8 points) gives 120 / 144 / 216 solves respectively, which is near version 1's "~100 cells" sketch. The sizing decision belongs to the operator and to the preregistration, not to this note.

## The output choice decides whether the design is informative

This is the sharpest lesson of the review and it is easy to get wrong.

- For a **price** output, four of the five labour-closure levels are identical (Lemma 1), so the labour factor collapses onto a near-binary "does a wage vector exist" indicator and its `Sf` is uninformative by construction. The substantive content sits in blocks B and C, not in block A.
- For an **employment or allocation** output the labour factor does carry genuine variation, and block A is the informative decomposition.

A single decomposition must therefore not be run and then read as a general statement. State the output in the design file and keep the two readings visibly distinct in any table.

## The GAMMA sub-variation

ADR-0021 makes the pinned wage **vector** a design parameter of the fixed-wage closure, and it has a structural zero built in:

- A **uniform rescale** of the pin is a numeraire change (measured invariant to 8.9e-16). It is not an instrument and must be excluded, or it enters as a factor with `Sf = 0` by theorem --- which is worse than useless, since it dilutes the reported shares.
- Only the **relative** structure is an instrument: the tilt *pattern* (categorical: programme-aligned push, programme-aligned restraint, broad tilt) times the tilt *magnitude* (scalar). That is block C.
- `eta_s -> infinity` is GAMMA-like **in prices only, and only asymptotically**: the direct 3N+1 formulation becomes singular (`eta_s = 1e6` fails to converge). GAMMA is therefore not a clean sixth level of the `eta_s` factor at finite values; it is the limit the family approaches.

## Measure and preregistration

A categorical factor needs declared design weights, not an implicit uniform measure. Uniform over matrix cells is the natural choice, but it is a choice and it must be written down, together with the factor levels and the assumed measure --- the pending `SobolResult` extension item in `ROADMAP.md`. Same discipline as any other design cell (ADR-0006).

# Machinery, cost, and what is missing

| Item | State |
| --- | --- |
| `SobolResult`, `variance_decomposition`, `summary_table` | exist in `src/core/diagnostics.jl` |
| `:sectoral_q` / `:sectoral_p` outputs | exist in the function signature, never exercised |
| Per-block factor design | **missing** --- the function takes scalar factor grids only |
| Runner | **missing** --- `experiments/run.jl` has no Sobol design type |
| Sectoral wage and employment vectors | **missing from the artifacts**: `solution.csv` carries `sector,price,quantity` only |
| `src/` change | none required --- no provenance invalidation of the kernel |

Cost is set by the stiff `eta = 0` cells (Jacobian condition number ~9.4e7 against ~52.8 for the mobile all-N system). The v9 generation was 33 cells; block A at full resolution is ~12 times that. The reduced core is the tractable entry point.

# Open questions the design would settle

- Does the reallocation friction matter for sectoral allocation even where it is aggregate-second-order? (`DOCS_ASSESSMENT.md` Stage 2 item 2; `definitive_guide.md` Phase 7.)
- The recombination of the two labour margins --- allocation x supply elasticity --- is recorded in the assessment as a deliberate omission. The sectoral family *is* that recombination, so a total-effect index over the rigid-group factor gives the omission a measurable end condition.
- Which sector blocks' wage-responsiveness drives the price headline? This is block B and it is the version-1 headline question, now well posed.
- Does the rigid-group **rule** dominate the rigid-group **share**? The executed ladder flips the employment sign across the two rules (programme sectors versus largest half by employment), so this is a sign question, not a magnitude question.

# Promotion path

ADR (the factor space and the design weights) + per-block factor design + preregistration (ADR-0006) + one design + a paper table citing the run ids (ADR-0004). An idea note carries no weight until then. The companion note `IDEA-0004` records the sectoral analyses that need **no** Sobol at all and should be executed first --- they are cheaper, they answer the descriptive half of the same question, and their output tells the operator whether block B is worth its 144 solves.

# Anchored in

- `docs/DOCS_ASSESSMENT.md` Stage 2 item 2; `docs/definitive_guide.md` §7.6 and Phase 7; `ROADMAP.md` (the `SobolResult` extension item).
- `docs/dead-ends/DE-0006` (renormalised shares --- do not repeat); `docs/dead-ends/DE-0010` (the retired cbase2 solver ladder).
- `docs/decisions/ADR-0022-sectoral-labour-markets.md` (the sectoral family, the exact nesting, the `eta_s -> inf` limit); `docs/decisions/ADR-0021-wage-structure-fixed-wage-closure.md` (the pinned vector, the numeraire invariance).
- `paper/equivalence.tex` Lemma 1 and the sectoral section (why single-wage closures are demand-free in prices).
- `docs/grouping_rule_evidence.md` (the input measure); `docs/WORKPLAN_SENSITIVE_PRICES.md` v7 decision points (the grouping rule and its sign flip).

# Revision Log

- **Version 1** (September 2026) --- the idea, the three scope constraints, the two routes.
- **Version 2** (September 2026) --- rewritten as a vantage point after the 2026-09-20 review: the aggregate Sobol recorded as built; the four structural reasons the version-1 factor set cannot be re-pointed at the five closures; the three-block segmented design with counts; the output-choice caveat; the GAMMA sub-variation with its numeraire zero-response control; the artifact gap (`solution.csv` carries no wage or employment vector); status raised from `unexplored` to `scoped`. Version 1's constraints and routes are retained, not superseded.
