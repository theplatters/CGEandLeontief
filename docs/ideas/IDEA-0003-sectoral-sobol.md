---
title: "IDEA-0003: Sobol decomposition of the sectoral labour response"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-09-20"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [ideas, sobol, sensitivity, sectoral, labour, variance-decomposition, follow-up]
last-updated: September 2026
---

**Version 1** (September 2026)

Status: `unexplored` (a design sketch exists; nothing measured).

# The idea

Variance-decompose the **sectoral** response of the ADR-0022 family — price
response `max |p-1|`, wage dispersion, employment, the consumption/welfare
index — over the elasticity parametrisation, using Sobol first-order and
total indices. Open question it answers: *which sectors' wage-responsiveness
actually drives the headline, and how much of the outcome sits in
interactions?*

The aggregate Sobol over `(theta, epsilon, sigma, eta)` is done (eta
first-order share ~15.7 %, ST ~ 0.166, little interaction; DE-0006 warns
against renormalising). The sectoral version was never run — the
`definitive_guide.md` Phase-7 item ("Sobol on sectoral quantities — does eta
matter for sectoral allocation?") and `docs/DOCS_ASSESSMENT.md` Stage 2 item
2 ("does the reallocation friction matter for sectoral allocation even where
it is aggregate-second-order?") are both parked asks.

Why it is not a garnish:

- The executed ladder is monotone and the magnitude is dominated by the
  rigid-group *rule* (which flips the employment sign); total-effect
  indices over sector groups would measure exactly the interaction the
  assessment records as a deliberate omission ("recombination of allocation
  x supply elasticity") — the sectoral family *is* that recombination, so
  the omission gains a measurable end condition.
- Its input measure should be the *evidence ranges* collected for the
  grouping rule (`docs/grouping_rule_evidence.md`), which is why it
  sequences after that work, not before it.

# Scope constraints (the design sketch)

- **Not 71-dimensional.** Full Sobol over `eta_s,i (i = 1..71)` is wrong:
  stiff 3N+1 solves are the cost unit, and the individual elasticities are
  weakly identified. Two workable routes: (a) grouped Sobol over K ~ 5-8
  sector blocks (programme incidence, employment size, import exposure,
  wage share ...) x financing, Saltelli sampling, first + total order; (b)
  Morris screening first for the cheap 71-sector map of "which sectors
  matter". Route (a) is the publishable headline.
- **Preregistered measure.** Sobol indices presuppose a probability
  measure on the inputs (ROADMAP's pending `SobolResult` extension item:
  factor levels, probability weights, assumed measure). The ranges come
  from the evidence doc, preregistered like any other design cell.
- **Monotone responses expected.** First-order indices should dominate;
  reporting ST ~ Sf is then an honest near-separability result, not an
  embarrassment.

# Machinery

`SobolResult` / `variance_decomposition` in `src/core/diagnostics.jl`,
`tests/test_variance_decomposition.jl`, the cbase2 Sobol grid runner
convention (32/32 points, 0 failures) — all existing; the work is a design
and a runner, not new infrastructure. `src/` untouched, so no provenance
invalidation.

# Anchored in

- `docs/DOCS_ASSESSMENT.md` Stage 2 item 2; `docs/definitive_guide.md`
  §7.6 and Phase 7; ROADMAP.md (the `SobolResult` extension item);
  `docs/dead-ends/DE-0006` (renormalised shares — do not repeat);
  `docs/WORKPLAN_SENSITIVE_PRICES.md` v7 decision points (the grouping rule
  and its sign flip); `docs/grouping_rule_evidence.md` (the input measure).

# Promotion path (if adopted)

ADR + per-block factor design + preregistration (ADR-0006) + one design
(~100 cells) + a paper table citing the run ids (ADR-0004). An idea note
carries no weight until then.