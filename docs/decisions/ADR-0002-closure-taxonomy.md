# ADR-0002 — Closure taxonomy and stable IDs

- **Status:** accepted
- **Date:** 2026-09-17
- **Related:** ADR-0001, `registry/closures.toml`, `docs/DOCS_ASSESSMENT.md`,
  `cbase2/review.md`

## Context

The word "closure" is overloaded in this project. It denotes at least three
different things:

- the **labour-market closure** (who adjusts: quantity, wage, or both);
- the **financing closure** for the investment shock (tax, external, or
  preference shift);
- the geometric **reallocation parameter** `η` of the BF interpolation, which
  is a market friction, not a labour-supply elasticity and not a closure in
  the CGE sense.

Documents have used the same symbols for different objects (`η` vs `η_s`),
variously labelled GAMMA as "the Keynesian closure" while implementing a
two-sided uncapped real-wage peg, and treated DELTA as a matrix row rather
than a corner. The external review (`cbase2/review.md`) documents these
mislabels. Without stable IDs, scenarios, runs, and paper text cannot refer
to the same objects.

## Decision

Use the following stable IDs. They are **project labels**, not a claim of
equivalence to any canonical closure in the literature; the qualifications
are recorded per entry in `registry/closures.toml`.

**Labour axis (`[labor.<ID>]`)**

| ID | Object | Not to be confused with |
| --- | --- | --- |
| `BF` | Geometric reallocation friction `L_i = L_fixed_i^(1-η)·L_costmin_i^η`, `Σ L_i = L̄` | a supply elasticity; BF 2019's immobile benchmark (`η=0` keeps one common wage) |
| `ALPHA` | Full-employment mobile labour, one wage, `Σ L_i = L̄` | an independent mechanism: it is the `η=1` limit of BF by construction |
| `BETA` | Elastic total labour supply `Σ L_i = L̄·((w/P)/(w₀/P₀))^{η_s}` | an earned labour–leisure model (no income effect) |
| `GAMMA` | Fixed real wage `w/P = w̄`, employment endogenous and uncapped | the canonical one-sided wage floor with rationed employment (that is `ZETA`) |
| `DELTA` | IO-type endpoint: GAMMA plus the Leontief limit (θ, ϵ, σ → 0⁺) at `η=1` | an independent matrix row; it is a corner, and it is a Type II-style multiplier rather than `(I−A)⁻¹` |
| `ZETA` | Unemployment complementarity `0 ≤ L̄−L ⊥ w/P−ω̄ ≥ 0` (future) | a copy of the BF 2022 sectoral sticky-labour structure |

**Financing axis (`[financing.<ID>]`)**

| ID | Object | Note |
| --- | --- | --- |
| `F1` | Preference reallocation: renormalized CES weights, budget-neutral | demand composition, not an autonomous investment multiplier |
| `F2` | Tax-financed public investment `Σ p_i g_i = T(p)` | balanced-budget counterparty required |
| `F3` | Externally/debt-financed programme `Σ p_i g_i = F` | accounting open-economy closure, not a monetary mechanism |

**Symbol discipline**

- `η` is always the BF reallocation parameter; `η_s` is always the BETA
  supply elasticity. They never appear unqualified as "the elasticity".
- Scenario and run ids use the closure IDs (`matrix_5x3-BETA-F2`, ADR-0004).
- DELTA is never listed as an independent row of the evaluation matrix.

## Consequences

- `registry/scenarios.csv`, `experiments/designs/*.toml`, run manifests, and
  the manuscript can refer to the same objects unambiguously.
- Contested mappings stay visible as registry `notes`/`open_gates` instead of
  being silently relabelled. Correcting a mapping is an ADR that supersedes
  this one, not a drive-by rename.
- Adding a closure means adding a new ID here (or a documented extension),
  never reusing an existing ID for a different formulation.

## Enforcement

- `registry/closures.toml` entries must use these ids and include the exact
  formulation string.
- `scripts/status.jl` renders the registry; unknown ids in scenarios.csv are
  warnings on `docs/status.md`.

## Amendment 2026-09-17 — DELTA in the matrix

The 5×3 evaluation matrix (`docs/DOCS_ASSESSMENT.md`) legitimately contains
DELTA as its IO-endpoint row. The statements above that DELTA is "never
listed as an independent row of the evaluation matrix" (Decision bullet and
taxonomy table, "not to be confused with …") are to be read as "**not an
independent mechanism**": DELTA's cells are reached as GAMMA plus Leontief
technology, and the paper reports the limit relations rather than hiding
them as redundancy. The original wording is kept per the append-only rule;
this amendment is authoritative where the two disagree.
