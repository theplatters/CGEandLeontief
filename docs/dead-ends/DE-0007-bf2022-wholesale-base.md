# DE-0007 — Adopting Baqaee–Farhi (2022) wholesale as the base model

- **Status:** dead end as a base-model replacement; the complementarity
  mechanism remains valid for `ZETA`
- **Recorded:** 2026-09-17
- **Origin:** closure-design assessment, 2026-09
- **Related:** closure `ZETA`; `roadmaps/vertdict.md`
- **Evidence:** `roadmaps/vertdict.md` ("Do not replace the model wholesale
  with Baqaee–Farhi (2022)"; the 2022 model has sector-specific sticky
  labour, capital, heterogeneous households, nominal rigidities and a
  monetary closure), `bf_replication2/src/network.jl:45` (N sticky labour
  factors rather than one mobile market)

## What was tried

Use the 2022 Keynesian production-network model as the common core, on the
grounds that it handles demand-driven unemployment natively.

## Why it seemed promising

It answers the referee's demand for a demand-constrained labour market, and
the port already exists (`bf_replication2/`).

## How it failed at the design stage

Its native labour structure is **sector-specific sticky labour**, not one
mobile aggregate market; importing it changes the paper's question from
"how does the closure move results between CGE and IO endpoints?" to "how do
sectoral wage rigidity and monetary conditions govern green-stimulus
effects?". Its capital, heterogeneous households, nominal rigidities, and
monetary closure are unnecessary for the static comparative exercise and
would obscure the IO endpoint.

## What replaced it

The 2019-style competitive CES production-network core as the common base,
with the 2022 complementarity machinery adapted only for the unemployment
closure `ZETA` (`ROADMAP.md` §5 Closure C; `roadmaps/vertdict.md`
recommendation).

## Revival conditions

Switch entirely to the 2022 base only if the paper's main question changes to
sectoral wage rigidity and monetary conditions — a separate ADR is required
before any such move.
