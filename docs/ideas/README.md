---
title: "Ideas register: unexplored routes recorded before they are lost"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-09-18"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [ideas, register, unexplored, adr-annexe, process]
last-updated: October 2026
---

**Version 3** (October 2026)
**Version 2** (September 2026)
**Version 1** (September 2026)

# What this register is for

A third register next to `docs/decisions/` and `docs/dead-ends/`, for the middle
case: routes that are worth remembering, that nobody has tested, and that would
be lost in a session log. An idea is not a decision: it carries no weight in the
paper and no entry in `registry/`. Promoting one follows the normal path --- an
ADR, a registry entry, a preregistered design --- and dropping one writes a
`DE-` record instead, so that the reason is not re-litigated later.

Status vocabulary for ideas, and nothing else:

| Status | Meaning |
| --- | --- |
| `unexplored` | Written down, never tested. The default state. |
| `scoped` | A specification exists (variables, equations, cost) but no measurement. |
| `promoted` | Adopted: an ADR and a registry entry exist. The idea note then points to them. |
| `dropped` | Abandoned: a `docs/dead-ends/DE-*` record exists. |

Every note is one file, named `IDEA-NNNN-short-slug.md`, numbered in the order of
writing. Do not renumber, do not reuse, do not delete a note whose idea was
dropped --- point it at its dead end.

# Index

| Id | Idea | Status | Anchored in |
| --- | --- | --- | --- |
| IDEA-0001 | Productivity shocks framed as climate change | unexplored | Door 2 of `docs/VariationinGamma.md`; the `eta_s` identifiability item in `docs/DOCS_ASSESSMENT.md` |
| IDEA-0002 | The five doors to demand-sensitive prices (labour, markup, capacity, external, technology), with the Verdoorn arm | promoted | Section 5 of `docs/ETAs.md`; Lemma 1 and Corollaries 1 and 3 of `paper/equivalence.tex`; `docs/ETAs.md` and `docs/VariationinGamma.md`; ROADMAP Phase 7; `experiments/probes/probe20_pricing_door.jl`; ADR-0029, ADR-0030, `registry/closures.toml` (P1, P2) |
| IDEA-0003 | Sobol decomposition of the sectoral labour response (grouped over the elasticity parametrisation) | scoped | `docs/DOCS_ASSESSMENT.md` Stage 2 item 2; `docs/archive/definitive_guide.md` Phase 7 (archived by ADR-0026); `docs/grouping_rule_evidence.md`; `OpenCGE@ADR-0022` (kernel `DECISIONS.md`, migrated by ADR-0028) |
| IDEA-0004 | Sectoral effects without a Sobol (contribution, incidence, network, similarity, wage-to-price) plus the case for reviving the submission's sectoral figures | scoped | `paper/submission-metro/submission.tex` §5-§6; `OpenCGE@ADR-0022` (kernel `DECISIONS.md`, migrated by ADR-0028); `docs/ideas/IDEA-0003-sectoral-sobol.md` |

| IDEA-0005 | Pricing-door workstream: handover of what is built, decided, measured and next (capacity/P1, markup/P2, Stage A retention, Stage B distributed profit) | scoped | `docs/decisions/ADR-0029-pricing-door.md`; `docs/decisions/ADR-0030-distributed-profit-income.md`; `registry/closures.toml` `[pricing.P1]`/`[pricing.P2]`; `ROADMAP.md` Phase 7; OpenCGE v0.2.0 (`e119f8d7`) |

All five notes stay in this register: each affects the paper's methods or
results, not only the kernel. The kernel-side ADR/DE records they cite were
migrated to OpenCGE's `DECISIONS.md` in the repo split (ADR-0028); their ids
are unchanged (`OpenCGE@ADR-…`, `OpenCGE@DE-…` in `registry/`).

# Notes on the process

The register is deliberately separate from `ROADMAP.md` and the ADRs: those carry
commitments, this carries candidates. It is also deliberately separate from
`docs/log/`, which records what happened. An entry here has not happened.

When an idea is promoted, its note keeps its number and gains a pointer to the
ADR and the registry entry, so that the paper trail runs from the first mention
to the executed run.

# Revision Log

- **Version 1** (September 2026)
- **Version 2** (September 2026) --- registered IDEA-0004 (sectoral effects without a Sobol, plus the revival verdict on the submission's sectoral figures); IDEA-0003 rewritten as a vantage point and raised from `unexplored` to `scoped`, with `docs/decisions/ADR-0022-sectoral-labour-markets.md` added to its anchors.
- **Version 4** (October 2026) --- IDEA-0002 promoted for the capacity and markup doors (ADR-0029, ADR-0030; closures P1/P2); the external, fixed-factor and technology doors remain ideas.
- **Version 3** (October 2026) --- IDEA-0002 raised from `unexplored` to `scoped` (ROADMAP Phase 7 specifies the markup and capacity doors; probe20 measures them), anchors extended.
- **Version 4** (October 2026) --- registered IDEA-0005 (the pricing-door handover: state, gaps, sequence and traps, written for a successor). It records that IDEA-0002 now meets this register's own `promoted` criterion --- `ADR-0029`/`ADR-0030` accepted and `[pricing.P1]`/`[pricing.P2]` registered --- and should be raised with pointers; left to the operator.
