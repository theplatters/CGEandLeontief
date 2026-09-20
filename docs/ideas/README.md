---
title: "Ideas register: unexplored routes recorded before they are lost"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-09-18"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [ideas, register, unexplored, adr-annexe, process]
last-updated: September 2026
---

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
| IDEA-0002 | The five doors to demand-sensitive prices (labour, markup, capacity, external, technology), with the Verdoorn arm | unexplored | Section 5 of `docs/ETAs.md`; Lemma 1 and Corollaries 1 and 3 of `paper/equivalence.tex`; `docs/WORKPLAN_SENSITIVE_PRICES.md` |
| IDEA-0003 | Sobol decomposition of the sectoral labour response (grouped over the elasticity parametrisation) | scoped | `docs/DOCS_ASSESSMENT.md` Stage 2 item 2; `docs/definitive_guide.md` Phase 7; `docs/grouping_rule_evidence.md`; `docs/decisions/ADR-0022-sectoral-labour-markets.md` |
| IDEA-0004 | Sectoral effects without a Sobol (contribution, incidence, network, similarity, wage-to-price) plus the case for reviving the submission's sectoral figures | scoped | `paper/submission-metro/submission.tex` §5-§6; `docs/decisions/ADR-0022-sectoral-labour-markets.md`; `docs/ideas/IDEA-0003-sectoral-sobol.md` |

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