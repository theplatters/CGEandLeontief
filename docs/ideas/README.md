---
title: "Ideas register: unexplored routes recorded before they are lost"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-09-18"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [ideas, register, unexplored, adr-annexe, process]
last-updated: October 2026
---

**Version 8** (October 2026)
**Version 7** (October 2026)
**Version 6** (October 2026)
**Version 5** (October 2026)
**Version 4** (October 2026)
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
| IDEA-0006 | The Verdoorn grid: which matrix closures to carry the capacity door (the dissolution of the three Lemma-1 equivalences, and why BETA rather than DELTA or BF) | scoped | `runs/pricing_*` (69 executed cells: GAMMA 45, ALPHA 24, BETA/BF/DELTA 0); `registry/closures.toml` `[pricing.P1]` (the dissolution claim); `experiments/designs/pricing_capacity_v2.toml`; Appendix C (prop. `gd`, `ab`, `nest`); Section 6.4 of the revision draft |
| IDEA-0007 | Restructuring the remaining appendices: four appendices by function (data, robustness, pricing door, traditions) instead of six by topic | scoped | The measured appendix block of 2026-10-18 (`revision/sections/`: 12,851 words, 59 percent of them in Appendix B; citation map counted from Sections 1-8); the A/B/C/H merge of 2026-10-18; Section 4's two accounting pointers into Appendix C; Section 6's five forwards into Appendix D |
| IDEA-0008 | One index convention for the kernel: `i` the sector the equation is about, `u` an input, `j` a user --- the rule, its scope, the rejected alternative and the six sites that still violate it | scoped | `revision/sections/(B)model-closures-proofs.tex` lines 315, 485, 493, 495; `revision/sections/(4)kernel-closures-financings.tex` (the $\Omega$ definition and the calibration sentence); `revision/sections/(H)pricing-door.tex` (aligned 2026-10-23) |
| IDEA-0009 | Pattern versus level in the sectoral response: the IO structure fixes the pattern, the labour closure the level | scoped | `runs/matrix_5x3-v7\|v8\|v9-*` (fifteen `solution.csv`, identical to the reported v11 runs on every headline metric); letter R2.0; Section 5.2 |


All eight research notes stay in this register: each affects the paper's methods
or results, not only the kernel. IDEA-0008 is the exception in kind --- a
presentation convention for the kernel's notation, registered because it changed
the appendix text and because the letter table must outlive the session that
derived it. The kernel-side ADR/DE records they cite were
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

- **Version 7** (October 2026) --- registered IDEA-0008 and its revised scope: local index roles, generic-index exceptions, baseline calibration shares and the proposed manuscript migration.
- **Version 6** (October 2026) --- registered IDEA-0007 (restructuring the remaining appendices: four by function instead of six by topic, with the measured block arithmetic, the three mismatches, the four moves and the F-placement option). Records adjacent findings: the appendix equation and table counters carry a fixed `A` prefix irrespective of the appendix; the register's own anchors still name the pre-merge appendix letters; and the version block at the top of this README had stopped at Version 3 while the log ran to Version 5, now brought level.
- **Version 5** (October 2026) --- registered IDEA-0006 (the Verdoorn grid: which matrix closures carry the capacity door; the two-gauge scope measured from `runs/`, BETA as the cell that retires Proposition `ab`'s dissolution, DELTA thin because its restriction lives on the substitution axis, BF optional and half-computational). Records the adjacent finding that Appendix J's "42 priced cells" is not reproducible from the manifests.
- **Version 1** (September 2026)
- **Version 2** (September 2026) --- registered IDEA-0004 (sectoral effects without a Sobol, plus the revival verdict on the submission's sectoral figures); IDEA-0003 rewritten as a vantage point and raised from `unexplored` to `scoped`, with `docs/decisions/ADR-0022-sectoral-labour-markets.md` added to its anchors.
- **Version 4** (October 2026) --- IDEA-0002 promoted for the capacity and markup doors (ADR-0029, ADR-0030; closures P1/P2); the external, fixed-factor and technology doors remain ideas.
- **Version 3** (October 2026) --- IDEA-0002 raised from `unexplored` to `scoped` (ROADMAP Phase 7 specifies the markup and capacity doors; probe20 measures them), anchors extended.
- **Version 4** (October 2026) --- registered IDEA-0005 (the pricing-door handover: state, gaps, sequence and traps, written for a successor). It records that IDEA-0002 now meets this register's own `promoted` criterion --- `ADR-0029`/`ADR-0030` accepted and `[pricing.P1]`/`[pricing.P2]` registered --- and should be raised with pointers; left to the operator.
- **Version 7** (October 2026) --- registered IDEA-0008 (one index convention for the kernel: `i` the sector the equation is about, `u` an input, `j` a user). Records the measured drift between the price index and the intermediate-demand equation, the alternative that was rejected and why, the six remaining violations with their line numbers, and the notation fixes made the same day (the transposed $\Omega$ in Section 4, Appendix H's letter pair, the dead `A^bill` symbol, and the scaling factor in the program bundle).
- **Version 8** (October 2026) --- registered IDEA-0009 (pattern versus level in the sectoral response). Records the
  first measurement of the sectoral solutions that the reported runs ship: the financing regime scales the sectoral
  response and leaves its ordering intact (rank correlations 0.97 to 1.00), the labour closure fixes the level and
  hence the sign of the typical sector (54 of 71 flip, none the other way), so the ordering breaks only where the
  flexible-price closures let the expansion run. Basis: `runs/matrix_5x3-v7|v8|v9-*`, identical to the reported v11
  runs on every headline metric. Written the day the letter's R2.0 dropped its surviving-closure claim.
