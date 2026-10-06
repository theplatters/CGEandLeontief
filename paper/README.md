---
title: "paper/ --- notes index and disposition"
project: "BFRep/(3)BeyondHulten"
created: 2026-10-04
updated: 2026-10-04
tags: [paper, notes, index, disposition]
---

Index of the paper-planning material after the 2026-10-04 tidy. The live
manuscript skeleton is in `revision/`; this folder keeps the provenance and
the still-relevant source notes.

# What lives here

| Path | What it is | Status |
| --- | --- | --- |
| `narrative_outline.md` (+ `.pdf`) | The master planning note (v3) --- the spine of the revision. | live |
| `submission-metro/` | **The validated original submission** (`submission.tex`, `submission-anon.tex`, `BandF.bib`, `pictures/`). The reference for the response letter and the baseline for "what is novel". | live, untouched |
| `to_evaluate/` | Manuscript-facing source notes (see below). | review pending |
| `tables/` | **Empty** --- the four current flow tables (`matrix_5x3_v10_flows.md` + the three `supply_etas*_flows.md`) moved to `revision/tables/`, the five older generations (`matrix_5x3_v{3,4,5,6,9}_flows.md`) to `paper/superseded/`, on 2026-10-04. | moved |
| `pictures/` | Figure pool + `archive/` (10 older renders). The 10 files that duplicated `submission-metro/pictures/` were removed on 2026-10-04. | live |
| `presentation/` | Beamer decks --- a separate deliverable. | **FROZEN** (`FROZEN.md`) |
| `superseded/` | Superseded drafts and old-generation flow tables. | **FROZEN** (`FROZEN.md`) |
| `README.md` | This index. | --- |

# `to_evaluate/`

| File | Disposition |
| --- | --- |
| `equivalence.tex` | **Ported** into `revision/sections/(A)proofs.tex` (App. A). |
| `closures_and_demand_shocks.md` | **Ported** into `revision/sections/(E)closures-mapping.tex` (App. E). |
| `framing_gamma.tex` (+ `.pdf`) | To port into Sections 4 and 6. |
| `construction_capacity.tex` (+ `.pdf`) | To port into Sections 5 and 6 (the rigid-group evidence base). |
| `possible_plots.md` | Candidate figures for the one heavy-lifting figure. |

# Removed 2026-10-04

| Path | Why |
| --- | --- |
| `paper/main.tex` | A near-duplicate of the validated submission: 730 of 732 lines identical, ratio 0.988. The three differences are a longer abstract variant (preserved in `revision/old-draft-notes.tex`, section `%% Longer abstract`), a missing backmatter block, and one paragraph merge --- no unique content. |
| `paper/BandF.bib` | Same 63 keys as `submission-metro/BandF.bib`; byte-identical to `revision/references.bib`. Its only consumer was `paper/main.tex`. |
| `paper/superseded/` contents | Moved out of the active folder on 2026-10-04 (cobbdouglas.tex, notes-rafi.tex, abstract.typ, matrix_5x3_v{3,4,5,6,9}_flows.md). |
| `revised_manuscript/` | Dissolved: `main.tex` + `chapters/` were an **early adaptation draft** (commit `2efbeb9`/`41144dd`), not a split of the submission; its nine novel passages and both title suggestions are preserved in `revision/old-draft-notes.tex` (six of the passages are also embedded as `CARRIED OVER FROM AN OLD DRAFT` comments in `revision/sections/(2)`,`(4)`,`(5)`,`(6)`, and both titles also sit above `\papertitle` in `revision/main.tex`). Pictures went to `paper/pictures/` (+ `archive/`). |

# Still open

- `pictures/`: 63 files + `archive/` (10). The 10 same-named files are
  byte-identical to `submission-metro/pictures/`; the 10 in `archive/` are
  different renders. Current figures = the 15 `panel_5x3_*`.
- `revision/tables/`: `supply_etas_flows.md` (100 ln) appears superseded by the
  `supply_etas_{prog,unif}_flows.md` pair (85 ln each).
- `to_evaluate/`: three notes still to port (see above).