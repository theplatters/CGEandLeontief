---
title: "Revision draft --- BeyondHulten / Metroeconomica"
project: "BFRep/(3)BeyondHulten"
created: 2026-10-04
status: "scaffold + ported appendices A and E; compiles in-container"
tags: [revision, manuscript, latex, reviewer-response]
---

Revision workspace for the Metroeconomica resubmission. Built from the
universal template `git/ai/RA/publishing/TeXTemplate/`; the section skeleton
follows `paper/narrative_outline.md` (v3); the response skeleton follows the
referee reports in `docs/reviews/` and the point maps in
`docs/DOCS_ASSESSMENT.md` (Stage 3, R1.1-R1.9 / R2.1-R2.14).

# Layout

| Path | What it is |
| --- | --- |
| `main.tex` | Manuscript root: preamble, title/abstract, section inputs, appendix. |
| `refresponse.tex` | Response root (mirrors the manuscript setup); inputs the three referee files. |
| `sections/(N)name.tex` | One file per section, numbered input order; `(A)`-`(E)` are the appendices. |
| `references.bib` | Bibliography (currently mirrors `paper/BandF.bib`; trim to the revision's citations). |
| `tables/` | The four publication-ready flow tables (`matrix_5x3_v10_flows.md` + the three supply-arm designs), moved from `paper/tables/`; linked from Sections 5, 6 and App. C. |
| `figures/` | Figure files; wired via `\graphicspath`. |
| `response/{editor,reviewer1,reviewer2}.tex` | One file per referee, prepopulated point by point. |

# Carried-over material

Passages that exist only in the dissolved old adaptation draft (not in the
validated submission) are embedded as comments at their point of use in
Sections 2, 4, 5 and 6, each under a `CARRIED OVER FROM AN OLD DRAFT` headline.
Two title suggestions sit as a comment above `\papertitle` in `main.tex`.

# Section map

| Sect. | Title | Source |
| --- | --- | --- |
| 1 | Introduction | stub |
| 2 | From transformation research to IO models | stub |
| 3 | Two primitive questions | stub |
| 4 | One kernel, five closures, three financings | stub |
| 5 | Application: the programme through the matrix | stub |
| 6 | Beyond demand-only shocks | stub |
| 7 | Discussion | stub |
| 8 | Conclusion | stub |
| A | Full proofs | ported from `equivalence.tex` v4 (842 lines, 8 subsections) |
| B | Kernel, calibration, A-bill accounting | stub |
| C | Matrix variants and gamma/beta robustness | stub |
| D | Vector construction and data provenance | stub |
| E | Closures and the neo-Keynesian mapping | ported from `closures_and_demand_shocks.md` (NEW vs. the outline's A-D) |

# Compile (in-container, verified 2026-10-04)

```{.bash}
cd revision
nix shell nixpkgs#biber -c latexmk -pdf -interaction=nonstopmode main.tex
nix shell nixpkgs#biber -c latexmk -pdf -interaction=nonstopmode refresponse.tex
```

Produces `main.pdf` (26 pp, no errors, 0 undefined references) and
`refresponse.pdf` (5 pp).

# Toolchain (container)

One-time setup, already done --- the user tree is `TEXMFHOME=/root/texmf`
(Mac-backed, persists):

```{.bash}
tlmgr --usermode init-usertree
tlmgr --usermode option repository \
  https://ftp.math.utah.edu/pub/tex/historic/systems/texlive/2025/tlnet-final
tlmgr --usermode --verify-repo=none install \
  biblatex endfloat todonotes framed logreq appendix preprint lscape booktabs
```

Notes: the live CTAN repo is TeX Live 2026 and is refused cross-release, so the
frozen 2025 archive is required; `authblk` ships inside the `preprint` package;
`biber` cannot be installed in user mode --- it comes from `nixpkgs#biber`. The
store's `texlive-combined-medium-2025-final` is the base and lacks
authblk/biblatex/todonotes/endfloat/biber.

# Open decisions (from the outline's external review, v3)

- `AGENTS.md` refresh (protected file; needs interactive approval).
- Grouping rule ratification (`docs/grouping_rule_evidence.md`) before the
  rigid-group cells in Section 6 are citable.
- psi composition (raw 2024-price incidence vs.\ constant-2019-price).
- Appendix E carries a reviewer flag: its closing dominance wording must be
  checked against the Stage-3 must-not-appear list.