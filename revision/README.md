---
title: "Revision draft --- BeyondHulten / Metroeconomica"
project: "BFRep/(3)BeyondHulten"
created: 2026-10-04
status: "compiles via build.sh (main.pdf 66 pp, refresponse.pdf 12 pp, 0 undefined refs); build script added 2026-10-22"
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
| `sections/(N)name.tex` | One file per section, numbered input order; `(A)`-`(I)` are the appendices. |
| `references.bib` | Bibliography (currently mirrors `paper/BandF.bib`; trim to the revision's citations). |
| `tables/` | All tables: the four flow tables (`matrix_5x3_v10_flows.md` + the three supply-arm designs, from `paper/tables/`; linked from Sections 5, 6, App. E) and `table1_closure-overview.tex` (Table 1, `\input` by the introduction, with its crafting comments). |
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
| Sect. | Title | Notes |
| --- | --- | --- |
| 1 | Introduction |  |
| 2 | From transformation research to input--output models | programme and its import content |
| 3 | Two primitive questions, one closure matrix |  |
| 4 | A shared machine room: the kernel, closures and financing |  |
| 5 | Application: the programme through the matrix |  |
| 6 | Beyond demand-only shocks | Kaldor--Verdoorn, markup door, no-go |
| 7 | Conclusion | absorbs the dissolved Discussion; no separate Sec. 8 |
| A | Equivalences, overlaps and scope of the closure matrix |  |
| B | The model, its closures and the formal results | B.1 kernel, B.2 labor closures, B.3 financing, B.4 closure systems, B.5 demand-only reduction, B.6 elastic supply, B.7 scope |
| C | From the input--output table to Omega: the data pipeline | calibration, A-bill accounting, consistency checks; C.9 = the composition of the programme's spending |
| D | Sectoral responses in every closure | figure panels only, no numbered subsections |
| E | The pricing door: derivations, guard and accounting |  |

# Compile (in-container and on the Mac, verified 2026-10-22)

`build.sh` builds both roots, strips stale aux files first, and prints a
verification line for each PDF:

```{.bash}
cd revision
bash build.sh                 # both PDFs
bash build.sh main            # manuscript only
bash build.sh refresponse     # response letter only
```

Both roots need `biber` on PATH: the manuscript and the response letter each
carry a biblatex bibliography over `references.bib`. Building happens in place,
so `main.synctex.gz` keeps matching `main.pdf`. Equivalent by hand:

```{.bash}
cd revision
nix shell nixpkgs#biber -c latexmk -pdf -interaction=nonstopmode main.tex
nix shell nixpkgs#biber -c latexmk -pdf -interaction=nonstopmode refresponse.tex
```

Last verified 2026-10-22: `main.pdf` 66 pp (0 errors, 0 undefined references,
3 overfull boxes) and `refresponse.pdf` 12 pp (0 errors, 0 undefined
references, 0 overfull). `build.sh` exits non-zero on any unresolved
reference, and on the response letter's own canary: a broken `Location.` line,
which an unbraced `\loc` call produces.

# The response letter

`refresponse.tex` inputs one file per addressee (`response/editor.tex`,
`response/reviewer1.tex`, `response/reviewer2.tex`). Each numbered point is one
cell, in this order:

- `\subsection*{R<n>.<m> --- Title}` for the point;
- `\refcom{...}` quoting the referee verbatim and in full: one joint block per
  point, whole paragraphs, nothing elided;
- `\response ...` for our answer;
- `\loc{...}` naming the location of the change. The braces are required: an
  unbraced `\loc Sec. 1;` swallows the first character into the bold label.

Citations use the shared `references.bib` and the manuscript's biblatex
options. Location lines are provisional until the manuscript is frozen.

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
- Appendix G carries a reviewer flag: its closing dominance wording must be
  checked against the Stage-3 must-not-appear list.