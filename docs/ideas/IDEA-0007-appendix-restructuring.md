---
title: "IDEA-0007: restructuring the remaining appendices --- four by function instead of six by topic"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-10-18"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [ideas, appendices, restructuring, revision, writing-plan]
last-updated: October 2026
---

**Version 1** (October 2026)

Status: `scoped`.

# The question

How to structure the six appendices that follow the two now settled ones. Appendix A
(statements: the equivalences, their mechanisms, their scope) and Appendix B (the
formal apparatus: notation, kernel, closures, proofs, numerical evidence) were merged
and rewritten on 2026-10-18 and are not part of this idea. What remains is letters C
to H, and three of those six are stubs, one is written but nearly unreachable, and the
most-cited appendix in the paper is 367 words long.

# What exists (measured 2026-10-18)

Counted from `revision/sections/`; citations counted from the main text (Sections 1-8)
and, where noted, from within the appendix block.

| Letter | Title | Words | Structure | Cited from | State |
| --- | --- | --- | --- | --- | --- |
| C | Kernel, calibration, A-bill accounting | 172 | no subsections, no tables | Section 4 (2), B (2) | stub, 2 TODO notes |
| D | Matrix variants and gamma/beta robustness | 367 | 1 subsection, 1 table input | Section 5 (1), Section 6 (5), B (2) | thin skeleton, 1 TODO note |
| E | Vector construction and data provenance | 95 | none | Section 6 (1) | placeholder, 1 TODO note |
| F | Closures and the neo-Keynesian mapping | 1,712 | 5 subsections, 1 table | Section 3 (1) | written; 2 ADR and 2 `\texttt` mentions left |
| G | From the input--output table to Omega: the data pipeline | 861 | 7 subsections | Section 4 (1) | scaffolded, 19 TODO notes |
| H | The pricing door: derivations, guard and accounting | 883 | 5 subsections | Section 6 (6, incl. sub-labels) | written; 3 ADR and 3 `\texttt` mentions left |

The block arithmetic is the diagnostic. The appendix block is 12,851 words and
Appendix B alone is 7,563 of them, 59 percent. The six remaining appendices total
4,090 words, and the three smallest of them --- C, D and E, 634 words together ---
absorb 13 references (9 from the main text). Writing effort and citation density run
in opposite directions across the block, and the letters are named for topics while
the readers' questions are functional.

# Three mismatches

- **Promise debt.** Section 4 twice sends the reader to the intermediate-bill and
  A-bill accounting (the leaks sentence in 4.1 and the seven-component deflator
  footnote) and once to the data pipeline. Appendix C, which carries the accounting,
  is 172 words with two TODO notes; the pipeline sits in a different letter.
- **Reach.** D is cited eight times and contains 367 words; F is the only fully
  written appendix of the six (1,712 words) and is reached by a single pointer. F is
  the R2-facing material --- the demand-shock cleavage between CGE and input--output
  practice, the FIDELIO comparison, the mapping of neo-Keynesian assumptions onto the
  closures --- and a referee currently has one route into it.
- **Function.** C, E and G all answer "where do the numbers come from"; D and G's
  last subsection both answer "what happens under a variant"; H answers "how does the
  pricing door work"; F answers "how does this relate to the traditions". The letters
  do not mark those four questions.

# The proposal: four appendices by function

| New letter | Appendix | Merged from | The question it answers | Estimated length |
| --- | --- | --- | --- | --- |
| A | Equivalences, overlaps and scope | A (settled) | what holds, without proofs | 1,200 |
| B | The model, its closures and the formal results | B (settled) | the formal apparatus | 7,560 |
| C | Data and calibration: from the published table to the model | C + E + G | where the numbers come from | 2,500 - 3,000 |
| D | Robustness: variants, profiles and the matrix | D, expanded | what survives a variant | 1,500 - 2,000 |
| E | The pricing door | H | how demand reaches prices | 880 |
| F | Closures and the neo-Keynesian mapping | F, plus a pointer | how this relates to the traditions | 1,700 |

Eight appendices become six, and the six remaining become four. Each letter then
carries one question, and the sequence follows the paper rather than the order in
which the material was written: statements, proofs, data, robustness, the door,
the traditions. B's internal composition is not part of this idea.

# Move 1 --- one data and calibration appendix (C + E + G)

Order the content as the reader would audit it: the source table and what is read from
it, the intermediate-use matrix, the normalisation to cost shares, value added, factor
shares and Domar weights, the final-demand and household shares, the consistency
checks, the construction of the impulse vector, and only then the intermediate-bill and
A-bill accounting with the calibration.

- Why: a single audit trail for Sections 4 and 5, in one letter, in the order the
  numbers enter the model.
- What it repairs: Section 4's two accounting pointers land in an appendix that
  contains the accounting; the 172-word stub and the 95-word vector note stop being
  appendices in their own right.
- What it costs: filling G's nineteen TODO notes (the largest single writing block in
  the restructure) and writing the A-bill core that C currently promises in two notes.
- Open decision: whether the sector-exclusion variants stay here or move to D. The
  test is which question the variant answers --- if it changes *which sectors enter the
  model* it is a data decision and stays in C; if it changes *what the matrix reports*
  it is a robustness check and belongs in D.

# Move 2 --- the robustness appendix (D, expanded)

D should hold the fifteen-cell matrix, the elasticity ladder, the pricing-door profiles
it already inputs, and the gamma/beta representation checks --- the places where the
matrix's nine economies and the window over the sectoral elasticities are shown to be
insensitive to representation choices.

- Why: Section 6 forwards five separate things here (the variant ladder, the mobile
  rows, the profiles, the sectoral cells, the elastic-supply benchmark). A referee who
  checks Section 6's claims lands in the appendix that carries the paper's central
  robustness argument, and it currently holds 367 words.
- What it buys: the robustness claim becomes auditable rather than asserted, and D
  stops being the thinnest appendix with the most inbound references.

# Move 3 --- the pricing door keeps its own letter (H becomes E)

No structural change, only the register purge (3 ADR and 3 `\texttt` mentions).

- Why: six references to 883 words is the best citation-per-word ratio in the block,
  and its derivations sit behind two subsections of Section 6. It is the one appendix
  of the six that already earns its place at its current size.
- Alternative worth weighing: fold it into the robustness appendix as one subsection.
  I would not: it answers a different question from D (a derivation, not a variant),
  and it is the one place where the paper's demand-sensitive-pricing claim is proved
  rather than measured.

# Move 4 --- give the traditions mapping a reach (F)

Add one pointer from Section 7 (and, if the discussion grows a literature passage,
from there rather than from 3). No writing needed.

- Why: 1,712 words answered by one citation is the largest substance-to-reach gap in
  the block, and this is the appendix that speaks to the reviewer's question about how
  the closure choice relates to existing modelling traditions.

# Optional: where F sits

The one genuine placement choice in this idea. Both options keep the same content and
the same letter-independent labels; only the reading order changes.

| F stays last (slot F) | F moves to slot C, directly after the formal apparatus |
| --- | --- |
| The block ends on the traditions comparison, after data and robustness, which is where a reader who wants a relation to the literature will look once the results are known. | The comparison reads as the conceptual companion to Section 3, immediately after the formal work it maps onto. |
| No section reference changes; the one pointer added in Move 4 works from either place. | The data and robustness appendices (the C and D of the new structure) move one letter later. |
| Risk: the "appendix tail" reading --- the traditions material looks like an afterthought. | Risk: the reader meets a literature comparison before the calibration and the robustness checks, so the block's middle is conceptual rather than evidential. |
| Works as it stands. | Recommended only if Section 7 gains a passage on the modelling traditions, so that the appendix is the reader's *destination* for a question the main text raises. |

# Writing order and cost

- **D first**: highest return (eight references), smallest volume (about 1,200 words to add).
- **C second**: largest promise debt; about 1,500 words, most of it filling G's notes.
- **Cheap fixes alongside**: F's pointer (one sentence), H's register purge, and the
  sector-exclusion decision.
- **What cannot be avoided in any structure**: the calibration Section 5 reports has to
  be auditable, so G's nineteen notes become prose whichever letter they live in.

# What this buys, and what it does not

- It buys: one question per letter; no promise the main text cannot keep; the audit
  trail contiguous; and an appendix block that reads in the paper's own order.
- It does not buy: any content. Nothing here writes the missing paragraphs, and the
  restructure changes no claim, no number and no measurement --- it is a filing
  decision that makes the existing debt visible and payable.
- It also does not reopen what was settled on 2026-10-18 (A's rewrite and the A/B/C/H
  merge), and it assumes the labels keep resolving as they do now, which the merge
  demonstrated is safe because every reference is label-based.

# Adjacent items found while writing

- **The appendix counters carry the wrong prefix.** Tables and equations in the
  appendix block are numbered with a fixed "A" prefix, set once in `main.tex`
  (`\theequation = A\arabic{equation}` and the matching table and figure counters).
  So `eq:zp` prints as (A5) inside Appendix B, and Appendix B's own equations run A3,
  A4, A5, ... A23 while Appendix B's sections are lettered B. Table A1 in Appendix A is
  correct; an equation labelled (A5) inside Appendix B is not. Per-appendix counters
  would fix it in one line each.
- **The ideas register predates the re-lettering.** IDEA-0006's anchors cite "Appendix C
  (prop. `gd`, `ab`, `nest`)" (now B), "Appendix E" (now D) and "Appendix J" (now H);
  IDEA-0005's body names the pre-merge file set. The same class of stale-comment defect
  was corrected in six appendix file headers and three internal pointers on 2026-10-18.
- **Two labels are carried without a reader.** `app:model` is referenced by nobody (kept
  deliberately, so any external reference to "the model appendix" still resolves), and
  the sub-labels of F, G and H (`app:cleavage`, `app:fidelio`, `app:matrix-map`,
  `app:neokeynesian`, `app:summary`; the seven `app:pipeline-*`; `app:pricing-foc`,
  `app:pricing-cd`, `app:pricing-accounting`) are referenced nowhere, so no reader is
  being routed into those subsections.
- **This register's own index is inconsistent.** The version block at the top of
  `README.md` stops at Version 3 while the Revision Log runs to Version 5, and the log's
  entries are out of order (5, 1, 2, 4, 3, 4).

# Open decisions for the operator

- Adopt the four-by-function structure, or keep six letters and only fill the stubs in place.
- Sector-exclusion variants: data appendix (C) or robustness appendix (D).
- Does the pricing door keep its own letter, or become a subsection of D.
- F's slot: last as it stands, or after the formal apparatus (conditional on Section 7
  gaining a traditions passage).
- Whether the ideas register's stale appendix letters are corrected now or left as a
  historical record of the pre-merge numbering.
