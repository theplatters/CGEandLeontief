---
title: "Motivation mindmap - machine-readable reprint"
project: "BFRep/(3)BeyondHulten"
date: 2026-09-20
version: 3
status: "working note - not a decision, no registry entry"
tags: [motivation, mindmap, paradigms, revision, beyondhulten]
---

# Source and purpose

- Board photo `upload_20260920_144445_1.jpg`, hand-drawn blue marker on a whiteboard.
- Purpose: machine-readable reprint of the mindmap motivating the new paper (the revision) to be drafted.
- Center node: `MOTIVATION`.
- Caveats: a white cutout at the right board edge may hide content; the `WHICH GAMMA?` / `WHICH BETA?` cluster is deliberately drawn as an invitation to reflect, not as a pinned wiring.
- Conventions used below: `origin: board` = drawn on the board; `origin: user-expansion` = wording supplied by the author when correcting this reprint (not drawn); `status: confirmed` = agreed across reads/corrections; `status: open` = intentionally unresolved; `status: unconfirmed` = single-pass sighting, not asserted; `preliminary_idea:` = agent-suggested tie-in to the project vocabulary, not author-confirmed - treat as preliminary idea.

# Canonical graph

```yaml
version: 3
date: 2026-09-20
center: "MOTIVATION"
center_meaning: "the motivation of the new paper / revision to be drafted"

hunch:
  statement: >
    Preliminary hunch for the PARADIGMS framing (A-E):
    speak about B (vision of the short-run) and C (theory of money),
    frame A (adjustment closure) as a summary of B & C,
    mention D (price-effects of demand) and E (theory of production) in passing at best.
  priority:
    foreground: [B, C]
    summary: [A]        # A presented as the summary OF B & C
    passing: [D, E]     # at best

nodes:
  - {id: N01,  origin: board,            label: "MOTIVATION",               role: center,          meaning: "hub of the motivation story"}
  - {id: N02,  origin: board,            label: "PARADIGMS",                position: top,         meaning: "the paradigmatic perspectives framing the contribution (A-E below)"}
  - {id: N02A, origin: user-expansion,   label: "A - which adjustment closure?",                   meaning: "closure dimension; hunch = summary of B & C",
     preliminary_idea: "ties to the labour-closure family (BF/ALPHA/BETA/GAMMA/DELTA) and the 5x3 closure matrix"}
  - {id: N02B, origin: user-expansion,   label: "B - which vision of the short-run?",              meaning: "hunch = foreground"}
  - {id: N02C, origin: user-expansion,   label: "C - which theory of money?",                      meaning: "hunch = foreground",
     preliminary_idea: "ties to the financing closures (F1/F2/F3) and the money discussion"}
  - {id: N02D, origin: user-expansion,   label: "D - which price-effects of demand?",              meaning: "hunch = passing mention at best",
     preliminary_idea: "ties to the demand-sensitivity / price-block discussion (p = 1 under demand-only shocks; ETAs arm)"}
  - {id: N02E, origin: user-expansion,   label: "E - which theory of production?",                 meaning: "hunch = passing mention at best",
     preliminary_idea: "ties to technology / accounting (A-bill, Leontief), the IO core"}
  - {id: N03,  origin: board,            label: "BF?",                      position: upper-left,  meaning: "question: is the Baqaee-Farhi contribution & prominence the vantage point, as in the original framing?"}
  - {id: N04,  origin: board,            label: "PAPER DRAFT",              position: upper-right, meaning: "the old version of the paper as originally submitted"}
  - {id: N05,  origin: board,            label: "OLD DOCS",                 position: right,       meaning: "notes on what kind of results we have"}
  - {id: N06,  origin: board,            label: "(CLOSURE MATRIX 5x3",      position: below N05,   meaning: "the 5x3 closure-matrix results; leading '(' as drawn"}
  - {id: N07,  origin: board,            label: "NOTES IN /PAPER",          position: below N06,   meaning: "supporting notes live in the paper folder"}
  - {id: N08,  origin: board,            label: "TRANSFORM <-> VECTOR",     position: lower-center, meaning: "socio-ecological transformation research <-> our empirically informed input vector"}
  - {id: N09,  origin: board,            label: "IO-MODELS?",               position: below N08,   plural: true, meaning: "IO models (plural) as the bridge"}
  - {id: N10,  origin: board,            label: "IMPORTANT!",               position: below N09,   meaning: "emphasis: the realistic-use-case argument for IO models"}
  - {id: N11,  origin: board,            label: "WHICH GAMMA?",             position: left,        meaning: "5x3 presentation decision: whether to show these in the 5x3 table, which version to show"}
  - {id: N12,  origin: board,            label: "WHICH BETA?",              position: lower-left,  is_child_of: N11, meaning: "same decision, subordinate: which version; which alternatives to present in what follows"}
  - {id: N13,  origin: board,            label: "[EQUIVALENCES]",           position: bottom-left, is_child_of: N12, meaning: "the ex-ante-equivalence notes; narrative/framing of the contribution"}

edges:
  # ---- confirmed (board drawings, all reads agree) ----
  - {from: N01, to: N02, status: confirmed}    # MOTIVATION        -> PARADIGMS
  - {from: N02, to: N03, status: confirmed}    # PARADIGMS         -> BF?
  - {from: N02, to: N04, status: confirmed}    # PARADIGMS         -> PAPER DRAFT
  - {from: N01, to: N05, status: confirmed}    # MOTIVATION        -> OLD DOCS
  - {from: N05, to: N06, status: confirmed}    # OLD DOCS          -> (CLOSURE MATRIX 5x3
  - {from: N06, to: N07, status: confirmed}    # (CLOSURE MATRIX.. -> NOTES IN /PAPER
  - {from: N01, to: N08, status: confirmed}    # MOTIVATION        -> TRANSFORM <-> VECTOR
  - {from: N08, to: N09, status: confirmed}    # TRANSFORM <-> VECTOR -> IO-MODELS?
  - {from: N09, to: N10, status: confirmed}    # IO-MODELS?        -> IMPORTANT!
  - {from: N11, to: N12, status: confirmed}    # WHICH GAMMA?      -> WHICH BETA?
  - {from: N12, to: N13, status: confirmed}    # WHICH BETA?       -> [EQUIVALENCES]
  # ---- user-expansion (A-E under PARADIGMS; not drawn) ----
  - {from: N02, to: N02A, status: user-expansion}
  - {from: N02, to: N02B, status: user-expansion}
  - {from: N02, to: N02C, status: user-expansion}
  - {from: N02, to: N02D, status: user-expansion}
  - {from: N02, to: N02E, status: user-expansion}
  # ---- open: intentionally unresolved ----
  - {from: N11, to: "<ROOT>", status: open,
     note: "attachment of the WHICH GAMMA?/WHICH BETA? cluster to the rest of the graph is intentionally left open.
            Candidate anchors seen on the board: direct off MOTIVATION, descending off BF?, or off the closure-matrix branch.
            Semantics (v3): the cluster is the 5x3 presentation decision - whether to show these in the 5x3 table,
            which version of them to show, which alternatives to present in what follows."}
  # ---- unconfirmed: single-pass sighting, not asserted ----
  - {from: N12, to: N08, status: unconfirmed,
     note: "long line from the WHICH BETA? area across the board toward TRANSFORM <-> VECTOR seen on one read only."}
```

# Node table

| ID | Label (as drawn) | Context | Fn |
|---|---|---|---|
| N01 | `MOTIVATION` | Hub - the motivation of the new paper / revision | - |
| N02 | `PARADIGMS` | The paradigmatic perspectives framing the contribution: question-dimensions A-E below | [1] |
| N02A | `A - which adjustment closure?` | Closure dimension; hunch = present as summary of B & C | [1] |
| N02B | `B - which vision of the short-run?` | Hunch = foreground | [1] |
| N02C | `C - which theory of money?` | Hunch = foreground | [1] |
| N02D | `D - which price-effects of demand?` | Hunch = mention in passing, at best | [1] |
| N02E | `E - which theory of production?` | Hunch = mention in passing, at best | [1] |
| N03 | `BF?` | Open: is the Baqaee-Farhi contribution & prominence the vantage point, as in the original framing? | [2] |
| N04 | `PAPER DRAFT` | The old version of the paper as originally submitted | - |
| N05 | `OLD DOCS` | Notes on what kind of results we have | [3] |
| N06 | `(CLOSURE MATRIX 5x3` | The 5x3 closure-matrix results (leading `(` as drawn) | [3] |
| N07 | `NOTES IN /PAPER` | Supporting notes live in the paper folder | [3] |
| N08 | `TRANSFORM <-> VECTOR` | Socio-ecological transformation research <-> our empirically informed input vector | [4] |
| N09 | `IO-MODELS?` | IO models (plural) as the bridge | [4] |
| N10 | `IMPORTANT!` | Emphasis: the realistic-use-case argument for IO models | [4] |
| N11 | `WHICH GAMMA?` | 5x3 presentation decision: whether to show these in the 5x3 table, which version to show | [5] |
| N12 | `WHICH BETA?` | Same decision, subordinate - which version, which alternatives to present in what follows | [5] |
| N13 | `[EQUIVALENCES]` | The ex-ante-equivalence notes; narrative/framing of the contribution | [6] |

Notes

Tie-ins inside the canonical YAML (A-E) are marked `preliminary_idea`: agent-suggested links to the project vocabulary, not confirmed by the author.

[1] Author hunch (2026-09-20): speak about B and C; frame A as a summary of B & C; mention D and E only in passing, at best. A-E are not drawn on the board - they are the author-supplied expansion of the `PARADIGMS` node (`origin: user-expansion`).

[2] From the author's correction: in the original framing the Baqaee-Farhi contribution and its prominence was presented as a vantage point; the board node questions whether that vantage stays.

[3] `OLD DOCS` = notes on what kind of results we have; the results are the 5x3 closure matrix; the notes live in the paper folder.

[4] Confirmed 2026-09-20: TRANSFORM = socio-ecological transformation research; VECTOR = our empirically informed input vector; the argument gives a realistic use-case for applying IO models to transformation research (hence `IMPORTANT!`).

[5] Confirmed 2026-09-20: WHICH GAMMA? / WHICH BETA? = whether to show these in the 5x3 table, which version of them to show, and which alternatives to present in what follows. The board wiring of the cluster is intentionally left open.

[6] `paper/equivalence.tex` (v3): ex ante equivalence of the labour-closure endpoints under demand-only shocks (ALPHA = BETA, GAMMA = DELTA, mobile F2 = F3); the notes that matter for the narrative and framing of the contribution.

# Revisions

- 2026-09-20 (v1): first machine-readable reprint from the board photo; two vision-pass reconciliation; open questions asked.
- 2026-09-20 (v2): author corrections: IO-MODELS plural; TRANSFORM = socio-ecological transformation research; VECTOR = empirically informed input vector; BF? = Baqaee-Farhi lineage; [EQUIVALENCES] = `paper/equivalence.tex`; gamma/beta = representation reflection.
- 2026-09-20 (v3, current): PARADIGMS expanded into A-E with hunch priorities; BF? re-read as the vantage-point question; PAPER DRAFT = the old submitted version; OLD DOCS = notes on the results we have; WHICH GAMMA? / WHICH BETA? = the 5x3 presentation decision.