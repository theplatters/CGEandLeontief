---
title: "Document Review Leftovers"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-09-14"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [archive, docs-audit, revision, beyondhulten, metroeconomica]
---

**Version 1** (September 2026)

\textcolor{revisionV2}{Archival note: the sections below were removed
from \texttt{docs/DOCS\_ASSESSMENT.md} at its Version 3 and retained
here as the chronological audit that produced the working plan. The
planning documents 1--8 discussed below live in \texttt{docs/archive/};
\texttt{definitive\_guide.md} and \texttt{labor\_closures.md} remain
active.}

# Chronological document map

| # | Document | Date | Status in one line |
|---|----------|------|--------------------|
| 1 | `REVISED_SALVAGE_PLAN.md` | 2026-07-29 | Foundational strategy; superseded in its results, still the best review-point map |
| 2 | `WORKPLAN.md` | 2026-07-29 | Superseded; historical record only |
| 3 | `WORKPLAN2.md` | 2026-07-29 (upd. later) | Superseded; its "GO / 88.4%" round was formally rejected |
| 4 | `OUTLOOK_SUMMARY.md` | 2026-08-02 | Superseded; abstract draft is built on invalid numbers |
| 5 | `PRELIMINARY-ASSESSMENT.md` | 2026-08-18 | The audit that killed the false claims; bug documentation still valuable |
| 6 | `milestone_D_plan.md` | 2026-08-19 | Executed and complete; historical |
| 7 | `WORKPLAN3.md` | 2026-08-18 (addendum 2026-09-03) | Repair record; its 2026-09-03 addendum holds canonical Part I numbers |
| 8 | `accounting_consistency_plan.md` | 2026-09-02 | Executed and integrated; now reference documentation |
| 9 | `definitive_guide.md` | 2026-09-03 | **Current authoritative narrative**; the manuscript should be written from it |
| 10 | `labor_closures.md` | 2026-09-05 | Current code semantics; contains an open design question that needs a decision |

External anchors:

- `ROADMAP.md` (2026-08-18/20): the governing plan, validation gates, and
  definition of done. Not superseded by anything in `docs/`.
- `roadmaps/vertdict.md` (2026-08-18): the independent audit verdict that
  rejected the first mobile-labor round. It is the arbiter WORKPLAN2 defers to.
- `docs/reviews/metro-rev1.docx` and `metro-rev2.docx`: the referee reports
  (R1 = constructive major revision; R2 = devastating reject). Extracted text
  is available at `docs/reviews/metro-rev1.txt` / `metro-rev2.txt` if needed.

\textcolor{revisionV1}{Note (v2): documents 1--8 of the map now live in \texttt{docs/archive/}; only \texttt{definitive\_guide.md}, \texttt{labor\_closures.md}, this assessment, and \texttt{reviews/} remain in \texttt{docs/} directly.}

# Evaluation, oldest first

## 1. REVISED_SALVAGE_PLAN.md (2026-07-29)

The foundational pivot document. It contains three things of unequal value.

**Still relevant:**

- **Section 2, the review-point map** (R1.1-R1.9, R2.1-R2.14) is the most
  complete mapping of referee demands to revision actions anywhere in the
  repository. Every response letter will be checked against it. Keep.
- The strategy insight of section 3.1 ("closure dominance claim") survives in
  modified form: acknowledge that the IO-CGE bridge is known (Robinson 2006,
  Rose 1995) and contribute measurement, not discovery. The corrected model
  qualifies the finding (see doc 9), but the *strategic posture* -- graceful
  concession to R2 plus a quantitative contribution -- is unchanged.

**Superseded / must not resurface:**

- The claim that "the core infrastructure exists and works" (sections 1.1
  and 7) was written before the mobile-labor audit and is wrong for the
  files it praises (`mobile_labor.jl`, `variance_decomposition.jl` at that
  date).
- The proposed labor-supply function `L = L_bar * (w/w_bar)^eta` (eta as
  labor-supply elasticity, eta in [0, inf)) is **not** what was finally
  implemented. Current code uses eta as a geometric intersectoral
  reallocation parameter with eta in [0, 1] (see doc 10). Any sentence in
  the manuscript drafted from this plan must be re-checked against
  `labor_closures.md`.
- The 6-10 week timeline and the phase plan are dead; the corrected critical
  path is in `ROADMAP.md` section 10.

**Dispositions for remaining work:** extract the review-point map into the
response-letter skeleton; otherwise treat as archive.

## 2. WORKPLAN.md (2026-07-29)

The three-stream coordination document (A: manuscript, B: model extension,
C: B&F replication).

**Still relevant:** only the *topology* -- three parallel streams with the
go/no-go gate between model work and writing -- and the Stream C pointer to
`bf_replication/REPLICATION_WORKPLAN.md`. The stream structure is inherited
by every later plan, so the file has organizational value as the origin of
the scheme.

**Superseded:** all statuses, all timelines, the priority matrix, and the
paths (it references `(1)Submission/revised/`, which now lives at
`revised_manuscript/` inside this repo). Its per-task efforts fed WORKPLAN2's
false "complete and verified" certification, which `vertdict.md` rejected.

**Disposition:** archive. Nothing unique except the stream topology.

## 3. WORKPLAN2.md (2026-07-29, updated)

The post-pilot workplan. Its headline claims -- "complete and verified", GO
decision, price invariance, eta = 88.4% -- were **formally superseded** by
the supersession note at the top of the document itself, which cites
`roadmaps/vertdict.md`. The document is internally annotated: it now reads as
a historical record of the rejected round.

**Still relevant:**

- The supersession note itself is a model of how the project handles
  retracted results; it should be cited in the response letter's
  methodological-correction section if the authors choose to disclose the
  internal audit trail.
- The Julia environment alignment record (1.12.6, container/host portability)
  remains operationally true.
- The destructuring-bug fix note (named-tuple destructuring of structs) is a
  real Julia pitfall worth remembering for future code review.

**Superseded:** every checkmark/GO/88.4% marker, the "Stream B results are
ready, Sections 5-6 can proceed immediately" instruction, and the
eta-as-labor-supply semantics throughout.

**Disposition:** archive, clearly marked as the rejected round. Do not use
for the manuscript.

## 4. OUTLOOK_SUMMARY.md (2026-08-02)

The response strategy essay ("From Searching for Bridges to Measuring the
Bridge") plus a draft abstract.

**Still relevant:**

- The rhetorical frame of the pivot is sound and still governs: retire the
  commensurability language, concede R2's point, contribute measurement.
- The priority matrix's pending Stream A tasks (introduction, Section 2
  rewrite, literature integration, figures, response letter) are an accurate
  list of what is *still* pending today -- Stream A has not advanced since.

**Superseded:**

- The draft abstract is quantitatively dead: it reports the eta = 88.4% and
  "order of magnitude" claims and the eta-as-labor-supply-elasticity reading
  of `mobile_labor.jl`. It must not be used even as a stylistic base without
  replacing every number and the eta interpretation.
- The claim "Done. The `mobile_labor.jl` module implements
  `L = L_bar * w^eta`" describes the since-rejected specification.

**Disposition:** keep the pivot rhetoric in mind; discard the abstract draft.
The abstract will be written anew from `definitive_guide.md` Part III.

## 5. PRELIMINARY-ASSESSMENT.md (2026-08-18)

The independent audit that caught the false certification. Verified by
execution that: the eta-sweep response was inverted (eta = 10 collapsed GDP
to 12%), equilibria were initialization-dependent, the variance decomposition
renormalized main effects, Cobb-Douglas grid points threw `DomainError`, and
the ROADMAP-vs-workplans contradiction existed.

**Still relevant -- genuinely, not just historically:**

- Its bug catalog is the reference description of *why* the old numbers were
  wrong. If the response letter includes a methodological-correction section
  (recommended by `definitive_guide.md` Result 1), this document supplies the
  technical content.
- The verification table mapping `vertdict.md` code-line claims to source
  lines is a completed audit trail.
- Its corrections section (overspend was 4.5% in the 71-sector model, not
  25.5%; the decisive fix is baseline-wage anchoring, not nominal-vs-real)
  prevents repeating two specific misreadings.

**Superseded:** nothing. It was the superseding document. The deleted
`tests/minimal_test/` diagnostics are noted; the organized test suite now
covers the same ground.

**Disposition:** keep as the audit record. Source for the response letter.

## 6. milestone_D_plan.md (2026-08-19)

A short task plan for the honest Sobol decomposition. Every step (first-order
and total-order indices, absolute shares, CSV output, robust missing-point
handling) was implemented and verified (WORKPLAN3 Milestone D plus the
2026-09-03 addendum).

**Still relevant:** only as the specification-of-record for what
`summary_table` and `SobolResult` are supposed to output. If the sectoral
Sobol extension (still pending, see doc 9 Part V) is built, this plan is the
template.

**Disposition:** archive; consult before extending the sensitivity code.

## 7. WORKPLAN3.md (2026-08-18, addendum 2026-09-03)

The repair workplan (Milestones A-F) with the critical 2026-09-03 re-run
addendum. This is where the *current canonical Part I numbers* live.

**Still relevant -- load-bearing:**

- **The 2026-09-03 wage-regime table** (flexible eta=0 turning negative with
  shock size; sticky `:fixed` positive and large, +19.3 pp at mult 10;
  eta=1 stalls from mult 0.5): this is Result 3 of the corrected story and
  the single most important empirical table for the revised paper.
- **The 2026-09-03 Sobol result** (theta = 0.395 dominant, sigma = 0.273,
  epsilon = 0.165, eta = 0.157 material-not-dominant): this is Result 2 and
  the numbers the manuscript will report -- *not* the 88.4%, *not* the "~2%
  negligible" reading, both of which are documented artifacts.
- The accounting canon: GDP P = I = 3,027,818; E = 2,864,724; residual
  5.387% documented valuation gap; `sum(lambda) = 2.1099`.
- Milestone A's list of equilibrium repairs (sector-1 zero-profit restored,
  normalized household demand, retcode checking, Tornqvist index) is the
  methodological-correction content again, from the fix side.

**Superseded:** the historical Milestone D (eta around 0.01%) and Milestone F
(74-112x ratios) numbers are explicitly superseded by the addendum; the
addendum's own supersession note is unambiguous.

**Disposition:** keep as the numerical record. The manuscript's empirical
section is built from this file plus `definitive_guide.md`.

## 8. accounting_consistency_plan.md (2026-09-02)

The section-4.1 accounting transformation (ROADMAP Phase 1): separate
imports, decompose value added, reconcile the three GDP sides, define shock
incidence, emit calibration artifacts. Executed, integrated into
`src/interface.jl`, and validated by the test suite.

**Still relevant:**

- It is the reference for the data pipeline: the verified source schema
  (rows 73-83 of the Destatis table), the proportional import-allocation
  assumption (stated explicitly), the 5.387% expenditure residual (a genuine
  raw-table valuation gap, not a bug -- after the off-by-one indexing fix),
  and the field mapping in the integration note.
- The manuscript's "Data and accounting" section (structure item 3 in both
  ROADMAP section 7 and definitive guide Part III) will be written almost
  directly from this document. It answers R2's "closed-economy miracle"
  demand with a documented, open-economy-consistent pipeline.
- The calibration-table demand of R1 is answered here (Step 6 artifacts,
  `output/AC_*.csv`, regenerable).

**Superseded:** nothing material. The ~0.8% residual figure was corrected in
place to 5.387%.

**Disposition:** keep as the data-pipeline reference.

## 9. definitive_guide.md (2026-09-03)

The current authoritative synthesis. It supersedes WORKPLAN/WORKPLAN2,
operationalizes the ROADMAP's Phase 6 go/no-go with evidence, and defines the
corrected three-part story:

1. **Methodological:** the original model was not an equilibrium; a
   corrected, verified baseline exists (machine-precision residuals).
2. **Negative but precise:** eta (intersectoral mobility) is a material but
   secondary channel (first-order Sobol share around 0.157, well below
   theta = 0.395); the mobility bridge is second-order (at most 0.07 pp);
   eta=1 solves only at multipliers 0.1/0.2.
3. **Positive:** the wage regime (sticky vs flexible) is the first-order
   labour-market margin (pp gaps +0.32 to +24.75); the question is "is the
   economy at full employment or in an unemployment regime?", not "how
   mobile is labour?".

**Still relevant -- this is the master document:**

- Part III's manuscript structure (9 sections), headline abstract results,
  the "remove/add" lists, and the ROADMAP compliance table.
- Part V's remaining-work list, which is the best available backlog:
  - Sobol on **sectoral** quantities (not yet run);
  - ROADMAP Phase 2 policy experiment (financing closure) -- **not started**;
  - Milestone C solver stability at the Cobb-Douglas limit -- open;
  - manuscript rewrite -- not started;
  - response letter -- not started;
  - `:fixed` with an explicit demand anchor (to show the sticky-wage result
    is not eta-driven) -- open.
- Part IV records that Closure D (IO endpoint) and Closure C
  (unemployment complementarity) are unimplemented.

**Superseded:** nothing; it internally supersedes the 2026-09-02 re-run
numbers with the 2026-09-03 ones.

**Disposition:** active. All manuscript writing starts here. Its remaining
work items should become the working task list.

## 10. labor_closures.md (2026-09-05)

The latest document, and the one that defines what the code *actually does*
now.

**Still relevant -- critical for correctness:**

- The three-dimension taxonomy: legacy exogenous `labor_slack` callback;
  geometric eta reallocation (`L_fixed^(1-eta) * L_costmin^eta`, eta in
  [0,1], beyond which is extrapolation); wage regime `:mobile`/`:fixed`.
  These must never be conflated in the manuscript -- and note that **the
  current eta is NOT the labor-supply elasticity the referees and the
  salvage plan talked about.** The paper's terminology must be chosen
  deliberately: either rename the model's parameter (e.g. "reallocation
  parameter") or implement the standard partial-mobility / elastic-supply
  formulation and map the referee's eta onto it.
- The proposed **standard partial-mobility formulation** (sectoral wages,
  wage-responsive sectoral supply, 3N system) is an *open design proposal*,
  not implemented. It is the main unresolved modeling decision: the current
  geometric-eta plus efficiency-penalty design is a project-specific reduced
  form whose allocative-wedge curvature is not established as the exact CES
  allocative-loss coefficient. This directly touches R2's demand for
  labor-supply elasticities ("the one obvious set of elasticities that
  really matters").
- The warning that `:fixed` is not a capped unemployment model (employment
  can exceed `labor_bar`; the gap is computed post-solve) -- essential for
  honest interpretation of the +19.3 pp sticky-wage result. The sticky-wage
  closure is currently a *fixed real wage with unconstrained employment*,
  i.e. closer to Closure D's spirit than to Closure C. The manuscript must
  not sell it as a calibrated unemployment closure without saying so.

**Disposition:** active. The eta-semantics decision and the
sticky-wage-interpretation caveat feed directly into Sections 4 and 6 of the
new manuscript.

# Synthesis: where the project stands

## Verified and usable today

- A corrected, equilibrium-consistent 71-sector model (Milestones A, B and
  the 2026-09-03 fixes), with organized tests in `tests/`.
- Section-4.1 accounting-consistent data pipeline integrated in
  `src/interface.jl`.
- Canonical Part I results: wage-regime table, Sobol decomposition
  (theta-led, eta material at about 15.7%), accounting reconciliation -- all
  from the 2026-09-03 re-run, regenerable via
  `julia --project=. rerun_results.jl`.

## Open items (candidates for the working plan)

1. **eta-semantics decision** (doc 10): keep geometric eta with honest
   renaming, or implement the standard partial-mobility / elastic-supply
   closure that R2 actually asked for. This blocks the framing of Sections
   4-6.
2. **Sticky-wage closure interpretation** (doc 10): `:fixed` = fixed real
   wage with endogenous employment, not a calibrated unemployment regime;
   decide whether Closure C (complementarity) is implemented or the
   limitation is stated.
3. **Financing closure** (ROADMAP Phase 2; definitive guide Part V.5): the
   demand shock is still unfinanced -- an explicit closure is required
   before CGE results are admissible by the ROADMAP's own rule. This is the
   largest *unstarted* modeling item.
4. **Sectoral Sobol** (definitive guide Part V.3): does eta matter for
   sectoral allocation even if not for aggregate GDP? R1 expects the
   aggregate/sectoral contrast to be quantified.
5. **Cobb-Douglas limit** (Milestone C): epsilon = 0.99 DomainError guard.
6. **eta=1 stall under autonomous demand**: needs the explicit
   demand/investment anchor or an honest limitation note.
7. **Manuscript rewrite** per definitive guide Part III / ROADMAP section 7
   -- `revised_manuscript/` chapters still contain only the old text
   ("% OLD TEXT FROM HERE ON").
8. **Response letter** (R1 + R2), to be checked against the R1.1-R1.9 /
   R2.1-R2.14 map in `REVISED_SALVAGE_PLAN.md` section 2.
9. Minor copy-edits from R1 (p.5 sentence, p.16 paragraph, Figure 3 axis).

## Suggested next decision

Item 1 (eta semantics) is the fork: it determines whether the paper says
"the reallocation parameter eta is material but secondary" (what the code
now supports) or "the labor-supply elasticity continuum behaves as follows"
(what R2's framework expects and what would require the standard
partial-mobility implementation from doc 10). Everything downstream --
figures, Sobol labeling, abstract -- inherits from this choice.

\textcolor{revisionV1}{(v2) These decisions are now laid out with full option tables in Foundation I and Foundation II at the top of this document; the workplan stages supersede the open-items list here.}
