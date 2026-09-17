# Archive

**Closed zone.** Everything in this directory is a historical planning or
assessment document, superseded as a statement of current status by
`docs/status.md` (generated from `registry/`, see
`docs/decisions/ADR-0003-registry-single-source-of-truth.md`). The governing
plan remains `ROADMAP.md` and `docs/DOCS_ASSESSMENT.md`.

Do not update files here and do not add new ones except to record what
superseded them. If a document here is still cited as current, fix the
citation instead of editing the document. Documents archived under ADR-0008
are also registered read-only in `registry/freeze.toml`, so an edit shows up
as a freeze-board warning and a `scripts/check_repo.jl` violation.

Contents include: `WORKPLAN.md`, `WORKPLAN2.md`, `WORKPLAN3.md`,
`REVISED_SALVAGE_PLAN.md`, `PRELIMINARY-ASSESSMENT.md`,
`OUTLOOK_SUMMARY.md`, `DOCS_BACKUPv2.md`, `document_review_leftovers.md`,
`accounting_consistency_plan.md`, `milestone_D_plan.md` (and PDF mirrors),
plus the archived root documents `selective_status_overview.md` (a superseded
status overview) and `varianten.xlsx` (the pre-Phase-0 closure-design
workbook, renamed from `varianten`; both moved here by ADR-0008), and the
root legacy notebooks `DemandShocks.ipynb`, `CompareModels.ipynb` and
`CobbDouglas.ipynb` in `notebooks/` (moved here by ADR-0009, still frozen).
