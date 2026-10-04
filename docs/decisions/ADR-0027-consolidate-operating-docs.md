# ADR-0027: Consolidate operating docs (one question, one document)

**Status:** Accepted (2026-10-04)
**Supersedes:** —
**Superseded by:** —

## Context

The repository accumulated ~14 operating documents (~5,100 lines) answering the
same questions in several places: plans lived in both `ROADMAP.md` and
`docs/DOCS_ASSESSMENT.md` §5 (and the now-archived `WORKPLAN_SENSITIVE_PRICES.md`);
results were split across `CONCISE_SUMMARY.md`, the archived
`CONCISE_SUMMARY_MERITOFSUPPLYSHOCKs.md`, `ETAs.md`, and the assessment; and the
closure prose in `labor_closures.md` duplicated `registry/closures.toml` (the
single source of truth, ADR-0005). Drift was already visible: `DOCS_ASSESSMENT`
§5 still lists Stage-1 items as "pending" / "v3 re-run" although they are
executed.

## Decision

Reduce the operating set to **one document per question**:

- **Contract:** `AGENTS.md`.
- **Plan and validation gates:** `ROADMAP.md` — the only plan (the house rule;
  `DOCS_ASSESSMENT` §5 is retired to a pointer).
- **Evidence, verdicts, reviewer maps:** `docs/DOCS_ASSESSMENT.md` (drops the
  duplicate workplan).
- **Results summary:** `docs/CONCISE_SUMMARY.md` — the supply-shock note is merged
  into it as a new section.
- **Supply-arm method:** `docs/ETAs.md` (annex).
- **Narrative plan:** `paper/narrative_outline.md` (paper layer).

Actions taken under this ADR:

1. `CONCISE_SUMMARY_MERITOFSUPPLYSHOCKs.md` merged into `CONCISE_SUMMARY.md` and
   moved to `docs/archive/`.
2. `WORKPLAN_SENSITIVE_PRICES.md` moved to `docs/archive/` (its subject is
   executed: `matrix_5x3_v9/v10` + probes 15–19; the one forward-looking piece,
   the wage-structure door, now lives in `docs/VariationinGamma.md`, ADR-0021).
   Its eight inbound citations are repointed to the executed evidence
   (`docs/VariationinGamma.md`, `docs/ETAs.md`, `docs/grouping_rule_evidence.md`,
   `paper/tables/matrix_5x3_v10_flows.md`).
3. `labor_closures.md` moved to `docs/archive/` (redundant with
   `registry/closures.toml`, ADR-0005); its one live status reference in
   `docs/DOCS_ASSESSMENT.md` is updated to "implemented".

**Pending (environmental conflict).** A concurrent draft session overwrote the
working-tree edits to `DOCS_ASSESSMENT.md`, `ETAs.md`, `ROADMAP.md` and
`paper/narrative_outline.md` while this consolidation ran — those four files are
currently reverted to the draft session's versions. The following (3) actions
are therefore still outstanding and must be redone once the draft session is
paused:

- `docs/DOCS_ASSESSMENT.md` §5 retirement (and its `labor_closures.md` active
  reference / "not implemented" status line).
- `docs/ETAs.md` dead-pointer fixes (`supply_etas.toml`, `probe5_supply_arm.jl`).
- `ROADMAP.md` Phase-4 ticks and thesis/title refresh.
- `paper/narrative_outline.md` cell-count fix (126 → 189) and the Stage-3
  must-not-appear import.

## Consequences

- The live operating set drops from 14 documents to ~10; one place to look per
  question, and no second plan.
- New documents must declare their host role and an expiry trigger (the event
  after which they are archived), so the set cannot regrow by drift.
