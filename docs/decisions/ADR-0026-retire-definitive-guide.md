# ADR-0026 — Retire `docs/definitive_guide.md`

- **Status:** accepted (operator decision, 2026-10-04)
- **Date:** 2026-10-04
- **Supersedes:** — (extends ADR-0008, whose Version-2 revision log explicitly
  left this file *active* alongside `labor_closures.md`)
- **Related:** ADR-0003, ADR-0008, DE-0006, `docs/DOCS_ASSESSMENT.md`,
  `paper/narrative_outline.md`, `docs/ideas/IDEA-0003-sectoral-sobol.md`

## Context

`docs/definitive_guide.md` (2026-09-03, no version stamp) was written between
the document review and Phase 0 of the registry, before any of the executed
generations existed. It was the strategy document of that moment and, per the
Version-2 entry of `docs/DOCS_ASSESSMENT.md`, it was deliberately kept active
when documents 1–8 moved to `docs/archive/`.

Three things made it stale in a way that matters for drafting:

1. **It presents retired results as live findings.** The `+19.3 pp` sweep is
   the retired unfinanced-shock mechanism (every shock is now financed,
   F1/F2/F3); the `88.4 %` / `100 %` variance shares and the
   `eta to infinity` reading are the rejected mobile-labour artefacts
   (DE-0006, `roadmaps/vertdict.md`); the Sobol ranking
   (theta 0.3951, sigma 0.2734, epsilon 0.1650, eta 0.1571) comes from a
   pre-registry run on code that has since changed (ADR-0012/0013
   calibration, ADR-0018 measurement, the C1 fix).
2. **It carries no version stamp**, so the revision scheme cannot show which
   of its claims were superseded.
3. **It is cited as evidence from the manuscript layer**:
   `paper/to_evaluate/closures_and_demand_shocks.md` asserted the wage-regime
   claim "as shown in `docs/definitive_guide.md` and our model runs".

## Decision

Move the file to `docs/archive/definitive_guide.md` and treat it as
superseded history. It is no longer an active reference:
`docs/DOCS_ASSESSMENT.md` drops it from its active-reference list (Version 9),
and no manuscript-facing document may cite it as evidence.

The move is compatible with the freeze register: `registry/freeze.toml` freezes
`docs/archive/` *per file* (the three legacy notebooks,
`selective_status_overview.md`, `varianten.xlsx`), so adding a file to the
directory changes no recorded tree hash.

## What was salvaged

The retirement was preceded by a claim-by-claim read of all 395 lines, because
a retired document can still hold a rule nobody else has written down.

**Carried forward (1 item):**

- The prohibition on **Type-I / Type-II multiplier framing that conflates
  intersectoral mobility with the wage regime** was not in the outline's
  must-not-appear import list; it now is — added to `paper/narrative_outline.md`,
  "What the revision retires".

**Checked and found already covered (no action):**

| Where it lived in the guide | Where it now lives |
| --- | --- |
| The "never say" list (88.4 %, 100 %, GO certification, eta to infinity, price invariance) | `paper/narrative_outline.md` "What the revision retires" and the external-review import bullet |
| "The original submission was not an equilibrium" (household overspend, dropped zero-profit) | `ROADMAP.md` §3 audit outcome, `roadmaps/vertdict.md`, DE-0001/DE-0002 |
| The 5.387 % production-vs-expenditure residual and `sum(lambda) = 2.1099` | the open-item list in `AGENTS.md` / `registry/closures.toml` (still open) |
| The Phase 1 accounting deliverables (`docs/accounting_consistency_plan.md`, `output/AC_*.csv`) | `docs/archive/accounting_consistency_plan.md` |
| The manuscript restructure (guide Part III) | `paper/narrative_outline.md` (v3) governs; the guide's nine sections are the retired thesis |
| The Phase 4 compliance tables (guide Part IV) | `ROADMAP.md` Phase 4, refreshed 2026-10-04 against the test suite |

**Not carried (retired with the document):**

- The Sobol variance decomposition numbers. The question itself is parked, not
  closed: `docs/ideas/IDEA-0003-sectoral-sobol.md` and
  `IDEA-0004-sectoral-effects-without-sobol.md` (both `scoped`) hold it, and
  DE-0006 forbids renormalising shares over main effects.
- The `+19.3 pp` / `fix_bridge` percentage-point gap tables (unfinanced shock).
- The `rerun_results.jl` headline table as an evidence source; the script and
  its tests remain in the tree, but the numbers are not part of the record.

## Consequences

- `docs/definitive_guide.md` no longer exists; the path moved by this ADR.
- Inbound references updated in the same pass:
  `docs/DOCS_ASSESSMENT.md` (active-reference list), `docs/ideas/README.md`
  (IDEA-0003 anchors), `docs/ideas/IDEA-0003-sectoral-sobol.md`,
  `paper/to_evaluate/closures_and_demand_shocks.md` (evidence re-pointed to the
  executed generation and `paper/to_evaluate/equivalence.tex`).
- The historical mentions inside `docs/DOCS_ASSESSMENT.md`'s inventory rows
  ("Mentioned in old docs") are left as written: they record where a claim was
  once stated, and the archive keeps the target readable.
- Anything wanting to revive a piece of it must go through the normal route:
  an idea record, or evidence from an executed generation (ADR-0004).
