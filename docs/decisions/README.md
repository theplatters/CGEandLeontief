# Decision records (ADRs)

Append-only record of design decisions about the model, the data pipeline,
and the way this repository is organized. An ADR explains **why** something
is the way it is, so that neither humans nor agents have to re-litigate a
settled question.

## Rules

- One file per decision: `ADR-####-short-slug.md`, numbered consecutively,
  never reused.
- A decision is immutable once `accepted`. To change it, write a new ADR with
  `Supersedes: ADR-####` and change the old record's status to
  `superseded by ADR-####` — do not rewrite its content.
- Status values: `proposed`, `accepted`, `superseded`, `rejected`.
- Cite related dead ends (`DE-####`), registry entries, and code paths.
- A commit that changes a closure formulation, a status rule, or the freeze
  policy must add or supersede an ADR (see `AGENTS.md`).

## Template

```markdown
# ADR-#### — Title

- **Status:** accepted
- **Date:** YYYY-MM-DD
- **Supersedes:** —
- **Related:** ADR-####, DE-####, registry entries

## Context
## Decision
## Consequences
## Enforcement
```

## Index

| ID | Title | Status |
| --- | --- | --- |
| ADR-0002 | Closure taxonomy and stable IDs | accepted |
| ADR-0003 | `registry/` is the single source of truth | accepted |
| ADR-0004 | Runs are immutable, manifest-backed records | accepted |
| ADR-0006 | Experiment entry point and preregistration | accepted |
| ADR-0007 | Repository gate and WIP limits | accepted |
| ADR-0008 | Archive superseded root documents as read-only | accepted |
| ADR-0009 | Archive the legacy root notebooks | accepted |
| ADR-0012 | A-bill calibration: domestic intermediate bill, clamp eliminated, s re-anchored | accepted |
| ADR-0013 | the intermediate-bill tax term: row 75 booked as an external leak | accepted |
| ADR-0018 | Real GDP measurement in the open economy | accepted |
| ADR-0023 | Supply-shock design schema for the identification arm | accepted |
| ADR-0024 | Programme total: G0 = 40,300 kept (the 2019-price horizon average) | accepted |
| ADR-0025 | Matrix visualisation: validated re-solve, headless data layer, figures are not runs | accepted |
| ADR-0026 | Retire `docs/definitive_guide.md`: moved to `docs/archive/`, superseded results, salvage review recorded | accepted |
| ADR-0027 | Consolidate operating docs: merge/archive MERITOF, WORKPLAN, labor_closures; one plan (ROADMAP) | accepted |
| ADR-0028 | Split into paper repo and OpenCGE kernel repo | accepted |
| ADR-0029 | The pricing door: an output-dependent unit cost (Foundation III) | accepted |
| ADR-0030 | Distributed profit income: Stage B of the markup door | accepted |
| ADR-0031 | Value-added-scaled capacity door (OpenCGE v0.2.1) and the services-zero δ profile | accepted |
| ADR-0032 | The pricing doors on every labour row (DELTA wired in the harness, sectoral BETA) | accepted |
| ADR-0033 | Financing formulations name the balancing variable (wage-rate tax, endogenous F, B_gov) | proposed |
| ADR-0036 | Print-ready matrix panels (OpenCGE v0.3.2) without a new generation | accepted |
| ADR-0037 | Consumption-unit price reporting | accepted |

**Migrated to OpenCGE DECISIONS.md** (2026-10-12, ADR-0028): ADR-0001,
ADR-0005, ADR-0010, ADR-0011, ADR-0014, ADR-0015, ADR-0016, ADR-0017,
ADR-0019, ADR-0020, ADR-0021 and ADR-0022 were kernel-layout and
model-kernel decisions. Their one-to-two-sentence summaries live in the
OpenCGE kernel's `DECISIONS.md` (`../opencge`, referenced as `OpenCGE@ADR-…`
from `registry/`); the full texts remain recoverable from this repository's
git history. The ids are never reused or renumbered (ADR-0021 is `proposed`
there, as it was here).
