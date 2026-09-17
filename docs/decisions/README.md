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
| ADR-0001 | One kernel, no living copies | accepted |
| ADR-0002 | Closure taxonomy and stable IDs | accepted |
| ADR-0003 | `registry/` is the single source of truth | accepted |
| ADR-0004 | Runs are immutable, manifest-backed records | accepted |
| ADR-0005 | Phase 2 layout and cbase2 backport | accepted |
| ADR-0006 | Experiment entry point and preregistration | accepted |
| ADR-0007 | Repository gate and WIP limits | accepted |
| ADR-0008 | Archive superseded root documents as read-only | accepted |
| ADR-0009 | Archive the legacy root notebooks | accepted |
| ADR-0010 | Port the `cbase2/review.md` pipeline fixes | accepted |
