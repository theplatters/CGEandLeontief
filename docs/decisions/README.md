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
| ADR-0011 | Conform the frozen cbase2 notebook 01 to the canonical kernel | accepted |
| ADR-0012 | A-bill calibration: domestic intermediate bill, clamp eliminated, s re-anchored | accepted |
| ADR-0013 | the intermediate-bill tax term: row 75 booked as an external leak | accepted |
| ADR-0014 | real-wage elastic supply and a verified scale-determinacy guard | accepted |
| ADR-0015 | monotone polish to ~1e-10 (tolerance-independent metrics) | accepted |
| ADR-0016 | CES-consistent valuation of the intermediate-bill leaks | accepted |
| ADR-0017 | Actual-matrix scale determinacy for the fixed-wage η = 1 system | accepted |
| ADR-0018 | Real GDP measurement in the open economy | accepted |
| ADR-0019 | All-N goods-market clearing with an explicit external account | accepted |
| ADR-0020 | The eta = 0 endpoint's external account: closure options | accepted (option C) |
| ADR-0021 | Wage structure in the fixed-wage closure (GAMMA) | proposed |
| ADR-0022 | Sectoral labour markets: demand-sensitive prices from the wage block | accepted |
