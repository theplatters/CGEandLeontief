# ADR-0003 — `registry/` is the single source of truth

- **Status:** accepted
- **Date:** 2026-09-17
- **Related:** ADR-0001, ADR-0004, `registry/README.md`, `docs/status.md`

## Context

Current status was stated in at least eight documents that contradict each
other. Examples from the freeze-day audit:

- `selective_status_overview.md` reports mobile labour and the variance
  decomposition as complete with a "GO" and an "88.4%" headline;
- `roadmaps/vertdict.md` and `ROADMAP.md` declare the same mobile-labour
  results invalid until regenerated;
- `cbase2/documentation.md` lists `src/validation.jl`, `scripts/run_*.jl`,
  and notebooks 04–08 as existing or pending, while only notebooks 01–03
  exist;
- `bf_replication/REPLICATION_ASSESSMENT.md` and `REPLICATION_REPORT.md`
  give different oil-shock headlines.

Agents and humans could not tell which document was current, and dead ends
were re-litigated because they lived in prose.

## Decision

- `registry/closures.toml`, `registry/scenarios.csv`, and
  `registry/freeze.toml` are the **only** place that states what exists and
  in which state. Schema and vocabularies: `registry/README.md`.
- `docs/status.md` is **generated** from the registry by
  `scripts/status.jl` and must never be edited by hand. It is the human
  entry point; the registry is the machine one.
- Other documents **link** to the registry or `docs/status.md`; they do not
  restate status. Superseded plans get a one-line banner pointing to
  `docs/status.md` and are moved to `docs/archive/` — never silently
  rewritten.
- Status transitions carry evidence requirements:
  - `implemented` — code exists in a listed file and runs;
  - `tested` — listed tests cover the contract and the suite passes;
  - `validated` — the Phase 4 gates of `ROADMAP.md` pass and a recorded run
    manifest is linked;
  - `rejected` — a `docs/dead-ends/DE-*` record exists and is linked;
  - `superseded` — the replacing id is named.
- Abandoned approaches require a dead-end record (`docs/dead-ends/`), and
  design choices require an ADR (`docs/decisions/`); both are append-only.
- `scripts/status.jl --check` (and later `scripts/check_repo.jl`) fail when
  the board is stale or references are broken, so review is part of the
  workflow rather than an afterthought.

## Consequences

- "Which document is right?" becomes "what does `docs/status.md` say?".
- The board shows warnings for dangling references (missing files, unknown
  ids, run directories that do not exist), making drift visible instead of
  hidden.
- Updating status in prose is a process violation; the fix is always a
  registry edit (plus ADR/DE/log where required).
- Paper claims should cite run ids from `registry/scenarios.csv` (see
  ADR-0004) rather than restating numbers in documentation.

## Enforcement

- `registry/README.md` documents the schema; `scripts/status.jl` validates
  vocabularies and references.
- `AGENTS.md` instructs agents to read `docs/status.md`, the relevant ADRs
  and dead ends before acting, and to update the registry as part of any
  status transition.
