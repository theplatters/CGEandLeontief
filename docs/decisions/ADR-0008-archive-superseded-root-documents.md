# ADR-0008 — Archive superseded root documents as read-only

- **Status:** accepted
- **Date:** 2026-09-17
- **Supersedes:** —
- **Related:** ADR-0001, ADR-0003, `registry/freeze.toml`,
  `docs/archive/README.md`, DE-0001, DE-0002, DE-0003, DE-0006

## Context

The repository root still carried two pre-Phase-0 documents that contradicted
the current record:

- `selective_status_overview.md`, an old status snapshot whose mobile-labour
  "GO" and "88.4%" variance share are rejected by `ROADMAP.md` §3,
  `roadmaps/vertdict.md`, and `docs/dead-ends/DE-0001`, `DE-0002`, `DE-0003`
  and `DE-0006`. Phase 0 marked it superseded with a banner but left it in the
  root, where it still read like a current status document.
- `varianten` (now `varianten.xlsx`), a January 2025 Excel workbook sketching
  labour-closure variants against the 2019 and 2022 Baqaee–Farhi models.
  Nothing in the repository references it, and the design space it sketches is
  recorded in `docs/DOCS_ASSESSMENT.md` and `registry/closures.toml`.

`docs/archive/` was already declared a closed zone in its README, but it was
invisible to `registry/freeze.toml` and `scripts/check_repo.jl`: files moved
there were protected by convention only, and the root kept carrying status
claims that the registry contradicts.

## Decision

- Superseded documents are moved, not deleted, from the repository root into
  `docs/archive/`, with a banner in the document naming what superseded it and
  the decision that archived it.
- Each archived file is registered as a read-only single-file zone in
  `registry/freeze.toml` with the archival commit and the blob hash, so
  `scripts/status.jl` renders it on the freeze board and `scripts/check_repo.jl`
  flags any later modification.
- `docs/archive/README.md` lists the arrivals; the original path stays
  reachable through git history.
- `archive/` keeps the orphan-source archive (`archive/src-orphans/`, DE-0009);
  deprecated documents go to `docs/archive/`.
- Later archival of a stale document follows this decision and records the
  corresponding read-only freeze entry; it does not require a new ADR.

## Consequences

- The root no longer contains documents that contradict `docs/status.md`.
- Archived files remain readable and citable; corrections go to the current
  document, never into the archive.
- Editing an archived file, unfreezing it, or moving it requires a superseding
  ADR and a new freeze record (ADR-0001).
- The repository gate skips registered read-only files in its artifact scan,
  exactly as it does for the legacy notebooks.

## Enforcement

- `registry/freeze.toml` entries `docs/archive/selective_status_overview.md`
  and `docs/archive/varianten.xlsx`, each with its recorded archival commit
  and blob hash.
- `scripts/status.jl` warns when an archived file differs from its freeze
  record; `scripts/check_repo.jl` must stay at 0 violations.
- `docs/archive/README.md` lists the closed zone.
