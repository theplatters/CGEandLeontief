# ADR-0009 — Archive the legacy root notebooks

- **Status:** accepted
- **Date:** 2026-09-17
- **Supersedes:** —
- **Related:** ADR-0001 (amended: the root notebooks leave their path),
  ADR-0008, `registry/freeze.toml`, `docs/archive/README.md`

## Context

The three notebooks `DemandShocks.ipynb`, `CompareModels.ipynb` and
`CobbDouglas.ipynb` are pre-package history: they predate the current `src/`
API, are committed with outputs, and are registered read-only in
`registry/freeze.toml` as `legacy-notebooks`. ADR-0001 kept all legacy
material at its path, so the three files still sat in the repository root,
mixed in with the live entry points (`README.md`, `ROADMAP.md`, `AGENTS.md`),
the reproduction script and the manifest sources.

Nothing in the repository references them besides historical records: the
kernel, tests, scripts and experiment tooling do not load them, and the
artefacts they produced were superseded by `rerun_results.jl`, the test
suite and the frozen `cbase2` snapshot. They are archive material, not live
root material — the same classification ADR-0008 applied to superseded root
documents.

## Decision

- The three root notebooks are moved to `docs/archive/notebooks/` (not
  deleted), next to the archived root documents of ADR-0008.
- Each notebook keeps its read-only registration in `registry/freeze.toml`,
  with the path and the archival commit updated and the blob hash unchanged;
  `scripts/status.jl` renders them on the freeze board as before.
- `docs/archive/README.md` lists the arrivals; git history preserves the
  original root paths.
- The `Notebooks/` directory and the other legacy zones named in ADR-0001
  stay at their paths and remain read-only. Only the "stays at its path"
  clause of ADR-0001 is amended; the one-kernel rule and the rest of
  ADR-0001 stand. ADR-0001 carries the amendment note.

## Consequences

- The repository root holds only live documents and configuration.
- The notebooks remain readable and citable; they still cannot be rerun
  against the current API, as already recorded in
  `docs/archive/REVISED_SALVAGE_PLAN.md`.
- Moving or editing an archived notebook again requires a new ADR and a new
  freeze record (ADR-0001).

## Enforcement

- `registry/freeze.toml` entries
  `docs/archive/notebooks/DemandShocks.ipynb`,
  `docs/archive/notebooks/CompareModels.ipynb` and
  `docs/archive/notebooks/CobbDouglas.ipynb`, each with its recorded
  archival commit and unchanged blob hash.
- `scripts/status.jl` warns when a notebook differs from its freeze record;
  `scripts/check_repo.jl` must stay at 0 violations.
- `docs/archive/README.md` lists the closed zone.
