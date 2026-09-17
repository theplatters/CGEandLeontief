# ADR-0004 — Runs are immutable, manifest-backed records

- **Status:** accepted
- **Date:** 2026-09-17
- **Related:** ADR-0003, `registry/scenarios.csv`, `ROADMAP.md` §8

## Context

Results currently live in notebooks, CSVs, plots, and prose, and solver
failures have been discarded silently (for example, the former variance
decomposition converted an incomplete design into a nonorthogonal one and
renormalized shares over main effects, producing the invalid "88.4%"
headline — `roadmaps/vertdict.md`, `docs/dead-ends/DE-0006`). The `cbase2`
pipeline already introduced preregistration and SHA-256 input provenance,
but only inside that folder.

Without run-level records it is impossible to trace a paper number to the
code version, calibration vintage, and gates that produced it.

## Decision

- Every experiment is a named **scenario** in `registry/scenarios.csv` with a
  unique `run_id` of the form `<design>-<labor>-<financing>[-<variant>]`
  (ADR-0002). The Phase-0 seed registers the planned 5×3 matrix and the
  historical `cbase2-v3` verification runs.
- From Phase 3, scenarios are executed through one entry point
  (`experiments/run.jl`) that writes `runs/<run_id>/`:
  - `manifest.toml` — git commit and dirty flag, scenario/design hash,
    data-vintage SHA-256 values, Julia and package versions (Manifest hash),
    solver settings and seed, per-gate residual results, actor
    (human/agent), and pointers to artifacts;
  - the raw results, logs, and figures for that run.
- A committed `runs/index.csv` holds one row per run: `run_id`, date, design,
  closures, status, gate summary, headline metrics, commit. `runs/` itself is
  gitignored except the index and manifests.
- **Failed and provisional runs are recorded, never dropped.** If a required
  factorial cell or matrix cell fails its gates, it stays visible with its
  evidence; the analysis must either fix it or report the feasible domain.
- Scenario status transitions (`provisional` → `executed`, `failed`, …) are
  made in `registry/scenarios.csv`; paper tables and figures cite `run_id`s.
- Historical runs that predate manifests (the `cbase2-v3` rows) are
  registered at their recorded commit with their caveats; they are not
  retrofitted with invented manifests.

## Consequences

- Every quantitative claim is traceable: claim → `run_id` → manifest →
  commit → closure/design → data vintage.
- Missing or failed cells are visible on `docs/status.md` instead of being
  omitted from a summary.
- Running an experiment becomes a defined operation (Phase 3) rather than
  ad-hoc notebook execution; notebooks may orchestrate and display, but
  solves happen in scripts.
- Storage discipline: run artifacts are not committed; the index is.

## Enforcement

- `registry/scenarios.csv` validation and, from Phase 3,
  `scripts/check_repo.jl`: every `runs/` directory has a manifest, every
  manifest has an index row, and no `runs/` artifacts are tracked by git.
- `AGENTS.md` requires a manifest and registry row for every executed
  scenario, and a dead-end record for every abandoned approach.
