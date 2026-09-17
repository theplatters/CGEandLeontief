# DE-0009 — Orphan `src/` files archived, never wired into the package module

- **Status:** dead end — never wired into the package module
- **Recorded:** 2026-09-17
- **Origin:** `src/` at commit `b05908d`; moved to `archive/src-orphans/` in Phase 1
- **Related:** ADR-0001, closure `BF` (`registry/closures.toml`), `docs/archive/REVISED_SALVAGE_PLAN.md` (historical mention)
- **Evidence:** `src/BeyondHulten.jl:41-50` (include list contains none of the six files), `git status` renames `R src/<file> -> archive/src-orphans/<file>`, reference greps below

## What was tried

Six files accumulated in `src/` alongside the canonical kernel, each a
plausible direction on its own:

- `ces_temporal.jl` — adds an intertemporal elasticity dimension
  (`CESTemporalElasticities.intertemporal_elasticity`) to the CES model.
- `goodwin.jl` — a Goodwin-cycle model (ModelingToolkit/OrdinaryDiffEq) with
  sectoral wages, profits, and employment dynamics.
- `loaders.jl` — MAT-file loaders (`loadInData`, `loadStfp`).
- `main.jl` — a one-off driver script running panels, elasticity sweeps, and
  writing `data/sector_names.txt` (`src/main.jl:41`).
- `mvnrnd.jl` — a 4-line multivariate-normal helper via Cholesky.
- `mobile_labor-JK.jl` — a 592-line alternative mobile-labour implementation
  exposing the same `MobileLaborCES` API as the canonical file (the "JK"
  suffix marks it as an author variant).

## Why it seemed promising

Each file looked like it could belong: a temporal CES extension, a dynamic
Goodwin model, data loaders, a runnable entry point, a sampling helper, and
a fuller mobile-labour implementation.

## How it failed

None of the six was ever `include`d by `src/BeyondHulten.jl` (verified:
the module includes exactly ten files, none of them these). Reference greps
over `*.jl`, `*.ipynb`, `*.md` found no live caller for any of them:

- `CESTemporal*`, `create_economic_model`, `loadInData`/`loadStfp` are
  referenced only inside their own files.
- The `ces_temporal.jl` residual stub references undefined symbols (`k`,
  `theta`, `A`, `Ind`) and could never have run.
- `goodwin.jl` executes solve-and-plot statements at top level — a script,
  not a module component — and pulls heavy dependencies (`ModelingToolkit`,
  `OrdinaryDiffEq`, `Plots`) for a different model class.
- `main.jl`'s only observable side effect, `data/sector_names.txt`, has no
  reader (the `cbase2` notebook reads the unrelated
  `data_raw/sector_names.txt`).
- `mvnrnd.jl`'s only consumers are the frozen legacy notebooks
  `Notebooks/Analysis.ipynb` and `Notebooks/Translation.ipynb` via
  `include("src/mvnrnd.jl")`; those notebooks are read-only and were not
  edited, so the reference now points at the archive path.
- `mobile_labor-JK.jl` is an older, less-hardened sibling of
  `src/mobile_labor.jl`: a diff shows the canonical file adds `_checked_eta`,
  `_positive_floor`, clamped log-space demand/interpolation helpers, and a
  floored efficiency wedge that the JK variant lacks.

## What replaced it

Nothing per-file was needed: the canonical kernel (`src/`) plus the
registered `BF` closure (`registry/closures.toml`) already cover the live
programme. CSV/XLSX pipelines replaced the MAT loaders; the mobile-labour
closure replaced both labour implementations' shared idea.

## Revival conditions

A file may return only as a registered closure/model with tests, an ADR, and
a registry entry — never as parallel `src/` code (ADR-0001). In particular,
`mobile_labor-JK.jl` must not be revived as a second kernel.
