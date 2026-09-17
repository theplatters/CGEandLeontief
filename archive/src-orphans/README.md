# `archive/src-orphans/` — closed, read-only zone

These six files lived in `src/` but were **never included by the package
module**: `src/BeyondHulten.jl` includes exactly `interface.jl`,
`solution.jl`, `cobbdouglas.jl`, `leontief.jl`, `ces.jl`, `mobile_labor.jl`,
`util.jl`, `impulses.jl`, `plots.jl`, `variance_decomposition.jl` — none of
the files below. They were moved here with `git mv` on 2026-09-17 (Phase 1)
from commit `b05908d` so that `src/` contains only the canonical kernel.
History is preserved, not deleted.

**Do not add new code here and do not `include` these files from the
package.** To rerun one you must create your own throwaway environment and
add its old dependencies: `Plots`, `ModelingToolkit`, `OrdinaryDiffEq` for
`goodwin.jl`; `MAT` for `loaders.jl` (see `Project.toml` at `b05908d`).

| File | Former role | Why archived |
| --- | --- | --- |
| `ces_temporal.jl` | Intertemporal CES sketch (`CESTemporal`, `CESTemporalElasticities` with an `intertemporal_elasticity` field) | No caller; no model uses it. The `problem` stub references undefined symbols (`k`, `theta`, `A`, `Ind`) and never ran. |
| `goodwin.jl` | Goodwin-cycle model built on `ModelingToolkit`/`OrdinaryDiffEq`, with top-level solve-and-plot statements | A different model class (dynamic Goodwin cycles), never integrated into the static closure programme. |
| `loaders.jl` | `MAT`-file loaders (`loadInData`, `loadStfp` reading `data/simulationData.mat`, `data/stfp.mat`) | Superseded by the CSV/XLSX data pipelines (`read_data`, `read_data_cb`). |
| `main.jl` | One-off driver script (`(@main)` running panels, sweeps, and writing `data/sector_names.txt`) | A script, not package code; nothing includes or calls it. |
| `mvnrnd.jl` | 4-line multivariate-normal helper (`mvnrnd(μ, Σ)` via Cholesky) | Only consumers were two legacy notebooks via `include("src/mvnrnd.jl")` (see below). `bf_replication/src/shocks.jl` carries its own independent port of `mvnrnd.m` (frozen; untouched). |
| `mobile_labor-JK.jl` | 592-line alternative mobile-labour implementation (same `MobileLaborCES` API as the canonical file) | Superseded by `src/mobile_labor.jl` (which adds log-space hardening: `_checked_eta`, `_positive_floor`, clamped helpers) and the registered `BF` closure (`registry/closures.toml`). Must not be revived as a parallel kernel (ADR-0001). |

**Broken reference left in place deliberately:** the frozen legacy notebooks
`Notebooks/Analysis.ipynb` and `Notebooks/Translation.ipynb` (plus their
`.ipynb_checkpoints` copies) contain `include("src/mvnrnd.jl")`. Those
notebooks are read-only frozen zones and were not edited; anyone rerunning
them must repoint the include to `archive/src-orphans/mvnrnd.jl`.

See `docs/dead-ends/DE-0009-orphan-src-files.md` and
`docs/decisions/ADR-0001-one-kernel-no-living-copies.md`.
