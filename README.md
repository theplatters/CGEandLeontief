# BeyondHulten

A production-network model of green public investment (Baqaee–Farhi lineage)
used to compare labour-market and financing closures, with the German housing
transformation as the application (71-sector IO calibration). The current
model is the canonical Julia package in `src/`.

Status — closures, scenarios, runs, frozen snapshots — lives in `registry/`
and the generated board `docs/status.md`; do not restate it elsewhere
(ADR-0003). Start with `docs/status.md`.

## Installation

[Julia](https://julialang.org/) ≥ 1.9 is required.

From the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # install dependencies (Project.toml)
julia --project=. -e 'using BeyondHulten'             # load check
julia --project=. -e 'using Pkg; Pkg.test()'          # test suite
```

The model reads the formatted German 2019 IO table from
`data/I-O_DE2019_formatiert.csv`; `data/` is gitignored, so the file must be
placed there locally (the frozen `cbase2/data_raw/` snapshot holds a copy).
`Manifest.toml` is machine-local and not committed.

Jupyter Lab (https://jupyter.org/) is only needed to run the legacy
notebooks, which are frozen and read-only (see below).

Plotting is optional: the GLMakie-based figures live in a package extension
(`ext/`) and load only when `GLMakie` is available
(`import Pkg; Pkg.add("GLMakie")`, then `using GLMakie, BeyondHulten`). The
package itself installs, loads, and tests without the heavy graphics stack.

## Commands

| Purpose | Command |
| --- | --- |
| Install dependencies | `julia --project=. -e 'using Pkg; Pkg.instantiate()'` |
| Load check | `julia --project=. -e 'using BeyondHulten'` |
| Test suite | `julia --project=. -e 'using Pkg; Pkg.test()'` |
| Regenerate status board | `julia --project=. scripts/status.jl` |
| CI gate (stale board) | `julia --project=. scripts/status.jl --check` |
| Repository gate | `julia --project=. scripts/check_repo.jl` (0 violations required) |
| List design cells | `julia --project=. experiments/run.jl --list <design>` |
| Preregister a design | `julia --project=. experiments/run.jl --preregister <design> [--actor NAME]` (commit the record first) |
| Run a design | `julia --project=. experiments/run.jl --design <design> [--cell <run_id>] [--cells a,b,c] [--runs-dir DIR] [--budget-seconds N] [--actor NAME]` |

Every experiment goes through `experiments/run.jl`; the workflow is register
the scenario row → preregister the design → run (the entry point refuses on a
preregistration mismatch before creating any run directory). See
`experiments/README.md` for the design/manifest/index schemas and
`registry/README.md` for the tracking schemas.

## Repository layout

| Path | What it is |
| --- | --- |
| `src/` | Canonical kernel — model, closures (`src/core/`, `src/closures/`), diagnostics. One kernel only (ADR-0001). |
| `ext/` | Package extension for optional plotting; loads only with `GLMakie`. |
| `tests/` | Test suite (entry `tests/runtests.jl`, shimmed by `test/runtests.jl`). |
| `scripts/` | `status.jl` (generates `docs/status.md`), `check_repo.jl` (repository gate, ADR-0007). |
| `registry/` | Machine-readable single source of truth: `closures.toml`, `scenarios.csv`, `freeze.toml`, `preregistration.toml`. |
| `experiments/` | Single run entry point `run.jl` and pinned designs in `experiments/designs/`. |
| `runs/` | One `runs/<run_id>/manifest.toml` + `log.txt` (`solution.csv` on success) per run; `runs/index.csv` is the committed index. |
| `docs/` | `status.md` (generated board), `decisions/` (ADRs), `dead-ends/` (DE register), `log/` (lab log), `closures/` (closure notes), `archive/`. |
| `data/` | Input data (mostly gitignored; the IO table has to be placed here). |
| `paper/`, `revised_manuscript/` | Manuscript sources. |
| `cbase2/`, `bf_replication/`, `bf_replication2/`, `Replication Files/`, `Dokumente/`, `Notebooks/`, `archive/` | Frozen/read-only snapshots and archives (ADR-0001); never edited in place. |

Key documents: `docs/status.md` (human entry point, generated from
`registry/`), `ROADMAP.md` (research plan and validation gates),
`docs/DOCS_ASSESSMENT.md` (documentation audit and closure catalogue),
`AGENTS.md` (operating contract), `docs/decisions/` (why things are the way
they are), `docs/dead-ends/` (abandoned approaches), `docs/log/` (narrative
history).

## Legacy material

The root notebooks (`DemandShocks.ipynb`, `CompareModels.ipynb`,
`CobbDouglas.ipynb`) and the `Notebooks/` directory (`Analysis.ipynb`,
`Translation.ipynb`, …) are the pre-package history: they predate the current
API, are committed with outputs, and are frozen/read-only for reference only.
The current model lives in `src/`.
