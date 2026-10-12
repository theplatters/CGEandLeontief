# BeyondHulten — paper and reproduction materials

Paper and reproduction materials for the green-public-investment CGE study: a
production-network model (Baqaee–Farhi lineage) used to compare labour-market
and financing closures, with the German housing transformation as the
application (71-sector IO calibration). The manuscript sources are `paper/`
and `revision/`; everything the paper reports is backed by run manifests under
`runs/` and by the pre-split reproduction baseline in `baseline/`.

Status — closures, scenarios, runs, frozen snapshots — lives in `registry/`
and the generated board `docs/status.md`; do not restate it elsewhere
(ADR-0003). Start with `docs/status.md`.

## Relationship to the OpenCGE kernel

The model itself is **not** in this repository. It is the separate
[`OpenCGE`](https://github.com/theplatters/opencge) Julia package (local
checkout `../opencge`; the GitHub repository may not be pushed yet), extracted
from this repo as tag `v0.1.0` with its git history preserved (ADR-0028). This
repo depends on it via the committed root `Manifest.toml`, which pins the
kernel exactly (`repo-rev = "v0.1.0"`, `git-tree-sha1`); every run manifest
records that pin in its `[provenance]` block. Kernel implementation changes
happen in the OpenCGE repository and land here only via a Manifest pin bump —
which starts a new run generation.

Kernel-side decisions are summarized in OpenCGE's `DECISIONS.md`; the
paper-side decision records are `docs/decisions/` (ADR-0028 split the two).

## Quick start

[Julia](https://julialang.org/) ≥ 1.9 is required. From the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # install dependencies (incl. the pinned OpenCGE)
julia --project=. tests/runtests.jl                    # paper/reproduction test suite
julia --project=. scripts/reproduce.jl --check         # reproduce baseline/ within 1e-8
```

The model reads the formatted German 2019 IO table from
`data/I-O_DE2019_formatiert.csv`; `data/` is gitignored, so the file must be
placed there locally (the frozen `cbase2/data_raw/` snapshot holds a copy).

The revision manuscript builds with `latexmk -pdf` in `revision/` (entry
`revision/main.tex`). Jupyter Lab (https://jupyter.org/) is only needed to run
the legacy notebooks, which are frozen and read-only (see below). Plotting is
optional and goes through the `OpenCGE` GLMakie extension
(`import Pkg; Pkg.add("GLMakie")`, then `using GLMakie`).

## Commands

| Purpose | Command |
| --- | --- |
| Install dependencies | `julia --project=. -e 'using Pkg; Pkg.instantiate()'` |
| Kernel load check | `julia --project=. -e 'using OpenCGE'` |
| Test suite | `julia --project=. tests/runtests.jl` |
| Reproduction baseline | `julia --project=. scripts/reproduce.jl --check` (1e-8 gate against `baseline/`) |
| Regenerate status board | `julia --project=. scripts/status.jl` |
| CI gate (stale board) | `julia --project=. scripts/status.jl --check` |
| Repository gate | `julia --project=. scripts/check_repo.jl` (0 violations required) |
| List design cells | `julia --project=. experiments/run.jl --list <design>` |
| Preregister a design | `julia --project=. experiments/run.jl --preregister <design> [--actor NAME]` (commit the record first) |
| Run a design | `julia --project=. experiments/run.jl --design <design> [--cell <run_id>] [--cells a,b,c] [--runs-dir DIR] [--budget-seconds N] [--actor NAME]` |
| Build the revision | `latexmk -pdf` in `revision/` |

Every experiment goes through `experiments/run.jl`; the workflow is register
the scenario row → preregister the design → run (the entry point refuses on a
preregistration mismatch before creating any run directory). See
`experiments/README.md` for the design/manifest/index schemas and
`registry/README.md` for the tracking schemas.

## Repository layout

| Path | What it is |
| --- | --- |
| `OpenCGE` (external, `../opencge`) | The model kernel — model, closures, diagnostics. Pinned in `Manifest.toml` at tag `v0.1.0`. |
| `tests/` | Paper/reproduction test suite (entry `tests/runtests.jl`). |
| `baseline/` | Pre-split reproduction baseline: the metrics, output CSVs and gate logs the paper pipeline produced before the kernel moved out; `scripts/reproduce.jl --check` re-derives them within 1e-8. |
| `Manifest.toml` | Committed pin of the `OpenCGE` dependency — deliberately tracked. |
| `scripts/` | `status.jl` (generates `docs/status.md`), `check_repo.jl` (repository gate, ADR-0007), `reproduce.jl` (reproduction baseline tool). |
| `registry/` | Machine-readable single source of truth: `closures.toml`, `scenarios.csv`, `freeze.toml`, `preregistration.toml`. Kernel paths/ids are `OpenCGE@…`-qualified. |
| `experiments/` | Single run entry point `run.jl` and pinned designs in `experiments/designs/`. |
| `runs/` | One `runs/<run_id>/manifest.toml` + `log.txt` (`solution.csv` on success) per run; `runs/index.csv` is the committed index. |
| `docs/` | `status.md` (generated board), `decisions/` (ADRs), `dead-ends/` (DE register), `ideas/` (idea register), `log/` (lab log), `closures/` (closure notes), `archive/`. |
| `data/` | Input data (mostly gitignored; the IO table has to be placed here). |
| `paper/`, `revision/` | Manuscript sources (`revision/main.tex` is the revision build). |
| `cbase2/`, `bf_replication/`, `bf_replication2/`, `Replication Files/`, `Dokumente/`, `Notebooks/`, `archive/` | Frozen/read-only snapshots and archives (ADR-0001, ADR-0009); never edited in place. |

Key documents: `docs/status.md` (human entry point, generated from
`registry/`), `ROADMAP.md` (research plan and validation gates),
`docs/DOCS_ASSESSMENT.md` (documentation audit and closure catalogue),
`AGENTS.md` (operating contract), `docs/decisions/` (why things are the way
they are), `docs/dead-ends/` (abandoned approaches), `docs/ideas/` (recorded
but untested routes), `docs/log/` (narrative history), `baseline/README.md`
(the reproduction contract).

## Legacy material

The root notebooks (`DemandShocks.ipynb`, `CompareModels.ipynb`,
`CobbDouglas.ipynb`) were moved out of the root into
`docs/archive/notebooks/` and stay frozen/read-only (ADR-0009). Together
with the `Notebooks/` directory (`Analysis.ipynb`, `Translation.ipynb`, …)
they are the pre-package history: they predate the current API, are
committed with outputs, and are kept for reference only. The current model
lives in the `OpenCGE` kernel package. The root `rerun_results.jl` is a
frozen legacy pipeline kept for the provenance of
`output/variance_decomposition_sobol.csv` (see `baseline/README.md`).
