# AGENTS.md — operating contract

## What this repository is

A production-network model of green public investment (BeyondHulten, Baqaee–Farhi
lineage) used to compare labour-market and financing closures. There is exactly
one canonical kernel: the root `src/` package (ADR-0001). Research plan and
validation gates live in `ROADMAP.md`; the documentation audit and closure
catalogue in `docs/DOCS_ASSESSMENT.md`.

## Layout (Phase 1 target)

| Path | What it is |
| --- | --- |
| `src/` | Canonical Julia package — model, closures, diagnostics. One kernel only: `src/core/{accounting,technology,equilibrium,diagnostics}.jl`, closure plugins in `src/closures/{labor,financing}/` plus `src/closures/registry.jl` (ADR-0005). |
| `ext/` | Phase 1: package extensions for heavy optional features (GLMakie plotting); they load only when the optional package is loaded. Wired via `[weakdeps]` / `[extensions]` in `Project.toml`. |
| `tests/` | Test suite (entry `tests/runtests.jl`, shimmed by `test/runtests.jl`). |
| `scripts/` | Repository tooling — `scripts/status.jl` generates the status board, `scripts/check_repo.jl` is the pre/post-batch gate (ADR-0007). |
| `registry/` | Machine-readable single source of truth: `closures.toml`, `scenarios.csv`, `freeze.toml`, `preregistration.toml`. Schema: `registry/README.md`. |
| `docs/` | `status.md` (generated status board), `decisions/` (ADRs), `dead-ends/` (DE register), `log/` (lab log), `archive/` (closed historical docs). |
| `experiments/` | Single run entry point `run.jl` plus pinned designs in `experiments/designs/` (schemas: `experiments/README.md`; ADR-0006). |
| `runs/` | One `runs/<run_id>/manifest.toml` + `log.txt` (`solution.csv` on success) per run; committed `runs/index.csv` holds one row per run (ADR-0004). Only the index and manifests are tracked. |
| `archive/src-orphans/` | Read-only archive of abandoned sources that were never wired into the module. |
| `.opencode/skills/` | Phase 1: agent workflow skills `tracking`, `closures`, `experiments`. |
| `data/` | Input data; mostly gitignored, a few small reference files are tracked. |
| `paper/`, `revised_manuscript/` | Manuscript sources. |

Frozen/read-only — never edit (ADR-0001; commits and tree hashes in
`registry/freeze.toml`): `cbase2/`, `bf_replication/`, `bf_replication2/`,
`Replication Files/`, `Dokumente/`, `Notebooks/`, the legacy root notebooks
(`DemandShocks.ipynb`, `CompareModels.ipynb`, `CobbDouglas.ipynb`), and `archive/`.

## Commands

From the repository root:

| Purpose | Command |
| --- | --- |
| Install dependencies | `julia --project=. -e 'using Pkg; Pkg.instantiate()'` |
| Load check | `julia --project=. -e 'using BeyondHulten'` |
| Test suite | `julia --project=. -e 'using Pkg; Pkg.test()'` |
| Regenerate status board | `julia --project=. scripts/status.jl` |
| CI gate (stale board) | `julia --project=. scripts/status.jl --check` |
| Repository gate | `julia --project=. scripts/check_repo.jl` (0 violations required, before and after every batch) |
| List design cells | `julia --project=. experiments/run.jl --list <design>` |
| Preregister a design | `julia --project=. experiments/run.jl --preregister <design> [--actor NAME]` (commit the record first) |
| Run a design | `julia --project=. experiments/run.jl --design <design> [--cell <run_id>] [--cells a,b,c] [--runs-dir DIR] [--budget-seconds N] [--actor NAME]` |
| Optional plotting | `julia --project=. -e 'using Pkg; Pkg.add("GLMakie")'`, then `using GLMakie`; extended functionality loads lazily. |

Julia ≥ 1.9 is required. `Manifest.toml` is gitignored — do not commit it.

## Operating contract

### Before changing anything

1. Read `docs/status.md` (human entry point) and the relevant entries in
   `registry/` (`closures.toml`, `scenarios.csv`, `freeze.toml`).
2. Read the ADRs in `docs/decisions/` and the dead ends in `docs/dead-ends/`
   that touch the area.
3. Load the matching skill from `.opencode/skills/` when the task matches:
   - `tracking` — status, decisions, dead ends, log;
   - `closures` — closure formulations, promotion, registry entries;
   - `experiments` — runs, scenarios, manifests.
4. Frozen zones are read-only (ADR-0001). Never edit a frozen file in place,
   not even to "fix" it: a change inside one requires a superseding ADR.
   Frozen hashes are recorded in `registry/freeze.toml`.

### Status transitions require evidence

| Status | Evidence required |
| --- | --- |
| `implemented` | Code exists in the listed file and runs. |
| `tested` | The listed repo tests pass. |
| `validated` | `ROADMAP.md` Phase 4 gates pass, with a recorded run manifest. |
| `rejected` | A linked `docs/dead-ends/DE-*` record exists. |
| `superseded` | The replacing id is named. |

Never upgrade a status without the evidence. Failed and provisional runs stay
visible (ADR-0004).

### After changing something

- Update `registry/*`: `closures.toml` status/symbols/tests, `scenarios.csv`
  rows; `freeze.toml` only via ADR.
- Regenerate `docs/status.md` (`julia --project=. scripts/status.jl`); never
  hand-edit it, and it must show **0 warnings**.
- Add an ADR for a decision, a DE record for an abandoned approach, and one
  entry in `docs/log/YYYY-MM.md` for the session.
- Never restate status in prose (ADR-0003); other documents link to `registry/`
  or `docs/status.md`.

### Runs

- Every experiment is a scenario in `registry/scenarios.csv`; ids follow
  `<design>-<labor>-<financing>[-<variant>]` (ADR-0002).
- Workflow per batch: register the scenario row (`planned`) → preregister the
  design (`run.jl --preregister`, commit `registry/preregistration.toml`) →
  run via `experiments/run.jl --design` (the only entry point; it refuses on a
  preregistration mismatch before creating any run dir) → `run.jl` updates the
  scenario row and `runs/index.csv` on every status transition.
- Run `scripts/check_repo.jl` before and after every batch; it must report
  0 violations. It enforces the board/warnings gate, registry bidirectionality,
  manifest/index/scenario consistency, preregistration integrity, and the WIP limit.
- WIP limit: at most one scenario row may be `running` at a time (ADR-0007).
- Failed and provisional runs stay visible with their caveats; never drop them.
- Paper tables and figures cite `run_id`s (ADR-0004).
- Current state: the `matrix_5x3` batch is blocked at the real-data reference
  continuation (first rung stalls at resid 3.6e-4) — a recorded model/solver
  open item (freeze `open_items`, closure `open_gates`), not a tooling gap.
  The batch refuses loudly instead of warm-starting from a bad root.

### Commits

- `closure(<id>): …` — closure formulation or implementation change
- `exp(<design>): run <id>` — executed run
- `adr(NNNN): …`, `deadend(NNNN): …`, `log(YYYY-MM): …` — tracking records
- Otherwise short imperative subjects. Keep commits focused; never mix
  generated figures with model changes.

## Code and tests

- Match the style of surrounding code.
- Tests live in `tests/`; new files are named `test_<feature>.jl` and included
  from `tests/runtests.jl`.
- Use `@testset` and `isapprox` tolerances, not exact floating-point equality.
- Experimental code lives outside `src/` (scratch, notebooks, `experiments/`
  smoke designs). It enters `src/` only promoted: registry entry + tests +
  docs, or it is archived with a DE record. No zombie files in `src/` — every
  file must be reachable from `src/BeyondHulten.jl` (`check_repo.jl` enforces this).
- Run `Pkg.test()` before finishing work.

## House rules

- No plan files in the repository. Plans live in `ROADMAP.md`, ADRs, or the
  skill/registry workflow.
- Do not commit proprietary or machine-local data or generated artifacts;
  `data/` and generated figures are largely gitignored.
