---
name: experiments
description: Use when defining, running, or recording an experiment or scenario in the BeyondHulten repo — run ids, registry/scenarios.csv rows and evidence paths, failed and provisional run retention, the ADR-0004/ADR-0006 run-manifest contract, the experiments/run.jl entry point and runs/index.csv, headless Julia runs that avoid GLMakie, and citing run ids in the paper.
---

# Experiments: scenarios and runs

## Scenario ids

`<design>-<labor>-<financing>[-<variant>]`, for example `matrix_5x3-BETA-F2` or `cbase2-v3-ALPHA-F1-mobile`. `labor`/`financing` are stable ids from `registry/closures.toml` (ADR-0002). `run_id` is stable and unique; never reuse or rename it.

## Register in registry/scenarios.csv

Header (ASCII only, `TBD` for unknown, double-quote fields that contain commas):

```
run_id,design,status,labor,financing,eta,eta_s,theta,epsilon,sigma,shock,magnitude,data_vintage,evidence,commit,actor,notes
```

Checklist:

1. Before running: append a row with `status=planned`, the parameter cells (`eta`, `eta_s`, `theta`, `epsilon`, `sigma`), `shock`, `magnitude`, `data_vintage`, `actor`, and the configuration plan in `notes`.
2. Run from the repo root in the project environment: `julia --project=.` or `julia --project=. <script>.jl`.
3. After running: update `status`, set `evidence` to `;`-separated repo-relative paths that exist, set `commit` to that recorded commit, and put caveats/residuals in `notes`.
4. Regenerate the board: `julia --project=. scripts/status.jl` → 0 warnings → `julia --project=. scripts/status.jl --check`.
5. Commit as `exp(<design>): run <run_id>` and log the session in `docs/log/YYYY-MM.md`.

Rules:

- Register before or with the run; update after. Evidence paths are validated by `scripts/status.jl` and must exist; `TBD` is allowed only for unknown cells.
- Failed and provisional runs stay visible — never delete a row and never drop a failed matrix/factorial cell silently (ADR-0004). `docs/status.md` renders all rows, including failures.
- Do not restate results in prose; link `docs/status.md` or the scenario row.

## Scenario status vocabulary

`planned | running | executed | provisional | failed | superseded | cancelled`

- `executed` — ran and passed the gates recorded with it.
- `provisional` — ran but the result is not gate-clean; the caveat is in `notes`.
- `failed` — ran and did not pass; evidence is kept.
- `superseded` — replaced by a later `run_id`.

## Run manifests (ADR-0004, ADR-0006)

Every quantitative claim is traceable: claim → `run_id` → manifest → commit → closure/design → data vintage.

- **Entry point `experiments/run.jl`** (register before running):
  `julia --project=. experiments/run.jl --list <design>` prints cells with pinned parameters;
  `julia --project=. experiments/run.jl --preregister <design> [--actor NAME]` writes `registry/preregistration.toml` (the only writer);
  `julia --project=. experiments/run.jl --design <design> [--cell <run_id>] [--cells a,b,c] [--runs-dir DIR] [--budget-seconds N] [--actor NAME]` executes cells in file order.
  `--design` refuses to start unless the design file's SHA-256 matches the preregistration record. Design files: `experiments/designs/<design>.toml` (schemas: `experiments/README.md`).
- Each run writes `runs/<run_id>/manifest.toml` (`running` → `executed`/`failed`; per-cell exceptions become `failed` with an `[error]` table), `log.txt`, and `solution.csv` (`sector,price,quantity`) on success. Existing run dirs are never overwritten — reruns register a `-v2` variant `run_id`. `runs/index.csv` rows: `run_id`, date, design, closures, status, gate summary, headline metrics, commit; one row per run, sorted by `run_id`, rewritten on status transitions. Only the index and manifests are tracked (`.gitignore`; verify with `git check-ignore`).
- Manifest contents: git commit and dirty flag, design hash, data-vintage SHA-256 values, Julia version and Manifest hash, solver settings and seed (`1234`), per-gate residual results (`residual`/`budget`/`labour`-or-`wage`, `overall`), actor, and artifact pointers. The `S = I + X − M` canary and `external_balance` are diagnostics, never gates.
- `--design` also updates the cell's `registry/scenarios.csv` row (status, pinned params, `evidence = "runs/<id>/manifest.toml; runs/<id>/log.txt"`, recorded HEAD).
- Historical runs predating manifests are registered with their caveats and are **not** retrofitted with invented manifests.

## Practicalities for this repo

- Environment: `julia --project=.` from the repo root. Package load check: `julia --project=. -e 'using BeyondHulten'`.
- Headless runs must not load GLMakie. Plotting is isolated in the lazy extension `ext/BeyondHultenGLMakieExt.jl` and loads only when both `GLMakie` and `BeyondHulten` are loaded, so `julia --project=. -e 'using BeyondHulten'` is headless-safe; if a script needs only the kernel, direct includes work too, as `rerun_results.jl` does:
  `include("src/core/accounting.jl"); include("src/core/technology.jl"); include("src/closures/labor/types.jl"); include("src/core/equilibrium.jl"); include("src/closures/labor/labor.jl"); include("src/closures/financing/financing.jl"); include("src/closures/registry.jl"); include("src/core/diagnostics.jl")`.
- Run `julia --project=. -e 'using Pkg; Pkg.test()'` before quoting any result.
- Paper tables and figures cite `run_id`s from `registry/scenarios.csv`, not numbers restated in documents (ADR-0003/ADR-0004).
