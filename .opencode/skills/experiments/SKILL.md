---
name: experiments
description: Use when defining, running, or recording an experiment or scenario in the BeyondHulten repo — run ids, registry/scenarios.csv rows and evidence paths, failed and provisional run retention, the ADR-0004 run-manifest contract, the planned experiments/run.jl entry point and runs/index.csv, headless Julia runs that avoid GLMakie, and citing run ids in the paper.
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

## Run manifests (ADR-0004)

Every quantitative claim is traceable: claim → `run_id` → manifest → commit → closure/design → data vintage.

- **From Phase 3 (planned; does not exist yet):** the single entry point `experiments/run.jl` writes `runs/<run_id>/manifest.toml` and appends to `runs/index.csv`. `runs/` holds raw results, logs, and figures and is not committed, except the index and manifests. Design files will live in `experiments/designs/<design>.toml` (planned).
- Manifest contents: git commit and dirty flag, scenario/design hash, data-vintage SHA-256 values, Julia and package versions (Manifest hash), solver settings and seed, per-gate residual results, actor (human/agent), and artifact pointers. `runs/index.csv` rows: `run_id`, date, design, closures, status, gate summary, headline metrics, commit.
- **Until then:** register runs manually in `registry/scenarios.csv` with existing evidence artifacts (the `cbase2-v3` rows cite `cbase2/scripts/verify_v3.jl`, `cbase2/process_comments.md`, `cbase2/review.md`), set `commit` to the recorded commit, and state the exact data vintage and configuration in `notes`.
- Historical runs predating manifests are registered with their caveats and are **not** retrofitted with invented manifests.

## Practicalities for this repo

- Environment: `julia --project=.` from the repo root. Package load check: `julia --project=. -e 'using BeyondHulten'`.
- Headless runs must not load GLMakie. Plotting is isolated in the lazy extension `ext/BeyondHultenGLMakieExt.jl` and loads only when both `GLMakie` and `BeyondHulten` are loaded, so `julia --project=. -e 'using BeyondHulten'` is headless-safe; if a script needs only the kernel, direct includes work too, as `rerun_results.jl` does:
  `include("src/interface.jl"); include("src/solution.jl"); include("src/ces.jl"); include("src/mobile_labor.jl"); include("src/util.jl"); include("src/variance_decomposition.jl")`.
- Run `julia --project=. -e 'using Pkg; Pkg.test()'` before quoting any result.
- Paper tables and figures cite `run_id`s from `registry/scenarios.csv`, not numbers restated in documents (ADR-0003/ADR-0004).
