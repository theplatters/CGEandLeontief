# ADR-0006 — Experiment entry point and preregistration

- **Status:** accepted
- **Date:** 2026-09-17
- **Supersedes:** —
- **Related:** ADR-0001, ADR-0002, ADR-0003, ADR-0004, ADR-0005,
  `registry/preregistration.toml`, `registry/scenarios.csv`,
  `experiments/run.jl`, `experiments/designs/matrix_5x3.toml`,
  `src/core/calibration.jl`

## Context

Phase 2 promoted the closures into the root kernel but left two Phase 3
gaps: the v3 open-economy calibration still lived only in the frozen
`cbase2/src/calibration.jl`, and there was no single entry point for runs
— `experiments/run.jl`, run manifests, and `runs/index.csv` were all
"planned" (ADR-0004, the `experiments` skill). The 15 `matrix_5x3` scenario
rows sat at `TBD` parameters with no pinned design.

## Decision

**Single entry point.** `experiments/run.jl` is the only way to execute a
scenario: `--list` (cells with pinned parameters), `--preregister`
(writes `registry/preregistration.toml`, the only writer), `--design`
(batch execution with `--cell`/`--cells`, `--runs-dir`, `--budget-seconds`,
`--actor`). It is `include`-safe (no `main()` on include) and every I/O
function takes explicit `root`/`runs_dir` so tests run against temporary
roots.

**Design + preregistration gate.** A design file
(`experiments/designs/<design>.toml`, `schema_version = 1`) pins data,
programme, reference continuation, gates, and every cell. `--design`
refuses to start — before creating any run dir — unless the file's SHA-256
matches `registry/preregistration.toml`. Re-pinning after an edit
overwrites the record; history stays in git.

**Manifest / index lifecycle.** `runs/<run_id>/manifest.toml` is written
`running` at cell start and `executed` (all gates pass) or `failed`
otherwise at the end, with `log.txt` progress lines and `solution.csv`
(`sector,price,quantity`) on success. Per-cell exceptions are caught and
recorded (`[error]` with type/message; `[gates] overall = "fail"` without
per-gate values). `runs/index.csv` (header `run_id,date,design,closures,
status,gate_summary,headline_metrics,commit`) is rewritten on every status
transition, sorted by `run_id`. Gates: residual / budget / labour (mobile)
or wage (fixed+DELTA); the `S = I + X − M` canary and `external_balance`
are diagnostics, never gates. Only `runs/index.csv` and
`runs/*/manifest.toml` are tracked (`.gitignore`, verified with
`git check-ignore`).

**Run-id immutability.** An existing run dir is never overwritten; a rerun
registers a `-v2` variant `run_id`. `--design` also updates the cell's
`registry/scenarios.csv` row (status, pinned params, `evidence =
"runs/<id>/manifest.toml; runs/<id>/log.txt"`, recorded HEAD).

**Calibration port scope.** `src/core/calibration.jl` ports
`cbase2/src/calibration.jl` (`DATASET_VARIANTS`, `drop_sectors`,
`dataset_coverage`, `recalibrate_open`, with header notes, finiteness
gate, saving-rate assertion, clamping report, and the `cbroot` argument),
adapted to the root `Data` layout. `read_data` gains a backwards-compatible
`datadir` keyword. One deliberate deviation: `drop_sectors` *subsets* the
v3 vectors with `[keep]` (keeping `saving_rate`) instead of zeroing them,
so a dropped-but-not-recalibrated dataset keeps a valid
`household_baseline`. Verification: the 71-sector recalibration reproduces
the recorded cbase2 values exactly (`s = 0.3979`, `tau0 = 0.2141`,
`X = 0.4219`, `I = 0.1624`); the 70-sector values (`s = 0.4259`,
`tau0 = 0.2167`, `X = 0.4271`, `I = 0.1644`) differ by the documented
sector-71 exclusion. No other kernel behavior changes.

**Matrix pins.** `experiments/designs/matrix_5x3.toml`: θ = 0.5, ϵ = 0.5,
σ = 0.9; BF η = 0.5, ALPHA η = 1.0, BETA η = 0.5 with η_s = 0.5, GAMMA
η = 0.5, DELTA `delta_epsilon = 1e-4` (η = 1.0). DELTA cells pin
theta/epsilon/sigma = 1e-4 (the actually solved values). Programme:
`cbase2/data_raw/impulses.csv` year 2024 columns 3:73,
`total_eur_m = 40300.0`, F1 shift `1 .+ ψ`. `matrix_5x3-DELTA-F1` is
expected to fail the scale-indeterminacy guard (F1 has no additive anchor;
kept visible).

**SHA-256 provenance repo-wide.** Manifests record `design_sha256` plus a
`data_sha256` table (design file, IO table, calibration artifact,
programme source), the git commit/dirty flag, the Julia version, and the
`Manifest.toml` hash (or `"absent"`). `cbase2` remains the frozen source of
record (ADR-0001): it is read, never edited.

## Consequences

- The matrix batch can run only after this preregistration is committed;
  `scripts/check_repo.jl` (next) enforces manifest/index/scenario
  consistency.
- Contract tests for the tooling land next (independent agent); the human
  operator runs the real batch afterwards. Schemas are fixed in
  `experiments/README.md`.
- Statuses stay `planned` until cells execute; nothing is `validated`
  without Phase 4 gates plus a recorded manifest.

## Enforcement

- `experiments/run.jl --design` enforces the preregistration gate itself.
- `scripts/status.jl` validates scenario ids/statuses/evidence; the board
  must show 0 warnings.
