# ADR-0025 — Matrix visualisation: validated re-solve, headless data layer, figures are not runs

- **Status:** accepted
- **Date:** 2026-09-20
- **Supersedes:** —
- **Related:** ADR-0004 (runs are immutable), ADR-0005 (closure plugins and layout),
  ADR-0006 (entry point and preregistration), ADR-0019 (external account),
  ADR-0020 (sectoral wages at η = 0), ADR-0022 (sectoral labour markets);
  `src/plots.jl`, `ext/matrix_plots.jl`, `ext/BeyondHultenGLMakieExt.jl`,
  `experiments/plot_matrix.jl`, `experiments/README.md` ("Plotting the matrix"),
  `tests/test_matrix_plots.jl`

## Context

The executed matrix generations record per-run `manifest.toml`
metrics/diagnostics and a `solution.csv` with only `sector,price,quantity`.
Exploratory analysis of the shock effect across closures needs per-sector real
wages (sectoral at η = 0 / ADR-0020 option C and ADR-0022), household
consumption, import content and exports, and it must be comparable across the
fifteen labour × financing cells. Runs are immutable (ADR-0004) and a `src/`
change normally re-mints a generation; the existing GLMakie extension
(`ext/BeyondHultenGLMakieExt.jl`) carried only cbase2-era plotting functions
and the package had no run-artifact-based plotting API.

## Decision

Figures for an executed generation are produced by **re-solving** each selected
cell with the run harness (`solve_cell` + `evaluate_gates` from
`experiments/run.jl`, the same reference continuation and programme vectors as
`run_design`) and **hard-validating** each cell against its recorded artifacts:
manifest status `executed`, gates `pass`, manifest metrics within 1e-8, and
`solution.csv` prices/quantities within 1e-8 (`validate_cell`). A cell that
does not validate is excluded from the figures and the driver exits non-zero.

The design must be **preregistered with a matching SHA-256**
(`preregistration_status`), tying the figures to the pinned design.

The data layer is **headless** in `src/plots.jl` (`MatrixBaseline`,
`MatrixCellData`, `MatrixDataset`, `matrix_baseline`, `sectoral_trade_flows`,
`matrix_cell_data`, `matrix_cell_ids`, `matrix_sectoral_frame`,
`matrix_summary_frame`, `validate_cell`, and the `plot_matrix_*` /
`save_matrix_figures` stubs); the drawing lives in the GLMakie extension
(`ext/matrix_plots.jl`, included from `ext/BeyondHultenGLMakieExt.jl`), so
`using BeyondHulten` stays graphics-free.

`experiments/plot_matrix.jl` is the single plotting entry point; `--no-figures`
runs validation + data export without GLMakie.

The figures are **generated artifacts** (written under `plots/`, gitignored)
and are not committed; the tidy datasets are written under `output/`
(gitignored). No new generation is minted for visualisation, and every figure
is traceable to the `run_id`s through the printed validation report.

Comparability conventions: canonical cell order (`matrix_labour_order()` ×
`matrix_financing_order()`), percent units with labelled axes, one shared
colour scale per metric per figure, baseline = the design reference
(no-programme) solution; the six figures are `overview`, `wages`, `prices`,
`quantities`, `consumption`, `trade`. On `matrix_5x3_v10` the driver re-solved
all 15 matrix cells and validated 15/15 with
`max|Δp| = max|Δq| = max|Δmetric| = 0` against the recorded artifacts
(bit-identical reproduction; command
`julia --project=. experiments/plot_matrix.jl --design matrix_5x3_v10 --no-figures`);
the six figures render headlessly in the GLMakie 0.13.13 environment.

## Consequences

Re-solving requires the local (gitignored) IO table and the programme source,
like the probes; a future generation could store richer per-cell artifacts
instead (recorded as an option, not adopted); the figures are exploratory (not
paper tables) and are not a substitute for the manifest-backed flow tables.

## Enforcement

The driver refuses to plot a non-preregistered design or a non-validating cell
(non-zero exit); `tests/test_matrix_plots.jl` covers the data layer, the
trade-flow identity against `gdp_components`, `validate_cell` and the headless
stubs; `experiments/README.md` ("Plotting the matrix") documents the contract.
