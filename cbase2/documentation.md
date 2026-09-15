---
title: "cbase2: Clean Notebook Pipeline for the 5 x 3 Evaluation Matrix"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-09-15"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [cbase2, pipeline, notebook, workplan, replication]
version: 1
last-updated: "2026-09-15"
---

This folder is the clean, self-contained pipeline for the Metroeconomica
revision (workplan Stages 1--2 of `docs/DOCS_ASSESSMENT.md`). It replaces
the ad-hoc top-level notebooks with one directed pipeline: raw data in,
5 x 3 evaluation matrix out, every stage leaving auditable artifacts. The
parent repo stays untouched; `cbase2` reads the parent only for kernel
drift control.

# Pipeline map

Run order is strict: each notebook reads only artifacts written by
earlier stages.

| Notebook | Workplan stage | Reads | Writes |
|---|---|---|---|
| `01_data_wrangling.ipynb` | Stage 1 (data, AC steps 0--3) | `data_raw/` | `data_processed/` |
| `02_accounting_consistency.ipynb` | Stage 1 (AC steps 4--7) | `data_raw/`, `01` output | `data_processed/AC_*.csv` |
| `03_financing_closures.ipynb` | Stage 1 item 1 (F1/F2/F3) | `data_processed/` | `data_processed/financing_*.csv`, kernel smoke tests |
| `04_labor_closures.ipynb` | Stage 1 items 2--3 (BETA, DELTA, CD guard) | `data_processed/`, `src/` | `results_intermediate/closure_tests.csv` |
| `05_preregistration.ipynb` | Stage 1 item 4 | `results_intermediate/` | `results_intermediate/preregistration.json` (+ md) |
| `06_residual_validation.ipynb` | Stage 1 item 6 | `src/`, `data_processed/` | `results_intermediate/validation_report.md` |
| `07_evaluation_matrix.ipynb` | Stage 2 item 1 | `05` preregistration, `scripts/` outputs | `results_final/matrix_5x3.csv`, inline figures |
| `08_sensitivity_robustness.ipynb` | Stage 2 items 2--4 | `07` output | `results_final/elasticity_table.csv`, `results_final/sobol_sectoral.csv`, inline figures |

Figures are produced inline in the notebooks, immediately next to the
evaluation commands that generate them, and saved to `plots/`. There is
no separate figures notebook.

# Directory layout

```{.text}
cbase2/
|-- documentation.md            this file
|-- 01..08*.ipynb               the pipeline (see map above)
|-- src/
|   |-- core/                   trimmed kernel, byte-identical copies of ../src files
|   |-- financing.jl            Stage 1.1: F1/F2/F3 closures (cbase2-only)
|   |-- closures.jl             Stage 1.2--3: BETA + DELTA corner (cbase2-only)
|   `-- validation.jl           Stage 1.6: residual gates (cbase2-only)
|-- scripts/
|   |-- diff_kernel.jl          drift check: core copies vs parent src/
|   |-- run_reference.jl        no-shock reference + calibration check (headless)
|   |-- run_matrix.jl           Stage 2 full 5 x 3 batch (headless)
|   `-- run_sobol.jl            sectoral Sobol batch (headless)
|-- data_raw/                   read-only inputs (never written by any stage)
|-- data_processed/             calibration artifacts (written by 01--03 only)
|-- results_intermediate/       tests, preregistration record, validation report
|-- results_final/              5 x 3 headline tables, appendix tables
`-- plots/                      figures saved by the notebooks
```

# Model kernel and drift control

`src/core/` holds byte-identical copies of the parent kernel files
(`interface.jl`, `solution.jl`, `ces.jl`, `mobile_labor.jl`,
`leontief.jl`, `util.jl`, `variance_decomposition.jl`), proven to run
without the heavy GLMakie/XLSX/Ipopt dependencies by the parent
`rerun_results.jl`. Copies carry no header comments so byte comparison
stays meaningful; provenance is this paragraph.

The hybrid rule: cbase2-specific code lives in `src/financing.jl`,
`src/closures.jl`, and `src/validation.jl`; the core is copied, not
referenced, so an assessor needs only this folder. Drift is flagged, not
silent:

```{.bash}
julia cbase2/scripts/diff_kernel.jl
```

`SAME` means the copy matches `../src/` byte-for-byte; `DIFF` is
legitimate once Stage 1 edits land in the kernel, but every `DIFF` must
be either backported to the parent or recorded in the notebook that
introduced it. Run this check before any batch run.

# Execution environment

- Julia 1.9 or newer. Notebooks are authored in the statistics container
  (`julia` at `/usr/local/bin/julia`); light cells (data wrangling,
  accounting, small solves) may be smoke-tested there.
- Heavy runs (full 5 x 3 matrix, Sobol batches) execute on the Mac via
  the headless `scripts/run_*.jl`, `julia --project=.` from the parent
  root for kernel dependency resolution, or with a local `Project.toml`
  if cbase2 is later cut loose. Batch outputs land in
  `results_intermediate/` / `results_final/`; notebooks then load and
  display them.
- No notebook performs a long solve in a cell. Solves belong in
  `scripts/`; notebooks orchestrate, display, and interpret.

# Conventions

- One-way data flow: `data_raw/` -> `data_processed/` ->
  `results_intermediate/` -> `results_final/`. No stage ever writes
  upstream.
- Pre-registration gate: `05_preregistration.ipynb` pins
  $\eta^{*} = 0.5$ (pre-registered before any Stage 2 inspection) and the
  expected-signature table of the assessment document into
  `results_intermediate/preregistration.json`, stamped with the git SHA
  and timestamp. `07` refuses to run without that record.
- The unfinanced autonomous demand shock of the original pipeline is
  retired here (Foundation II): every matrix cell is financed via F1
  (preference reallocation), F2 (tax-financed $g_i$ with
  $\sum_i p_i g_i = T$), or F3 (external debt, $\sum_i p_i g_i = F$).
- DELTA is computed as the GAMMA + Leontief corner, never as an
  independent equilibrium row.
- Seeds and pinned package versions are recorded in the header cell of
  each notebook that draws randomness (Sobol) or relies on solver
  tolerance (validation).
- Generated artifacts are gitignored except the small CSVs under
  `results_final/` needed for the paper.

# Status

- Scaffold created (2026-09-15): directory tree, `data_raw/` copies
  (SHA-256 recorded below), `src/core/` kernel copies verified `SAME` by
  `diff_kernel.jl`, drift script tested.
- Not yet present: notebooks 01--08, `financing.jl`, `closures.jl`,
  `validation.jl`, batch runners.
- Next step: build `01_data_wrangling.ipynb` and
  `02_accounting_consistency.ipynb` from
  `Notebooks/AccountingConsistency.ipynb` (its steps 0--7 map directly),
  then the financing core per workplan Stage 1 item 1.

# Raw-input provenance (SHA-256)

```{.text}
e9c3299731940f7ce683f0e508495c62ad255f0a5536cf1da63c02478058914d  data_raw/I-O_DE2019_formatiert.csv
5545c748d848ec121fa5aff21d0c46a55c8d00462518718d3d54374c7f0b8c5f  data_raw/impulses.csv
6e4a11289c81ba6c3aeae94617e2086732202e59f08e7191985df7fadc1f53e5  data_raw/sector_names.txt
```
