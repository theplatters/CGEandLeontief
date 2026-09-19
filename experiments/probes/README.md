# Read-only diagnostics (`experiments/probes/`)

These scripts are read-only diagnostics: no run dir, no `runs/index.csv` row,
no `src/` change. They rebuild the calibration exactly as `experiments/run.jl`
does and solve cells directly for measurement purposes only.

Run from the repo root with `julia --project=. experiments/probes/<script>.jl`.

The scripts need the local (gitignored) input table and the stored `runs/matrix_5x3-v5-*/solution.csv` warm starts; a checkout that carries only the tracked manifests cannot run them.

- `bf_identification.jl` — eta = 0 quantity identification (stored vs jittered
  start) and the BF factor-market gap.
- `supply_identification.jl` — supply-shock (+20% sector-1) BETA
  identification over the eta_s grid.
- `numeraire_invariance.jl` — demand-only programme scaling (k = 1/2/5/10);
  gauge invariance and F2/F3 neutrality.
- `analytic_headlines.jl` — Leontief/closed-form headline re-measurement on
  the v5 runs (capacity, neutrality, curvature).

These four are the measurements of the 2026-09-18 v5-kernel session (DE-0011,
BF identification, capacity, closed forms; the session entry in
`docs/log/2026-09.md`). The `probe1..probe11` scripts are the later lineage:
after ADR-0020/0021 and matrix_5x3_v6 the live registry citations are
`probe5*/5b/5c`, `probe7` and `probe11`; DE-0011 and the DELTA/F1/F2
closed-form gates still cite the four scripts above.
