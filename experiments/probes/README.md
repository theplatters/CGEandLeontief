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

- `probe12_segmented_wages.jl` — Step 1 of `docs/WORKPLAN_SENSITIVE_PRICES.md`:
  the segmented-wage (sticky set + one free wage) prototype. Measures the
  endpoint nesting (S = empty bit-exact onto ALPHA/BETA; S = all onto GAMMA),
  the sticky-share ladder x financing, the tilted-pin case, the pure-pin level
  control, and the segmented pin-level behaviour. Result: the restricted
  segmented closure nests, closes the account and moves employment, but its
  prices and its free wage are demand-free (the free wage is pinned by zero
  profit + the CPI numeraire), so it delivers Door-1 heterogeneity, not
  demand-sensitivity.
- `probe13_sectoral_labour.jl` — the successor: the **general**
  sectoral-labour-market closure (workplan v4; unknowns [p; y; w(1:N); F],
  N sectoral supply conditions). Measures the `eta_s,i = 0` corner against the
  executed v6 `BF` row, the uniform-`eta_s` ladder, the two-group variant, the
  demand sensitivity across the financing cells, and the conditioning.
- `probe14_s1_s5_comparison.jl` — the five scenario slots `S1`..`S5` of the
  workplan measured side by side on full-71: `S1` the rigid `eta = 0` row,
  `S2` the uniform-`eta_s` closure, `S3` the two-group variant, `S4` the pinned
  wage vector (ADR-0021), `S5` the utilization/Verdoorn externality
  (`delta = +-0.5`, cost-hook reduced form). Discriminator: whether
  `max abs(p-1)` differs across F1/F2/F3 (demand-sensitive) or is identical
  (exogenous heterogeneity). Result: `S4` is the only heterogeneous-only slot;
  `S1`, `S2`, `S3` and `S5` are demand-sensitive, and `S5` reaches the same
  price magnitude as `S2` with a single parameter and no wage or employment
  movement.
