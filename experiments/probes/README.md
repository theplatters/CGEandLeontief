# Diagnostics (`experiments/probes/`)

These scripts are read-only diagnostics: no run dir, no `runs/index.csv` row,
no `src/` change — with one documented exception,
`probe6_reconstruct_run_journals.jl` below, which writes only missing
`runs/<id>/log.txt` journals (never manifests, index rows or solutions).
They rebuild the calibration exactly as `experiments/run.jl`
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
- `probe15_sectoral_supply.jl` — the sectoral closure under a **supply** shock
  (`A[1] = 1.2`), the route to identifying `eta_s`. Three row sets: supply
  shock only, supply shock + the programme under F2, and the nesting canary.
  Result: the `eta_s` rungs **separate** under the shock (employment
  1.000000 / 1.002187 / 1.008218 and the real-wage spread narrowing
  1.215 / 1.105 / 1.043 at `eta_s` = 0 / 0.5 / 2), so `eta_s` is identifiable
  only with a supply-side scenario — the executed demand-only ladder is a
  sensitivity band. The shock also **doubles the price response** and flips the
  sign of the welfare effect (supply-only `+1.20 %`, supply + programme
  `-0.40 %`, demand-only `-1.57 %` at `eta_s = 0.5`), so the two channels do not
  add linearly. The nesting canary holds **bit-exactly under the shock**
  (`max|dp| = max|dq| = max|dw| = |dF| = 0.00e+00` at `A[1] = 1.2`), confirming
  that Proposition 3 never uses `A = 1`. Caveat recorded: under a supply shock
  `max abs(p-1)` is no longer a demand-sensitivity measure (the technology shock
  moves prices with no programme at all, 0.226 at `eta_s = 0.5`); the demand
  channel must be read from the F1/F2/F3 differences.
- `probe16_rigidity_share.jl` — how much the **grouping rule** matters, by
  sweeping the rigid share `s` in {0, 0.25, 0.5, 0.75, 1} under two rankings and
  three rankings at `s = 0.5`. Result: the rule dominates the level. With the
  programme ranking the price response **saturates immediately** (0.2731 at
  `s = 0.25` against 0.2790 at `s = 1`, because the seven programme sectors are
  already inside the rigid group), while with the employment ranking it crawls
  (0.1057 to 0.1262) and only jumps when the whole economy is rigid. At
  `s = 0.5` the price response reads 0.2734 (programme), 0.1234 (largest
  employers) and 0.1068 (least exposed), against 0.1057 uniform — so *which*
  sectors are rigid matters far more than *how many*. Employment is
  non-monotone under the programme ranking (1.001255 / 0.996911 / 0.997300 /
  0.999370 / 1.000000) and its sign is rule-dependent: the band is
  `[-0.31 %, +0.16 %]` in employment and `[0.0969, 0.2790]` in the price
  response. Both rankings coincide at `s = 1` (0.279009 = the executed `BF`
  value), a consistency check on the sweep.

- `probe17_gamma_capacity.jl` — **GAMMA + the capacity/Verdoorn channel**, the
  one fixed-wage variant with no measurement. Applies the workplan's reduced
  form (`A_eff_i = (y_i/lambda_i)^(-delta)` in the unit-cost hook) inside the
  *fixed-wage* system, which is a simultaneous solve because the price block
  stops being separable; it is solved by the fixed-point iteration the reduced
  form implies (`A <- A .* (y/lambda)^(-delta)`, damped, re-solve), whose fixed
  point is the solution of the simultaneous system. Canary: `delta = 0`
  reproduces the executed `matrix_5x3-v9-GAMMA-*` cells exactly. Result:
  `delta = +0.5` (capacity pressure) gives a demand-sensitive price response
  with the wage still pinned — `max abs(p-1)` 0.098457 / 0.110761 / 0.128463
  across F1/F2/F3, deflator up to 1.021461 (F3), employment 1.004167 / 1.008489
  / 1.038814, consumption +0.46 % / -1.08 % / +2.20 %; `delta = -0.5` (the
  Kaldor-Verdoorn sign) *reverses* the price sign, the F3 deflator reading
  0.985336, with consumption -1.21 % / -2.46 % / +2.27 %. Iterations to
  convergence: 35-39 at `delta = +0.5` against 70-76 at `delta = -0.5` — the
  slower rate is the positive feedback of increasing returns (proximity to
  instability, not numerics). Contrast with probe14's `S5`, which applied the
  same externality to the *mobile* system (employment pinned at 1.0); here the
  wage is pinned and employment absorbs.
- `probe18_gamma_bf_eta.jl` — **GAMMA carrying the BF allocation elasticity**
  (`eta = 0`): the fixed wage *and* the frozen allocation, i.e. the
  double-rigidity corner. The BF closure's parameter is the allocation
  elasticity (ADR-0010), and no cell of the matrix combines it with the fixed
  wage. Result: employment is a datum (`L = 1.000000` in every financing
  column), prices and the deflator are pinned at 1, and the programme moves only
  the output composition (`max|dy/lambda - 1|` ~ 0.21) and the tax mix; under F3
  household consumption is *literally unchanged* (`+0.000 %`), because the
  programme is externally financed and real income is pinned on both sides. The
  intermediate values `0 < eta < 1` are not available — they were retired with
  the B&F reallocation wedge.
- `probe19_gamma_dual_labour.jl` — **GAMMA with a dual labour market**: the
  pressured sectors' employment is frozen (insiders at `Lbar_i`) while the rest
  adjusts at the pinned wage. The constraint enters through the constrained CES
  unit cost (labour fixed, the firm substitutes towards intermediates; the
  output ceiling `y_max = A alpha^(1/(eps-1)) Lbar` is finite because
  `eps = 0.5 < 1`). Canary: the empty mask reproduces the executed GAMMA-F2 cell
  exactly (`0.000000` / `1.000000` / `L 1.001260` / `cons -1.533 %`). Result: a
  **19x dualism signature in prices** (insider `max|p-1|` 0.0278 against
  outsider 0.0015 at F2 — the wages are identical in both segments), and the
  effect is **contractionary**: total employment `-1.04 % / -0.88 % / +0.60 %`
  and consumption `-1.34 % / -2.88 % / +0.69 %` across F1/F2/F3, because the
  bottleneck raises intermediate costs economy-wide. The constrained sectors sit
  at about half of their capacity ceiling (the baseline fraction `alpha_i`), so
  the constraint bites on the substitution margin rather than at a hard ceiling.
  Comparison: the sectoral family's rigid-programme variant (employment
  `-0.34 %`, consumption `-2.69 %` at F2) is *milder* on employment — a vertical
  supply curve still lets the wage absorb, a quantity constraint forces the cost
  up.

- `probe6_reconstruct_run_journals.jl` — the one probe that writes into
  `runs/`: rebuilds missing `runs/<id>/log.txt` journals from the committed
  manifests (`--verify` checks the journals that exist, `--write` writes the
  missing ones; journals only, never manifests/index rows/solutions). Used to
  reconstruct the v4-v9 journals after generations executed on another working
  copy (journals are gitignored working artefacts by design, ADR-0004).
