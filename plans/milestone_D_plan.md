# Plan for Milestone D — Honest Variance Decomposition

## Current state
- `variance_decomposition.jl` computes one-way partial R² on a full factorial grid.
- `summary_table` renormalizes over main effects (double-counts), masking the true interaction residual.
- ~30% of grid points fail at ε=0.1 (complementarity), producing NaN entries that are silently dropped.
- No CSV output.
- The original "η=100%" result was an artifact of the broken GDP collapse + renormalization.

## What Milestone D requires
1. **Sobol first-order S_f** — same as current SS_f/SS_total (correctly computed).
2. **Sobol total-order S_Tf** — 1 − SS_{-f}/SS_total, where SS_{-f} is variance explained by all factors except f.
3. **Report absolute shares** — no renormalization over main effects.
4. **CSV output** — the result survives container resets.
5. **Fix summary_table** — show S_f, S_Tf, interaction strength S_Tf−S_f, and residual = 1−sum(S_f).
6. **Handle missing grid points robustly** — clearly report how many points failed and which regions.

## What's blocking
- Solver stalls at ε=0.1 (complementarity) for many θ/σ combos. Need more robust solver settings or continuation.
- The 5×3×3×3 = 135-point grid takes ~5 min. Larger grids or Saltelli sampling would take much longer.

## Implementation plan
- **Step 1**: Improve solver robustness — increase maxiters, use a closer warm start (run ϵ=0.5 first, use its solution as init for ϵ=0.1).
- **Step 2**: Extend `SobolResult` with total-order indices (`VarianceDecompositionResult` remains a deprecated alias).
- **Step 3**: Rewrite `_compute_partial_r2` → `_compute_sobol_indices` that returns both S_f and S_Tf.
- **Step 4**: Rewrite `summary_table` → absolute shares, no renormalization, S_f + S_Tf + interaction strength, residual, CSV export.
- **Step 5**: Write `diag_vd_honest.jl` — runs the full decomposition, validates, saves CSV.
- **Step 6**: Update WORKPLAN3.md Milestone D to VERIFIED.
