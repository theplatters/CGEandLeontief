# ADR-0010 — Port the `cbase2/review.md` pipeline fixes

- **Status:** accepted
- **Date:** 2026-09-17
- **Supersedes:** —
- **Related:** ADR-0001, ADR-0002, ADR-0005, ADR-0006, ADR-0007,
  `cbase2/review.md`, `registry/closures.toml`, `registry/freeze.toml`,
  `experiments/designs/matrix_5x3.toml`

## Context

`cbase2/review.md` (pinned to `17c5730..3dd0d60`) reported concrete pipeline
bugs. The `revisefinal` branch (`644ba37`, `3fa6d45`, `b46912d`) fixed the
retained-sector accounting and the mobile equation count, and those fixes had
to be ported into this branch's canonical root kernel (`cbase2/` is frozen at
`dfd60a2`, ADR-0001).

Two things surfaced during the port that changed the plan:

1. **The `revisefinal` all-N mobile formulation is not well posed for open
   economies.** It enforces all N clearing equations and drops the CPI
   numeraire. The residual is then homogeneous of degree 1 in `(p, w)`, so the
   equation count exceeds the effective unknown count by one. Closed fixtures
   solve because the N-th market is functionally dependent there; open
   economies with additive demand do not. Measured on the v3 contract fixture
   at η = 1: multi-start Gauss–Newton converges to a least-squares stationary
   point with `‖J'r‖ = 3e-12` and residual floor `9.0e-3` (F3) / `4.7e-4` (F2)
   — no root exists. On the real 70-sector calibration at the
   `matrix_5x3` reference rung (`exo_scale = esc0`, θ = 2.0) the all-N solve
   stalls at `4.41e-4`, the same order as the pre-existing N−1+CPI stall
   (`4.27e-4`). The review's own recommendation was to *either* impose the
   external balance *or* define it residually and assert the canary; the
   `revisefinal` commits did neither.
2. **Only the BF endpoints η ∈ {0, 1} are kept** (project decision). The
   geometric interpolation `0 < η < 1` and the ad hoc B&F (2019)
   "allocative-efficiency wedge" that carried it are retired. Without the
   wedge, η = 0 and η = 1 differ only in the reported sectoral allocation and
   the omitted-market residual; the wedge was never derived as a CES
   allocative-loss coefficient (review §1, §2.4).

## Decision

1. **Retained-sector pipeline (review findings 2.2/2.3).** Replace the
   `drop_sectors` matrix slicer with `retained_io_table` +
   `retained_dataset`: slice the RAW IO table, re-run `generate_data`
   (`number_sectors` keyword, final-demand columns located by name) and assert
   `Σ_j Ω_raw[u,j] = 1` and `Σ labor_share = 1`. `recalibrate_open(data;
   exo_scale)` now reads the open-economy blocks from `data.io` directly (no
   `cbroot`/`drops` arguments, no `AC_domestic_final_demand.csv`).
2. **η ∈ {0, 1} only.** `_checked_eta` rejects intermediate values with a
   `DomainError`; the `eta_sweep`/`variance_decomposition` grids and the
   matrix design use `{0, 1}` only. The allocation wedge
   (`_allocation_efficiency_wedge` and its use in `problem`) is deleted.
3. **Mobile formulation: N−1 clearing + CPI = 1 (the pre-port form), with the
   omitted market made explicit.** The N-th clearing is not Walras-redundant
   once imports leak; it is the residual external account.
   `market_clearing_residuals` exposes the full N-vector and
   `external_balance_canary` computes `S − (I+X−M)`; acceptance tests and
   `experiments/run.jl` assert the two agree at mobile η = 1 solutions. The
   all-N `revisefinal` formulation is rejected as unsound for open economies.
4. **Finding 2.5 (baseline clamp) stays open.** `recalibrate_open` continues to
   report the clamp mass; no reconciliation is attempted here.

## Consequences

- The calibration numbers change for the 70s variant: `s = 0.410965` (was
  0.4259), `E_h0 = 0.783259`, clamp mass `0.086197`; the full table is
  unchanged (`s = 0.397878`). The reduced variant recalibrates to
  `s = 0.386970`, clamp mass `0.091363`.
- Mobile kernel goldens move: the wedge removal changes the η = 0 solution,
  and η = 0/η = 1 no longer coincide in aggregate.
- The `matrix_5x3` design is re-pinned: BF cells use η = 0 (immobile endpoint),
  BETA/GAMMA/DELTA/reference use η = 1; `registry/preregistration.toml` is
  re-registered.
- The real-data mobile stall (`≈ 4.3e-4` at the reference rung) is *not*
  caused by the equation count; it remains the review §4
  scaling/conditioning/fold open item and is recorded in
  `registry/freeze.toml` / `registry/closures.toml`.
- Findings 2.6–2.9, 2.10, 3.x and 4 remain as recorded in
  `registry/closures.toml` `open_gates`; this ADR does not claim them fixed.

## Enforcement

- `tests/test_calibration.jl` asserts the retained-pipeline normalization and
  the real-data numbers (guarded on the gitignored IO table).
- `tests/test_kernel_regression.jl` pins the endpoint goldens and rejects
  intermediate η; `tests/test_promoted_closures.jl` asserts the
  external-account canary identity.
- `tests/test_eta_sweep.jl` / `tests/test_variance_decomposition.jl` reject
  retired η values; `experiments/run.jl` refuses non-`{0,1}` η grids and
  asserts the canary for mobile η = 1 cells and the reference.
- `scripts/check_repo.jl` enforces preregistration integrity for the re-pinned
  design.

## Amendment 2026-09-17 — `revisefinal` `0f33ad6` converges on the N-1+CPI form

`revisefinal` commit `0f33ad6` restored the same N-1+CPI mobile form and
bypassed the allocation wedge, so both branches have converged on the
formulation this ADR records. The reorg root kernel stays stricter: the wedge
is deleted (not bypassed), η ∈ {0,1} is validated by `_checked_eta`, and the
omitted market is exposed as `market_clearing_residuals` /
`external_balance_canary` with test and experiment assertions. The `0f33ad6`
F1 preference tilt is ported into the experiment pipeline under the new
`f1_shift = "tilt_g0_over_c0"` value. The `0f33ad6` `autodiff` keyword fix
applies only to the frozen `cbase2/src/solvers.jl` and is intentionally not
ported (ADR-0001/ADR-0005). The `0f33ad6` θ-consistent reference fix is not
needed by the current `matrix_5x3` design because cells and reference share
θ = 0.5; if the intended headline spec is θ = 1.0, that remains a design
question. The original wording above is kept per the append-only rule; this
amendment is authoritative where the two disagree.
