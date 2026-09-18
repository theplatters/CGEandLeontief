# ADR-0017 — actual-matrix scale determinacy for the fixed-wage η = 1 system

- **Status:** accepted (user instruction, 2026-09-18)
- **Date:** 2026-09-18
- **Corrects:** the ADR-0014 closed-form round-gain trigger
  (`max(A_bill/λ + (1−m)(1−s)·fs) < 1`) for the fixed-wage η = 1
  admissibility guard. ADR-0014 is NOT superseded as a whole: its real-wage
  elastic-supply decision (DE-0004) stands.
- **Related:** ADR-0010 (kept BF endpoints η ∈ {0, 1}; the demand hook the
  guard shares), ADR-0012/ADR-0013 (A-bill calibration the clearing matrix is
  built on), ADR-0014 (the corrected trigger), ADR-0015 (the polish the
  admitted cells solve under), `src/core/equilibrium.jl` (`problem_fixed`,
  `_fixed_clearing_affine`, `_solve_fixed`), `tests/test_fixed_closure.jl`,
  `tests/test_kernel_regression.jl`, `tests/test_promoted_closures.jl`;
  found in the 2026-09-18 review of commit range `5ea54465..901ec07`

## Context

The ADR-0014 guard admitted the fixed-wage η = 1 system exactly when the
closed-form "round-gain" column sums
`A_bill_u/λ_u + (1−m_u)(1−s)·fs_u` peaked below 1. The review showed the
trigger is mathematically wrong twice, with executed counterexamples on the
`tiny_fixture` (fixed η = 1, elasticities (.5, .5, .9), manna
`autonomous_demand = [0.1, 0.0]`):

1. `maximum(colsums) >= 1` does not imply a unit root. With
   `A_bill = [0.5, 0.4]`, `M_int = [0.0, 0.1]` (s = m = 0) the actual
   clearing matrix is `G = [0.75 0.25; 0.25 0.65]`: column sums `[1.0, 0.9]`,
   eigenvalues ≈ {0.445, 0.955}, `det(I−G) = 0.025`, and an exact equilibrium
   at `p = [1, 1]`, `y = [0.7, 0.5]` with residual 5.6e-17 — yet the old guard
   rejected it.
2. The closed-form sums are not the actual matrix. Household income is spent
   across sectors, so the household feedback column sum is
   `(1−s)·fs_u·Σ_i((1−m_i)·b_i)` with `b_i` the normalized effective
   consumption weights — import margins enter at the SPENDING sector, not the
   receiving one. With `import_margin = [0.0, 0.2]` the actual column sums are
   `[0.95, 0.95]` while the old formula reads `[1.0, 0.9]` and falsely
   rejects; the exact equilibrium is `y = [0.6, 0.4]` (residual 5.6e-17).

User ruling: the strict column-sum bound is only a SUFFICIENT contraction
shortcut; otherwise the guard must assess the actual clearing matrix and
distinguish contraction, nonsingularity, and positive-solution existence.

## Decision

`problem_fixed` now evaluates its clearing block through
`_mobile_market_demand` (w = 1) — the single demand kernel shared with the
mobile system — so the guard's matrix cannot drift from the solver's
equations. The refactor is bit-identical (max residual difference 0.0 at
three (p, y) points on both the v3-style and `tiny_fixture` fixed models).

`_fixed_clearing_affine(model, p)` extracts the exact affine representation
`intermediary_demand + total_final_demand = G·y + c` at fixed prices (the
block is exactly affine in y) by unit differences at `y = λ` through the same
hook. `_solve_fixed` assesses determinacy on `(G, c)` at the calibration
baseline `p = ones(N)`:

1. CONTRACTION shortcut: `max(colsum(G)) < 1 − 1e-12` admits (G ≥ 0 at the
   reference, so this is sufficient for ρ(G) < 1);
2. NONSINGULARITY: otherwise `(I − G)` is tested by SVD; a singular
   `(I − G)` (σ_min/σ_max ≤ 1e-10) is a unit root — a continuum of solutions —
   and the cell is rejected with the legacy `"scale-indeterminate"` surface.
   The message keeps the `"autonomous or investment"` substring asserted by
   the contract tests, and states that additive demand (manna or a
   TaxFinanced / ExternalDebt bundle) is a constant that cannot remove a unit
   root;
3. POSITIVE-SOLUTION EXISTENCE: a nonsingular `(I − G)` admits only if the
   unique candidate `y* = (I − G)⁻¹c` is positive at the model resolution
   (min y* > 1e-12·max(1, |y*|_max)); otherwise the `"positive"` rejection
   fires.

Near-one η is snapped to the η = 1 endpoint for the assessment (the hook
admits only the kept endpoints), preserving the legacy surface that near-one
η on a closed fixture throws before any endpoint validation.

## Consequences

- Committed generations are invariant: on the full-71 A-bill calibration the
  ACTUAL max column sum is 0.831378781888721 (≈ 0.8314 < 1; the old formula
  gave 0.8559345674224728), and 0.8320 with the matrix design's F1 tilt — so
  every committed `matrix_5x3_v2`/`matrix_5x3_v3` fixed-wage cell takes the
  shortcut-admit path. No `-v4` generation is minted; v3 remains the cited
  generation (ADR-0004: v1/v2 stay as history).
- The test suite pins all three branches (`tests/test_fixed_closure.jl`):
  closed-fixture rejection (singular `(I − G)`), both counterexample
  admissions (residuals < 1e-6 at the exact quantities), and the
  non-positive-candidate rejection (`det(I − G) = 0.025` but `y* = 0`).
  `tests/test_kernel_regression.jl` and `tests/test_promoted_closures.jl`
  keep their assertions with actual-matrix wording (v3 fixture actual max
  column sum 0.8905 < 1).
- Scope: the assessment is at the baseline reference prices. Supply-shock
  equilibria with p ≠ 1 are admitted on the same structural grounds; the
  guard is a fixed-η = 1 corner check, not a general rank certification away
  from the reference.
- Open items unchanged: BETA `η_s` identification needs a supply-side
  scenario; `T_int` recycling and the raw table's 5.387 % residual stand.
