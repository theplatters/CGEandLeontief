# ADR-0032 — The pricing doors on every labour row

- **Status:** superseded in part by ADR-0037 (price reporting units)
- **Date:** 2026-10-18
- **Supersedes:** for citation, `pricing_markup_v1`, `pricing_profit_v1` and
  the markup cells of `pricing_mobile_v1` (OpenCGE v0.2.0), which stay on
  record (ADR-0004)
- **Related:** ADR-0029 (the pricing door), ADR-0030 (distributed profit
  income), ADR-0031 (value-added-scaled capacity door), ADR-0022 (sectoral
  labour markets), ADR-0028 (pin bumps are generation events);
  `registry/closures.toml` `[pricing.P1]`, `[pricing.P2]`; `ROADMAP.md`
  Phase 7

## Context

The pricing doors ran on two labour rows only: GAMMA (`pricing_capacity_v2`,
`pricing_markup_v1`, `pricing_profit_v1`) and ALPHA (`pricing_mobile_v1/_v2`).
ADR-0029's Consequences state that a door dissolves GAMMA ≡ DELTA and makes
BETA ≠ ALPHA reappear under demand-only shocks. Section 6.4 stated both
claims from the scope of Lemma 1 without measuring them, because the harness
rejected doors on scalar BETA and on DELTA.

Neither rejection was a kernel limitation. `delta_model` and `solve_beta`
are thin wrappers around `mobile_labor_model`, which accepts `pricing` for
every closure, and `pricing_admissibility` takes any
`Model{MobileLaborCES}`. The sectoral BETA form (ADR-0022) already accepted a
door.

## Decision

1. **DELTA takes the door in the harness.** `build_cell_model` builds the
   DELTA cell as `mobile_labor_model(…, delta_elasticities(ε)…; closure =
   :fixed, pricing)`, which is the `delta_model` call bit for bit under
   `NoPricing`. No kernel change and no pin bump.
2. **The BETA row of the pricing designs is the sectoral form** at uniform
   `η_s = 0.5` (`eta_s_rigid_group = "none"`, the BETA-nesting value of
   `matrix_5x3_v10`). Scalar BETA (`solve_beta`) keeps rejecting doors.
3. **Every door runs on all five labour rows** (BF, ALPHA, BETA, GAMMA,
   DELTA) × F1–F3:
   - `pricing_capacity_rows_v1` (48 cells): BF, BETA and DELTA on the full
     capacity grid (δ ∈ {0, ±0.5}, homogeneous and services-zero), plus the
     δ = 0 ALPHA canaries. Same pin as `pricing_capacity_v2` /
     `pricing_mobile_v2`, so the capacity table is one generation.
   - `pricing_markup_v2` (45) and `pricing_profit_v2` (60): the markup door
     re-minted on OpenCGE v0.2.1 across all five rows, so the markup table is
     one generation as well. v0.2.1 leaves P2 untouched; the GAMMA and ALPHA
     cells must reproduce their v1 manifests.
4. **Collisions are kept.** Under demand-only shocks with the door closed,
   DELTA reproduces GAMMA (Lemma 1). The δ = 0 / κ = 0 cells record that
   collision; the priced cells measure its dissolution.
5. **Manuscript.** Table 6 carries the capacity door by labour row and the
   new Table 7 the markup door by labour row (`make_pricing_door.py`).

## Consequences

- The equivalence-dissolution claims of Section 6.4 are now measured, not
  only stated from Lemma 1's scope.
- The markup results are cited from the v2 generation; the v1 manifests stay
  visible.

## Results

See `registry/closures.toml` `[pricing.P1]`, `[pricing.P2]` and the run
manifests of `pricing_capacity_rows_v1`, `pricing_markup_v2` and
`pricing_profit_v2`.
