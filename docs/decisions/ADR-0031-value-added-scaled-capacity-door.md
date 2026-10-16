# ADR-0031 — Value-added-scaled capacity door and the services-zero profile

- **Status:** accepted
- **Date:** 2026-10-16
- **Supersedes:** the `(y/λ)^(−δ)` exponent of ADR-0029 Decision 2 for the
  capacity door (P1); `pricing_capacity_v1` and the capacity cells of
  `pricing_mobile_v1` stay on record as the superseded generation
- **Related:** ADR-0029 (the pricing door), ADR-0004 (runs are immutable
  manifests), ADR-0028 (pin bumps are generation events); `OpenCGE@ADR-0031`;
  `registry/closures.toml` `[pricing.P1]`; `ROADMAP.md` Phase 7, Route 1

## Context

ADR-0029 put the capacity door on gross-output technology,
`A_eff_i = A_i (y_i/λ_i)^(−δ)`, and ran it at `δ = ±0.5` in
`pricing_capacity_v1`. A Verdoorn coefficient of that size is an estimate
for (labour or value-added) productivity, mostly in manufacturing, and close
to zero in services. Under the model's Hicks-neutral gross-output
technology, a gross-output TFP change `d ln A` raises value-added
productivity by `d ln A / v_i` (`v_i` the value-added share of gross
output). Applied to gross output, `δ = −0.5` therefore implied a value-added
elasticity of `0.5 / v_i`: above 1 for every sector with `v_i < 0.5`,
about 1.3 for glass (`v ≈ 0.38`, the largest v1 response) and 11.6 for
coke and refined petroleum (`v ≈ 0.04`).

OpenCGE v0.2.1 (`OpenCGE@ADR-0031`) scales the exponent:
`A_eff_i = A_i (y_i/λ_i)^(−δ_i v_i)`, `v_i = factor_share_i`. `δ` is now the
value-added (Verdoorn) elasticity (`β = −δ`), and it may be a per-sector
vector.

## Decision

1. **Pin.** The paper repo pins OpenCGE `v0.2.1`. `scripts/reproduce.jl
   --check` reproduces `matrix_5x3_v10` with zero difference (15/15) on the
   new pin, so the bump changes P1 only.
2. **Generation.** The capacity door is re-run as a new generation:
   `pricing_capacity_v2` (GAMMA, F1–F3, 15 cells) and `pricing_mobile_v2`
   (the capacity cells of the mobile row, ALPHA, F1–F3, 12 cells). The
   markup and profit designs (P2) are not affected by v0.2.1 and are not
   re-run.
3. **δ values.** The ADR-0029 bracket `δ ∈ {0, +0.5, −0.5}` is kept. It now
   reads as the value-added elasticity, a bracket rather than a calibration.
4. **Two δ profiles per cell.**
   - *Homogeneous* (`-dp05`, `-dm05`): every sector carries `δ`.
   - *Services-zero* (`-s0`): the goods sectors (CPA sections A–F, indices
     1–35 of the full-71 table: agriculture, mining, manufacturing, energy
     and water supply, construction) carry `δ`. The services sectors (CPA
     G–T, indices 36–71) carry `δ = 0`. This follows the Kaldor reading of
     Verdoorn's law as a law of goods production; services estimates are
     near zero.

   The services block is listed explicitly in the design
   (`[pricing] services`, `n_sectors = 71`) and is selected per cell by
   `delta_services = 0.0` (`cell_pricing` in `experiments/run.jl`). Each
   manifest records `delta`, `delta_services` and `services`.
   `tests/test_pricing_design.jl` asserts that the list is CPA G–T of the
   source table.
5. **Animations.** The P1 animations sweep the homogeneous δ only. Sectoral
   cells are not animated, and frames are validated only against
   homogeneous runs of the same generation.

## Consequences

- At the same `δ`, the gross-output exponent is weaker by the factor
  `v_i ∈ [0.04, 0.77]`, so P1 price and quantity responses shrink.
- The v1 capacity results are superseded for citation. They stay visible
  (ADR-0004).
- A sectoral `β_i` beyond the two-block profile (for example, estimated
  Verdoorn coefficients by sector) needs only a new `[pricing] services`-style
  list or an explicit vector; the kernel already accepts any per-sector `δ`.

## Results

See `registry/closures.toml` `[pricing.P1]` and the run manifests of
`pricing_capacity_v2` and `pricing_mobile_v2`.
