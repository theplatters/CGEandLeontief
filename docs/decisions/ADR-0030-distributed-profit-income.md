# ADR-0030 — Distributed profit income: Stage B of the markup door

- **Status:** accepted
- **Date:** 2026-10-15
- **Supersedes:** —
- **Related:** ADR-0029 (the pricing door; this is its Decision 6 follow-up),
  ADR-0018 (GDP measurement), ADR-0019 (the external-account acceptance
  identity), ADR-0004 (runs are immutable manifests); `OpenCGE@ADR-0030`;
  `registry/closures.toml` `[pricing.P2]`; `ROADMAP.md` Phase 7, Route 2,
  steps 8–14

## Context

ADR-0029 introduced the markup door `p_i = μ_i c_i` and RETAINED its rent
`Π = Σ_i (1 − 1/μ_i) p_i y_i` (Stage A): booked as a leak beside saving, the
household budget untouched. Its Decision 6 named domestic profit income as
the developed target and required a separate decision, because paying the
rent to the household adds a third income foundation and forces a tax-base
and a welfare-metric choice. OpenCGE v0.2.0 implements the mechanics behind
`MarkupPricing(κ; s_pi)`; this record is the decision those mechanics carry.

## Decision

1. **Budget.** With a profit retention rate `s_π ∈ [0, 1]` the household
   receives the distributed rent:

       E = (1 − τ) w ΣL + F + (1 − s_π) Π

   and `s_π Π` is the retained leak. `s_π = 1` is Stage A EXACTLY (the
   `s_π = 1` branch is arithmetically separate in the kernel, so the nesting
   is bit-for-bit, not to round-off).
2. **Tax base (v1).** `τ` stays on wage income; the distributed rent is
   untaxed. The two-base variant (a one-line change to the `tau_rate`
   denominator) is a robustness row, not part of v1.
3. **Demand side.** `(1 − s_π)Π` enters household expenditure, so it reaches
   domestic demand through the import margin like any other household
   spending; the household then saves `s·E` of it. Both residual systems see
   it through the shared demand block.
4. **Accounting.** The external-account identity books saving as
   `s·E + s_π Π` (with `E` including the distributed rent); GDP is factor
   income `w·ΣL + Π` (ADR-0018's income side extended by the rent), so the
   companion identity `gdp_wedge = −canary_diff` still holds.
5. **Welfare.** The consumption index (`real_consumption`) is the
   consumption of the AUGMENTED household; a companion
   `real_income_index` (labour income plus distributed rent, CPI-deflated)
   and the profit share `Π / (w·ΣL + Π)` are recorded in every priced
   manifest.
6. **Design values.** `s_π ∈ {0.5, s}` with `s = 0.1199` the household saving
   rate (profit saved like wages), crossed with `κ ∈ {0.25, 0.5}` on F1–F3
   in the fixed-wage gauge: design `pricing_profit_v1`. The Stage A cells
   (`s_π = 1`) are `pricing_markup_v1` of the same generation.

## Consequences

- The model becomes "one kernel, three income foundations" only where a
  markup cell sets `s_π < 1`; every other cell keeps the one-factor income
  identity. Section 4's parsimony claim must be scoped accordingly.
- Distribution feeds back into demand, so Stage B is the variant in which
  the markup door is a distribution channel and not only a relative-price
  channel; the retention variant stays the cheap closure that keeps the
  household budget fixed.
- An untaxed rent understates the fiscal response under F2 relative to a
  two-base tax; the v1 tables must say so.

## Enforcement

- Kernel: `OpenCGE@tests/test_pricing_closure.jl` ("Stage B nests Stage A
  (s_π = 1) and the baseline (κ = 0)"; "Stage B budget, account and
  diagnostics").
- Paper: `tests/test_pricing_design.jl` — the real-table external account at
  `s_π = 0.1199`, and `evaluate_gates` on a real-table Stage B cell, whose
  budget gate includes the distributed rent.
- Every Stage B number cites a `pricing_profit-v1-*` run id (ADR-0004).
