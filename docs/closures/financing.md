# Financing closures — F1, F2, F3

Implementation: `src/closures/financing/financing.jl` (supertype
`AbstractFinancing` in `src/core/accounting.jl`); hooks consumed by
`src/core/equilibrium.jl` (`problem`, `problem_fixed`, `solve`). Status of
each closure lives only in `registry/closures.toml` (`[financing.F1/F2/F3]`)
and renders on `docs/status.md` (ADR-0003); this page carries no status.

## Formulations (ADR-0002)

- F1: `β̃ᵢ = βᵢdᵢ / Σⱼβⱼdⱼ`; `Σ pᵢcᵢʰ = Eʰ`; budget-neutral composition shift.
- F2: `Σ pᵢgᵢ = T(p)`; lump-sum / balanced-budget tax.
- F3: `Σ pᵢgᵢ = F`; external balance.

F1 is demand composition, never an autonomous investment multiplier; F2 needs
a balanced-budget counterparty; F3 is an accounting open-economy closure, not
a monetary mechanism.

## Design notes (preserved from the frozen `cbase2/src/financing.jl` header)

Programme: real bundle `g_i = m · G0 · ψ_i` (numeraire units), with `ψ` the
sectoral incidence (year-1 impulse shares) and G0 = 40,300 EUR m scaled into
model units (1 model income unit = GDP at basic prices).

- F1 `PreferenceReallocation(d)`: preference weights become
  `demand_shock .* d`; the household CES normalizer keeps `Σ p_i c_i = E`
  exactly. NO additive demand. Design tilt (ports
  `cbase2/scripts/verify_v3.jl`): `d_i = 1 + G0·ψ_i/c0_i` with `c0` the
  household baseline and `ψ` restricted to positive-baseline sectors and
  renormalised; the retired `1 .+ ψ` stand-in was a non-reference tilt
  (ψ added directly instead of `G0 .* ψ1 ./ c0`).
- F2 `TaxFinanced(g)`: additive real demand with endogenous lump-sum tax
  `T(p) = Σ p_i g_i`; household expenditure `E = w·ΣL − T`.
- F3 `ExternalDebt(g)`: identical additive demand, household untaxed; the
  external balance `F = Σ p_i g_i` is recorded post-solve.

Hooks: `preference_weights`, `tau_rate` (homogeneous budget: real purchases
financed at CURRENT prices), `household_expenditure` (`E = (1−τ)wΣL`),
`additive_demand` (g under F2/F3, zeros otherwise),
`public_budget`, `external_balance`, `has_additive_anchor` (only F2/F3 anchor
the fixed-wage η ≈ 1 scale).

## What was promoted (Phase 2, ADR-0005)

- All four types and all seven hooks verbatim, modulo the root-module paths.
- Compatibility (documented, intentional deviation from cbase2): the legacy
  unfinanced autonomous/investment manna (A/G) is RETAINED alongside the
  financed demand so the characterization goldens and `rerun_results.jl`
  stay reproducible; with `NoFinancing` and zero v3 fields the demand is
  numerically identical to the pre-Phase-2 code.
- `cbase2/src/calibration.jl` was NOT promoted (Phase 3).

## Remaining gates and caveats

Source: `registry/closures.toml` (`[financing.F1/F2/F3]` `open_gates`) and
`cbase2/review.md` §§1–2.

- Contract tests for the promoted code land in the next step (kept in
  `open_gates`; status stays `implemented`); matrix cells not yet run.
- `public_budget` prices at `sum(gov_demand)` while taxes are
  `dot(p, gov_demand)` (review §2.7): government budget reports fail as soon
  as p ≠ 1.
- `external_balance` returns only the programme's import content (review
  §2.8); the recorded F figure is ambiguous.
- F1 must never be described as an autonomous investment multiplier.
