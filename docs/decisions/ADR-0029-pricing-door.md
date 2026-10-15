# ADR-0029 — The pricing door: an output-dependent unit cost (Foundation III)

- **Status:** accepted
- **Date:** 2026-10-14 (proposed); accepted 2026-10-15 with OpenCGE v0.2.0
- **Supersedes:** —
- **Related:** ADR-0002 (closure taxonomy and ids), ADR-0004 (runs are
  immutable manifests), ADR-0012 / ADR-0013 (the domestic-bill calibration and
  the `T_int` leak), ADR-0016 (CES valuation of the intermediate-bill leaks),
  ADR-0017 (scale-determinacy guard), ADR-0018 (GDP measurement), ADR-0019
  (the external-account acceptance identity), ADR-0022 (sectoral labour
  markets — the labour door), ADR-0028 (paper/kernel split);
  `OpenCGE@ADR-0016`, `OpenCGE@ADR-0019`; DE-0001, DE-0004 (dead ends not
  revisited here); `IDEA-0002` (five doors to demand-sensitive prices);
  `ROADMAP.md` Phase 7; `experiments/probes/probe17_gamma_capacity.jl`,
  `experiments/probes/probe20_pricing_door.jl`

## Context

Foundations I (labour) and II (financing) are demand-side. Under a
demand-only shock (`A ≡ 1`) the zero-profit block is homogeneous of degree one
in `(p, w)` and contains no demand term, so it pins prices at one in both
gauges and the real wage at its anchor; the quantity block then reduces to an
affine multiplier system in which none of `(θ, ϵ, σ, η_s)` appears. That is
Lemma 1 and Lemma 2 of `paper/equivalence.tex`, and it is why ALPHA ≡ BETA and
GAMMA ≡ DELTA hold as propositions rather than as numerical coincidences.

Closure variety inside the demand-shock story therefore requires demand to
enter the PRICE block — and the kernel's cost line says where it can:

    cost = _ces_unit_cost(A_eff, factor_share, w, intermediate_price, ϵ)

Exactly three arguments are demand-injectable: `A_eff` (capacity / Verdoorn),
`w` (sectoral labour markets — the labour door, ADR-0022, already executed) and
`ip` (endogenous import prices). A demand-sensitive markup enters the same line
as a multiplier but creates a rent, which is a different KIND of change because
it needs an income decision. The substitution elasticities `(θ, ϵ, σ)` are NOT
a door: Lemma 2 makes the aggregate invariant to them at every value, so
elasticity variation can never carry the demand-shock punchline.

`ROADMAP.md` Phase 7 records two routes: Route 1 (capacity / Kaldor–Verdoorn)
and Route 2 (distribution, i.e. the markup). Both are specified in
`docs/ETAs.md`; `IDEA-0002` is `scoped` and ranks the full menu.
`probe17_gamma_capacity.jl` measured Route 1 read-only.
`probe20_pricing_door.jl` (commit `79d5cc8`) measured both routes on the real
retained-sector calibration — 75 cells across F1/F2/F3 — and it falsified one
canary and exposed one FOC site, both of which this decision fixes in the
kernel rather than leaving to the probe.

## Decision

**1. The pricing axis becomes Foundation III.** Three closures, in the kernel
(`../opencge`), file `src/closures/pricing/pricing.jl`:

- `NoPricing()` — the baseline: demand-invariant unit cost, `p = c`.
- `CapacityPricing(δ)` — `A_eff_i = A_i (y_i/λ_i)^(-δ)`: `δ > 0` capacity
  pressure (cost rises with own output), `δ < 0` the Kaldor–Verdoorn sign.
  Rent-free: the gain passes into cost under zero profit, so NO income closure
  is required.
- `MarkupPricing(κ)` — `p_i = μ_i c_i` with `μ_i = 1 + κ (y_i/λ_i - 1)`.
  Creates a rent.

All three nest exactly: `δ = 0` and `κ = 0` reproduce the baseline bit for bit,
so every existing number is untouched.

Closure ids (ADR-0002 rule: adding an id is an ADR-level change): `P1` =
`CapacityPricing`, `P2` = `MarkupPricing` (Stage A and, via `s_pi`, Stage B of
ADR-0030). The baseline `NoPricing` carries no id, as `NoFinancing` carries
none on the financing axis. Scenario ids keep the ADR-0002 shape
`<design>-<labor>-<financing>[-<variant>]`; the door parameter is the variant
suffix (`-dm05` for `δ = -0.5`, `-k025` for `κ = 0.25`, `-spi05` for
`s_pi = 0.5`).

**2. `Model` gains a `pricing` field defaulting to `NoPricing()`.** The 3-arg
and 4-arg constructors keep working, so no existing cell, test or client
changes.

**3. One hook set, one definition per expression:** `pricing_technology`,
`pricing_markup`, `pricing_foc_price`, `pricing_profit`, `pricing_bill_factor`.

**4. `p` has TWO roles, and telling them apart is normative, not stylistic.**
The purchaser price (household and intermediate demand take it as given) stays
`p`; the multiplier of the first-order condition `w = c ∂Φ/∂L` is the MARGINAL
COST `c = p/μ`. Using `p` there leaves the price response right and every
labour demand inflated by `μ` — the dangerous failure mode. The FOC pattern
`p^ϵ A^(ϵ-1)` occurs at FOUR sites, and all four must follow the closure:

- `_cost_minimizing_labor` — the wage-FOC inversion;
- the INNER `p^ϵ` of `intermediary_demand` (the OUTER `p^(-θ)` is the demanded
  input's purchaser price and stays);
- `external_balance_canary`'s bill factor `k = p^ϵ A^(ϵ-1) ip^(1-ϵ)`
  (ADR-0016);
- `gdp_components`' `k` — previously a hand-copied duplicate of the third,
  now calling the shared `pricing_bill_factor`.

A fifth occurrence, `economy_wide_wage` (the wage FOC inverted to build solver
starting points), never enters a residual; it reads `pFOC` and `A_eff` through
the same hooks all the same, so the "no fifth site" grep below holds without
an exemption.

Sites 3 and 4 are reached by the accounting gate rather than by the solve, and
a synthetic fixture CANNOT see them: the leaks they value (`M_int`, `T_int`)
are zero in anything the compact `Data` constructor builds, so every identity
stays green at machine precision while the real calibration is wrong. The
identity must therefore be asserted on a fixture whose purchaser-price bill is
split by hand (`A_bill + M_int + T_int = (1-fs) λ`).

**5. Capacity needs no FOC change; the markup does.** The capacity externality
moves COST, so `p = c` survives and the FOC multiplier stays `p`. That
asymmetry is why Route 1 is a one-line hook and Route 2 is cross-cutting.

**6. Rent routing: retention (Stage A) is the default.** The rent is booked as
a leak beside the saving-rate leak `s E`; the household budget is unchanged and
all effects run through relative prices. Foreign owners (fits the external
account) and a government wedge (the state is already the F2/F3 instrument) are
the cheap alternatives. DOMESTIC profit income (Stage B,
`E = (1 - τ) w ΣL + (1 - s_π) Π`) adds a THIRD income foundation, forces a
welfare-metric and a tax-base decision, and requires its own ADR.

**7. Retention is feasible only where the economy has injections.** The
aggregate add-up is

    (μ - 1) c y  +  (s + m - s m) w L  =  p (1 - m)(gG + exo) + pX

so with no injections and no leakages the fixed-wage system has NO positive
solution — the classical demand-shortfall condition, not a solver failure. The
mobile gauge is always well posed because the free external transfer absorbs
the rent, which means the retained rent is then effectively externally
financed. The retention variant is therefore confined to designs carrying
injections; the paper's application does.

**8. The promoted implementation is the KERNEL hook, not the probe's outer
fixed point.** The probe iterates the map `A <- A (y/λ)^(-δ)` and can diverge:
at `δ = -0.5` under Cobb–Douglas it diverges at damping 0.6 and converges only
at 0.2. With `A_eff` endogenous inside the residual the kernel solves one
simultaneous system. The probe stays as the cheap instrument for measuring
before an ADR; it is not the promoted form.

**9. Canaries and acceptance readings.**

- `δ = 0` and `κ = 0` reproduce the baseline bit for bit.
- The uniform-markup closed form is the CASCADE,
  `p = μ̄^((I-B)^(-1) 1)` with `B = diag(1-fs) Ω_raw` (cascade ≈ 1.8 on the
  German table), NOT `p = μ̄`. A canary written as `p = μ̄` rejects a CORRECT
  implementation.
- The external-account gate has two documented readings: with zero legacy manna
  (the matrix designs) the identity
  `S + T_int + M - (I+X) = F + B_gov` holds AS WRITTEN with the rent term
  booked; on a legacy-manna design the gap legitimately reads `p (A + G)`.
  Assert against the documented value, never against zero.
- The rent's placement is confirmed both ways: adding it restores the
  documented reading, and removing it leaves the identity short by exactly
  `Π_ret`.

**10. Presentation decisions for sections 6.3 / 6.4.**

- Routes 1 and 2 are ONE family with two signs and a residual:
  `p_i = (1 + κ (y_i/λ_i - 1)) A_i (y_i/λ_i)^(-δ) c`. They are presented as a
  price-response BRACKET — a result — not as two unrelated robustness arms.
- The results table carries a SIGNED price statistic; a lone `max|p - 1|`
  hides the financing-conditional sign.
- The bracket is financing-conditional: at `δ = -0.5` the GDP deflator falls
  only under F3 (0.9853) and RISES under F1 (1.0035) and F2 (1.0013), so
  "Verdoorn lowers prices" is an F3 statement, not a property of the door.
- The scope clause is load-bearing: ALPHA ≡ BETA, and Lemma 1 itself, survive
  only on the slice "`A ≡ 1` AND price-invariant technology".

**11. Settled at acceptance (2026-10-15).** The items the proposal deferred:

- **The `μ ≥ 1` domain.** `μ = 1 + κ (y/λ - 1)` stays SYMMETRIC about the
  baseline: a contracting sector receives a discount (measured
  `μ = [1.0497, 0.9480]` at `κ = 0.5`, probe20) and its sectoral rent is
  negative. This keeps the exact nesting at `κ = 0` and a smooth Jacobian at
  the baseline, which the one-sided form `1 + κ·max(y/λ - 1, 0)` would kink.
  The discount reading is stated, not hidden: the kernel guards `μ > 0`
  (`PricingFeasibilityError`) and every priced manifest records `markup_min`
  / `markup_max`.
- **The `δ < 0` uniqueness guard.** A gate in the paper harness
  (`experiments/run.jl`, `evaluate_gates`), on every priced cell: at least two
  of the `pricing_admissibility` starts converge, all within
  `[pricing] multistart_tol = 1e-6` of the solution, and the residual
  Jacobian has full rank. The tolerance is the residual gate's: the measured
  solver-noise spread is 0.5-8e-9 at condition numbers ~1e7, a distinct
  equilibrium would sit ~1e-2 away. Each priced cell is reached by a
  deterministic continuation in its door parameter from the NoPricing twin
  (`continuation_steps = 4`).
- **Stage B** is ADR-0030 (distributed profit income).
- **The mobile-gauge ladder** for both routes is the robustness design
  `pricing_mobile_v1` (ALPHA × F1-F3 × {δ = ±0.5, κ ∈ {0.25, 0.5}}).

## Consequences

- Any demand-sensitive pricing door breaks Lemma 1's `p ≡ 1` and therefore
  DISSOLVES GAMMA ≡ DELTA and makes BETA ≠ ALPHA reappear under demand-only
  shocks. The paper's headline equivalences now depend on an explicit scope
  clause rather than holding unconditionally.
- Promotion is a generation-level event: an OpenCGE rev bump invalidates the
  provenance of existing runs (ADR-0028), so the door must be minted as a new
  generation rather than mixed into `matrix_5x3_v10`.
- Route 1 is cheap in code and budget; Route 2 is cheap in code and expensive
  in budget, because it forces the routing decision of Decision 6 and the
  feasibility condition of Decision 7.
- The kernel's one-factor, zero-profit, no-profit-income identity survives only
  under retention or an external/state routing. A domestic profit block would
  make the model "one kernel, three income foundations" and dilute the
  parsimony claim section 4 rests on.

## Enforcement

- Kernel tests, `OpenCGE@tests/test_pricing_closure.jl` (OpenCGE v0.2.0): the
  three nesting canaries bit-identical; a demand-only shock pins prices at one
  WITHOUT a pricing door and moves them with one; `w = c ∂Φ/∂L` holds at a
  markup solution while `w = p ∂Φ/∂L` does not; the external-account identity
  holds as written under F2/F3; the legacy-manna reading is unchanged by the
  rent term; and a hand-split-bill fixture exercises the bill-factor site. The
  full kernel suite is green at the `v0.2.0` tag (68 testsets, 0 failures,
  re-run 2026-10-15 for this acceptance).
- Paper tests, `tests/test_pricing_design.jl`: the `[pricing]` design schema,
  its preregistration pinning, and — because the synthetic fixtures carry
  `M_int = T_int = 0` — the external account and the companion identity
  `gdp_wedge = -diff` with each door open ON THE REAL TABLE, plus the Stage B
  budget gate through `evaluate_gates`.
- The pin bump to `v0.2.0` is verified on the real pipeline:
  `scripts/reproduce.jl --check` reproduces every `matrix_5x3_v10` cell with
  zero difference.
- `pricing_bill_factor` is the single definition for the pattern; a grep for
  `p .^ ϵ` must not find a fifth site.
- Promotion checklist before any citable run: this ADR accepted ->
  `registry/closures.toml` entry -> design + `run.jl --preregister` -> run ->
  manifests -> the paper repo's Manifest pin bump as a deliberate generation
  event.
- Probe numbers (`probe17`, `probe20`) are NOT citable: no preregistration and
  no manifest (ADR-0004). They motivate this decision; they do not evidence a
  result.
- On acceptance, a one-to-two-sentence summary is appended to the kernel's
  `DECISIONS.md` (`OpenCGE@ADR-0029`).
