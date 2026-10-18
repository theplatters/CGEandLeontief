---
title: "IDEA-0005: pricing-door workstream --- handover of what is built, decided, measured and next"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-10-15"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [ideas, pricing, handover, verdoorn, markup, adr-annexe, phase-7]
last-updated: October 2026
---

**Version 1** (October 2026)

Status: `scoped`.

# What this note is, and what it deliberately is not

A HANDOVER of the pricing-door workstream (ROADMAP Phase 7): what was built,
what was decided, what was measured, what is still open, and in what order the
remaining work has to happen. Written for whoever picks this up next --- a
collaborator, a later session, or the author after a break.

Three neighbours already own parts of this material, and this note restates
none of them:

| Owns | Document |
| --- | --- |
| The DECISIONS | `docs/decisions/ADR-0029-pricing-door.md` (accepted), `docs/decisions/ADR-0030-distributed-profit-income.md` (accepted) |
| The PLAN (coding + writing steps) | `ROADMAP.md` Phase 7, Routes 1 and 2 |
| The DOOR RANKING (why these and not others) | `docs/ideas/IDEA-0002-five-doors-demand-sensitive-prices.md` |

This note is the map over the three: state, gaps, sequence, traps.

Register caveat, as the register requires: an idea carries no paper weight and
no `registry/` entry. The commitments live in the ADRs and the registry; this
note only points at them. The status vocabulary has no "handover" state ---
`scoped` is the nearest fit, because the route it hands over is specified but
UNEXECUTED at the batch level (the 42 cells are preregistered and not run).

**Register housekeeping this note implies.** IDEA-0002 now satisfies this
register's own definition of `promoted` --- "an ADR and a registry entry
exist": `ADR-0029`/`ADR-0030` are accepted and `[pricing.P1]` / `[pricing.P2]`
exist in `registry/closures.toml`. It should be raised from `scoped` with
pointers to the ADRs and the registry entries. Left to the operator.

# The obstruction this axis exists to remove

Foundations I (labour) and II (financing) are demand-side. Under a demand-only
shock ($A \equiv 1$) the zero-profit block is homogeneous of degree one in
$(p, w)$ and contains no demand term, so it pins prices at one in both gauges
and the real wage at its anchor; the quantity block then reduces to an affine
multiplier system in which none of $(\theta, \epsilon, \sigma, \eta_s)$
appears. That is Lemma 1 and Lemma 2 of `paper/equivalence.tex`, and it is why
ALPHA $\equiv$ BETA and GAMMA $\equiv$ DELTA hold as propositions.

Closure variety inside the demand-shock story therefore requires demand to
enter the PRICE block. The kernel's cost line is the door map ---

    cost = _ces_unit_cost(A_eff, factor_share, w, intermediate_price, eps)

--- and exactly three of its arguments are demand-injectable: `A_eff`
(capacity/Verdoorn, the door taken), `w` (sectoral labour markets, ADR-0022,
already executed) and `ip` (import prices, not taken). A markup enters the same
line as a multiplier but creates a rent, which is a different KIND of change
because it needs an income decision. The substitution elasticities are NOT a
door (Lemma 2 makes the aggregate invariant to them at every value).

Cross-cutting consequence, and it is load-bearing: any pricing door breaks
$p \equiv 1$, so it DISSOLVES GAMMA $\equiv$ DELTA and makes BETA $\neq$ ALPHA
reappear under demand-only shocks. ALPHA $\equiv$ BETA survives only on the
slice "A $\equiv$ 1 AND price-invariant technology".

# State of the build (verified 2026-10-15)

| Item | Id / path | State |
| --- | --- | --- |
| Pricing door decision | `ADR-0029` | **accepted** 2026-10-15 with OpenCGE v0.2.0 |
| Distributed profit (Stage B) | `ADR-0030` | **accepted** 2026-10-15 |
| Kernel | OpenCGE **v0.2.0**, tree `e119f8d778632b10f4312c9f0fa37b7d827ed721` | pinned in the paper repo's `Manifest.toml` |
| Kernel record | `OpenCGE@ADR-0029`, `OpenCGE@ADR-0030`, `OpenCGE@DE-0012` | in the kernel's `DECISIONS.md` |
| Registry | `[pricing.P1]` capacity, `[pricing.P2]` markup | status **`tested`** (not yet `validated`) |
| Designs | `pricing_capacity_v1`, `pricing_markup_v1`, `pricing_profit_v1`, `pricing_mobile_v1` | in `experiments/designs/` |
| Preregistration | 42 scenario rows | **all `planned`; none executed** |
| Kernel tests | `OpenCGE@tests/test_pricing_closure.jl` | pass (full suite 454 tests, 58 testsets) |
| Paper tests | `tests/test_pricing_design.jl` | asserts the external account on the REAL table |
| Manuscript | App. B `model-closures-proofs` (merged), appendices relettered A--H, `revision/proposals/` | done (merged 2026-10-18) |

Kernel commits that make up v0.2.0, oldest first, all in `../opencge`:
`1be5809` hooks + `Model.pricing` (no behaviour change); `c40cc4d` `P1`
`CapacityPricing` solved INSIDE the residual; `7dfe555` `P2` `MarkupPricing`
Stage A with the FOC at marginal cost; `a08e40e` `P2` Stage B distributed
profit; `1265a70` registry id map; `7cf8ed0` the ADR/DE lines; `75116d6`
release.

Door parameters: `CapacityPricing(delta)` with `A_eff_i = A_i (y_i/lambda_i)^(-delta)`.
`MarkupPricing(kappa; s_pi)` with `p_i = mu_i c_i`, `mu_i = 1 + kappa (y_i/lambda_i - 1)`;
rent `Pi = sum_i (1 - 1/mu_i) p_i y_i`; `s_pi = 1` is Stage A EXACTLY, `kappa = 0`
nests the baseline, `delta = 0` nests the baseline.

Design schema that carries them (per design, overridable per cell):

```{.toml}
[pricing]
kind = "none" | "capacity" | "markup"
delta = 0.0          # capacity
kappa = 0.0          # markup
s_pi  = 1.0          # markup retention; 1.0 = Stage A
multistart_tol = 1e-6
continuation_steps = 4
```

# The single most important fact in this handover

**Nothing in this workstream is citable yet.** The 42 cells are preregistered
and unrun, so there is no manifest, so P1/P2 are `tested` (the listed tests
pass) and not `validated` (a recorded run manifest exists). The probe numbers
(`probe17`, `probe20` --- 75 cells on the real calibration) MOTIVATE the
decision and must never be cited as results: no preregistration, no manifest
(ADR-0004).

Everything downstream --- Section 6.3/6.4 prose, the price-response bracket
figure, the tables --- starts by minting the generation.

# What remains, in order, with its gate

1. **Run the batch.** `experiments/run.jl --design pricing_capacity_v1` (then
   `pricing_markup_v1`, `pricing_profit_v1`, `pricing_mobile_v1`). One scenario
   may be `running` at a time (ADR-0007); `scripts/check_repo.jl` before and
   after, 0 violations required. PRE-FLIGHT FIRST: the two earlier generations
   (v7, v8) each burned 33 runs on plumbing --- `include` `run.jl` from a
   scratch script and drive `load_design` -> `build_reference` ->
   `programme_vectors` -> `build_cell_model` -> `solve_cell` ->
   `evaluate_gates` over ONE cell of every type before opening the design.
2. **Mint the generation and flip the status** --- P1/P2 from `tested` to
   `validated`, with the run ids recorded in `registry/scenarios.csv`.
3. **Write Sections 6.3/6.4 from the manifests**, numbers by run id, never
   hand-typed. The presentation decisions are ADR-0029 #10: the KV/markup pair
   as ONE family with two signs and a residual, presented as a
   PRICE-RESPONSE BRACKET (a result, not two robustness arms); a SIGNED price
   statistic (a lone `max|p-1|` hides the financing-conditional sign); the
   financing-conditional reading; and the scope clause.
4. **Housekeeping owed:** a `docs/log/2026-10.md` entry for the 2026-10-15
   session (the register and the log are separate obligations; the log entry
   for the promotion/naming/registration pass is not written); IDEA-0002 raised
   to `promoted` with pointers (see above).

# Open gates and live questions

Carried in `registry/closures.toml` as `open_gates`; listed here with the
decision each one needs.

| # | Open gate | What it needs |
| --- | --- | --- |
| 1 | `delta < 0` is increasing returns (positive feedback) | every priced cell carries the ADR-0029 #11 uniqueness gate: `pricing_admissibility` multi-start spread `< 1e-6` plus a full-rank residual Jacobian. Non-convergence at extreme `delta` is a documented FINDING, not a bug |
| 2 | The ADR-0017 scale guard runs on a NoPricing clone at the reference | exact there (`A_eff(lambda) = A`); determinacy AWAY from the reference rests on the local Jacobian check --- a stated limitation, not a computation |
| 3 | The markup is symmetric about the baseline | a contracting sector carries `mu < 1` (a DISCOUNT, hence negative sectoral rent). Accepted as stated (ADR-0029 #11), guarded by `mu > 0` (`PricingFeasibilityError`), min/max `mu` recorded per cell. The MANUSCRIPT must state the discount reading rather than hide it |
| 4 | Retention needs injections in the fixed-wage gauge | with none the system has no positive solution (demand shortfall: `(mu-1) c y + (s + m - s m) w L = p(1-m)(gG+exo) + pX`). The application carries injections; say so |
| 5 | Stage B v1 leaves the rent untaxed (`tau` on wage income only) | the untaxed rent UNDERSTATES the F2 fiscal response; the v1 tables must say so. The two-base tax variant is a robustness row, not run |
| 6 | `s_pi in {0.5, s}` with `s = 0.1199` | crosses kappa in `{0.25, 0.5}` on F1-F3 (ADR-0030 #6); Stage A is `s_pi = 1` of the same generation |
| 7 | Stage C (two-class household) | parked |
| 8 | The measured SIGN is financing-conditional | at `delta = -0.5` the GDP deflator falls only under F3 (0.9853) and rises under F1 (1.0035) / F2 (1.0013) --- "Verdoorn lowers prices" is an F3 statement, not a property of the door |
| 9 | The collapse count depends on the beta on display | 7 distinct economies under the single-wage beta, 9 under the sectoral column Table 4 shows. The abstract's collapse sentence, Table 4's caption and the "why only N economies" box must state the SAME number for the SAME beta |

# Doors not taken, and why they stay open

- **Import prices / terms of trade** (`ip`): the best accounting fit to the
  existing external machinery and no domestic income closure, but the most
  assumption-laden --- a distinct `p^M_i`, an isoelastic foreign supply curve,
  a MANDATORY Armington aggregation choice, an import-content assumption for
  final demand, a two-price CPI, an export treatment, and a revaluation of the
  external-account identity that the repo already flags as easy to get wrong by
  exactly one term. Spec: `IDEA-0002` and the door reference in the project
  skill.
- **Capital / fixed factor with a rental rate**: the structural
  MICROFOUNDATION of `delta`, but the largest surgery (second factor, rental
  schedule, rental-income closure) and it collides with the stated one-factor
  design boundary. Present it in the open-doors discussion, never as a
  headline: the paper's contribution is a map of closures, not a new model
  class.

# Traps that cost time here (read before touching the code)

- **A synthetic fixture cannot see the bill-factor sites.** The FOC pattern
  `p^eps A^(eps-1)` lives at FOUR places, not two: the wage-FOC inversion, the
  INNER `p^eps` of `intermediary_demand`, and BOTH CES bill factors (the
  canary's `k` and its hand-copied duplicate in `gdp_components`). The leaks
  those bill factors value (`M_int`, `T_int`) are ZERO in anything the compact
  `Data` constructor builds, so every identity stays green at machine precision
  while the real calibration is wrong. Cure: one shared `pricing_bill_factor`,
  plus a fixture whose purchaser-price bill is split BY HAND
  (`A_bill + M_int + T_int = (1-fs) lambda`), plus the real-table assertion in
  `tests/test_pricing_design.jl`.
- **`p` has two roles.** The purchaser price stays `p`; the FOC multiplier
  becomes marginal cost `c = p/mu`. Getting this wrong leaves the price
  response right and every labour demand inflated by `mu` --- the dangerous
  failure mode, because the headline numbers look plausible.
- **The uniform-markup canary is the CASCADE**, `p = mu_bar^((I-B)^(-1) 1)`
  with `B = diag(1-fs) Omega_raw`, about 1.8 on the German table --- NOT
  `p = mu_bar`. A canary written as `p = mu_bar` rejects a CORRECT
  implementation.
- **The external-account gate has two documented readings.** With zero legacy
  manna (the matrix designs) the identity holds AS WRITTEN with the rent term
  booked; on a legacy-manna design the gap legitimately reads `p (A + G)`.
  Assert against the documented value, never against zero --- an `approx 0`
  assertion there reports the NEW term as the culprit.
- **The outer fixed point is not the promoted form.** Iterating
  `A <- A (y/lambda)^(-delta)` around the solver diverges at damping 0.6 for
  `delta = -0.5` under Cobb-Douglas (`OpenCGE@DE-0012`). `A_eff` is solved
  inside the residual instead. The probe loop stays as the cheap instrument.
- **`pricing_technology` needs the generic `AbstractPricingClosure` fallback**
  returning `A` unchanged, because a markup alters the PRICE, not the
  technology: a per-type definition silently omits it and every markup solve
  dies inside the demand block.
- Environment: `TMPDIR=/tmp JULIA_DEPOT_PATH=/opt/julia-depot`, or `Pkg` dies
  on a non-existent temp dir; `Pkg.instantiate()` takes about 10 minutes
  (SciML precompilation) and leaves the committed `Manifest.toml` untouched.
  Verify a kernel checkout by `git rev-parse <tag>^{tree}` against the paper
  repo's Manifest `git-tree-sha1` BEFORE trusting a number from it.
- Background shells report the WRAPPER's exit code: take a suite's status from
  its log file, not from the completion notice.

# Reproduction one-liners

```{.bash}
# kernel (../opencge): suite
cd ~/Git/BFRep/opencge && TMPDIR=/tmp JULIA_DEPOT_PATH=/opt/julia-depot \
  julia --project=. tests/runtests.jl

# paper repo: board, gate, tests
cd "~/Git/BFRep/(3)BeyondHulten" && TMPDIR=/tmp JULIA_DEPOT_PATH=/opt/julia-depot \
  julia --project=. scripts/status.jl
julia --project=. scripts/check_repo.jl        # 0 violations required
julia --project=. tests/runtests.jl

# the batch (after the pre-flight)
julia --project=. experiments/run.jl --design pricing_capacity_v1
```

# Pointers

- Decisions: `docs/decisions/ADR-0029-pricing-door.md`,
  `docs/decisions/ADR-0030-distributed-profit-income.md`
- Plan: `ROADMAP.md` Phase 7 (Routes 1 and 2, shared decisions, outputs)
- Ranking and per-door specs: `docs/ideas/IDEA-0002-…`,
  `references/demand-sensitive-price-doors.md` (project skill)
- Registry: `registry/closures.toml` `[pricing.P1]` / `[pricing.P2]`;
  `registry/preregistration.toml`
- Designs: `experiments/designs/pricing_{capacity,markup,profit,mobile}_v1.toml`
- Evidence (not citable): `experiments/probes/probe17_gamma_capacity.jl`,
  `experiments/probes/probe20_pricing_door.jl`
- Kernel: OpenCGE v0.2.0, `../opencge`, `src/closures/pricing/pricing.jl`,
  `DECISIONS.md`, `DE-0012`
- Log: `docs/log/2026-10.md` (2026-10-14 probe20 entry; the 2026-10-15 entry
  is owed)
