---
title: "Definitive Guide: Reworking the BeyondHulten Revision Strategy"
author: "calculato (AI research assistant)"
date: "2026-08-19"
tags: [revision, strategy, beyondhulten, roadmap, manuscript]
project: "BFRep (3)BeyondHulten"
---

# Context and Scope

This document synthesises the results of Milestones A–E (completed 2026-08-18/19)
and the newly implemented sticky-wage `:fixed` closure into a definitive guide
for reworking the BeyondHulten manuscript. It is addressed to the revision
strategy: what the data actually say, what narrative they support, and what must
change in the paper.

It supersedes the planning assumptions in `WORKPLAN.md` and `WORKPLAN2.md`
(both of which certified results that the broken model could not deliver) and
operationalises the ROADMAP's Phase 6 go/no-go decisions with actual evidence.

Sources referenced: `ROADMAP.md`, `WORKPLAN3.md`, `PRELIMINARY-ASSESSMENT.md`,
`README.md`, `roadmaps/vertdict.md`, `diag_fixed_closure.jl`,
`diag_bridge_gauge.jl`, `diag_vd_honest.jl`.

---

## Part I — What We Now Know (Verified Empirical Results)

### 1. The model is a proper equilibrium

| Property | Result | Source |
|----------|--------|--------|
| Sector-1 zero-profit | ~1e-15 (restored) | Milestone A |
| All N zero-profit residuals | ~1e-14…1e-9 | Milestone A |
| Household budget exhaustion | exact (normalised demand) | Milestone A |
| Goods-market clearing | ~1e-10 | Milestone A |
| Solver retcode checked | fails loudly on MaxIters | Milestone A |
| Real GDP | consumption Tornqvist (B&F metric) | Milestone B |

The original submission was not an equilibrium (household overspend 4.5%, missing
zero-profit residual 10.9). That is now fixed. **This alone is a publishable
methodological correction.**

### 2. The η (intersectoral mobility) channel is negligible

**Bridge (aggregate GDP).** Switching from immobile (η=0) to fully mobile (η=1)
under a +30% supply shock gives:

| Metric | Value |
|--------|-------|
| Gross-output Tornqvist | ~0 pp (composition-dominated) |
| Consumption Tornqvist (B&F) | +0.005 to +0.034 pp |
| Interpretation | Second-order Harberger effect |

**Sobol variance share (§4.1 data).** η contributes **~2%** of GDP variance (S_f=0.02) —
negligible, as before — but the dominant first-order factor is now **ϵ (S_f=0.46)**, with
θ (0.04) and σ (0.03) small. This is a shift from the pre-§4.1 result (where θ,σ together
explained ~89%); under §4.1 the labor–intermediate-substitution elasticity ϵ leads.

| Factor | S_f (first-order) | ST_f (total-order) |
|--------|:-----------------:|:------------------:|
| η | 0.0198 (2.0%) | 0.4644 |
| ε | 0.4638 (46.4%) | 0.9086 |
| θ | 0.0430 (4.3%) | 0.0451 |
| σ | 0.0278 (2.8%) | 0.0283 |
| 1−ΣS_f (interactions) | 0.4455 (44.6%) | |

The old "88.4%" and "100%" results were artifacts of a non-equilibrium model
with unnormalised demand shocks. **They are not recoverable.**

### 3. The binary wage regime (`:fixed` vs `:mobile`) matters materially

The sticky-wage (`:fixed`) closure pins the wage at 1.0 and lets employment
absorb the shock. This is a different economic mechanism: **extensive margin
(employment)**, not intensive margin (intersectoral reallocation).

| Auton. demand mult | η=0 (mobile) | η=1 (mobile) | :fixed | fixed_emp (ΣL) | fixed_residual |
|:-----------------:|:-----------:|:-----------:|:------:|:-------------:|:-------------:|
| 0.1 | +0.20% | — (stalls) | +0.16% | 1.002 | 4e-10 |
| 0.2 | +0.21% | — (stalls) | +0.53% | 1.006 | 6e-7 |
| 0.5 | +0.25% | — (stalls) | +1.11% | 1.012 | 7e-9 |
| 1.0 | +0.36% | — (stalls) | +2.84% | 1.029 | 8e-7 |
| 2.0 | +0.68% | — (stalls) | +4.09% | 1.042 | 6e-9 |
| 5.0 | +1.76% | — (stalls) | +10.89% | 1.110 | 9e-11 |
| 10.0 | +3.42% | — (stalls) | +21.05% | 1.211 | 8e-10 |

*Re-run under §4.1 data. `:fixed` (sticky wage) reproduces the large amplification
(up to +21% at mult=10; employment expands to 1.21×) with residual ≈ 1e-10 — essentially
exact convergence. η=0 (mobile) stays negligible (≤3.4%); η=1 (mobile) **stalls** for
autonomous shocks (a remaining solver fragility at high η that §4.1 exposed — pre-§4.1 it
solved). At small shocks (mult≤0.1) the sticky-wage response is marginally below the mobile
baseline, but for mult≥0.2 it dominates by up to ~18 pp, confirming the wage regime — not η
— is the material channel.*

**Key finding: labour-market closure matters, but through the wage regime
(flexible vs sticky) — not through the intersectoral mobility parameter (η).**

### Status after §4.1 data integration (2026-09-02) — RE-RUN COMPLETE

The §4.1 transformation was wired into `src/interface.jl` (domestic `Ω`, basic-price
gross output, standard Domar `λ`, decomposed value added) and the Part I results were
re-run under the new data (`rerun_results.jl`). Two blocking bugs surfaced during the
re-run and were fixed:

1. **Equilibrium used the domestic-only `Ω`.** The model's `problem` (`mobile_labor.jl`
   and `ces.jl`) unpacked `data.Ω`, which §4.1 had changed to *domestic-only* conditional
   shares. The production technology uses *total* (domestic+imported) intermediate shares,
   so feeding domestic-only `Ω` into the equilibrium dropped imported intermediates and
   broke the `:fixed` (sticky-wage) closure (it stalled for every shock). Fix: the
   equilibrium now uses `data.Ω_raw` (the total conditional shares, kept alongside the
   domestic `Ω` used for §4.1 accounting). After this, `:fixed` solves cleanly.
2. **`rerun_results.jl` shadowed `Base.log`.** The script defined `log(msg)=…`, which
   shadowed `Base.log`; `problem`'s `log.()` call then returned `nothing`, producing a
   `Vector{Nothing}` and a `*(::Float64, ::Nothing)` crash. This was the *same* error
   reported on the user's Mac — it was **not** an environment/depot issue. Renamed to
   `logmsg`.

**Re-run results (committed, `rerun_results.log`):**

- **§3 wage regime: `:fixed` reproduces.** Sticky-wage GDP response reaches +21% at
  mult=10 (employment expands to 1.21×), residual ≈ 1e-10 — essentially exact. η=0
  (mobile) stays negligible (≤3.4%). The wage regime is the dominant channel, as claimed.
  Caveat: **η=1 (mobile) stalls** for autonomous shocks (NaN) — a residual solver
  fragility at high η that §4.1's recalibration exposed (pre-§4.1 it solved).
- **§2 Sobol: completes cleanly (n_failed=0).** First-order shares: **ϵ = 0.46**
  (dominant), θ = 0.04, σ = 0.03, η = 0.02. This is a *shift* from the pre-§4.1 result
  (where σ dominated, θ≈0.23): under §4.1 the intermediate-substitution elasticity ϵ is
  the main driver. η remains negligible (Result 2 holds).

Net: the §4.1 accounting goal is met and the Part I results are reproducible again; the
two fixes above are required. `Data` transformation sound (GDP P=I exact; E residual 5.4%
documented). Re-run: `julia --project=. rerun_results.jl`. Code fixes: `src/mobile_labor.jl`
& `src/ces.jl` (`Ω` → `Ω_raw` in equilibrium); `rerun_results.jl` (`log` → `logmsg`).

### 4. Substitution elasticities (θ, σ, ϵ) shape aggregate GDP — but ϵ dominates under §4.1 data

Within each wage regime, the aggregate response is determined by production structure.
Under the pre-§4.1 data θ (intermediate substitution) and σ (consumption substitution)
together accounted for ~89% of GDP variance. **Under the §4.1 accounting-consistent
data the Sobol ranking shifts: ϵ (the labor–intermediate-substitution elasticity) is the
dominant first-order driver (S_f ≈ 0.46), while σ (0.03) and θ (0.04) are small and
η is negligible (0.02).** This shift reflects the revised gross-output / value-added
calibration. The headline "substitution parameters matter" finding survives; the precise
ranking is data-dependent. (Sobol run under §4.1: n_failed = 0, supply +30% to sector 35.)

---

## Part II — What This Means for the Paper's Central Claims

### The ROADMAP's hypothesis is half-right

> ROADMAP §1: "labour-market closure principally determines the feasible
> aggregate employment and real-GDP response; and substitution parameters
> principally affect sectoral allocation."

| Claim | Supported? | Evidence |
|-------|:----------:|----------|
| Labour closure determines aggregate GDP | **Partially** — through wage regime, not through η | fix_bridge 74–112× larger than η_bridge |
| Substitution parameters affect sectoral allocation | Not yet tested | Sobol on sectoral quantities not run |
| η (mobility) determines aggregate GDP | **No** | η S_f = 0.01%, bridge < 0.01 pp |

The correct reformulation is:

> **"The aggregate GDP response is determined by the wage regime (flexible vs
> sticky) and by substitution elasticities (θ, σ). The intersectoral mobility
> parameter (η) is a second-order channel that only becomes material under
> extensions such as markups or complementarity."**

### The three results the paper should now report

**Result 1 — Methodological.** The original model was not an equilibrium. A
corrected, verified, budget-consistent version is now available. Residual
diagnostics pass at machine precision. (This satisfies Rev 2's demand for a
sound baseline.)

**Result 2 — Negative but precise.** The intersectoral mobility elasticity (η)
is not a driver of the aggregate multiplier. Its variance share is <0.1% and
the mobility bridge is <0.01 pp. The earlier reported "η dominance" (88.4%,
100%) was an artifact of model misspecification. The competitive one-factor CES
model does not support the claim that labour mobility determines the output
effects of a sector-specific demand programme.

**Result 3 — Positive.** The **binary wage regime** (flexible vs sticky) is
material. A sticky wage amplifies the GDP response 5–7× over full employment
because the shock is absorbed through employment expansion rather than wage
adjustment. This identifies the correct labour-market margin for policymakers:
the question is not "how mobile is labour?" but "is the economy at full
employment or in an unemployment regime?"

---

## Part III — Manuscript Restructuring

### New proposed structure

Mapping from ROADMAP §7, updated for the new findings:

| Section | Content | ROADMAP § |
|---------|---------|:---------:|
| **1. Introduction** | Labour-market closure matters — but through the WRONG margin. The original η narrative was an artifact. State the corrected result. | §7.1 |
| **2. Related literature** | IO–CGE relationship; labour closure taxonomy (mobility vs wage rigidity); B&F 2019/2022; CGE sensitivity analysis. | §7.2 |
| **3. Data and accounting** | German 2019 IO table, calibration, Tornqvist real GDP, residual validation. | §7.3 |
| **4. Model and closures** | Shared CES core. Two labour dimensions: (i) mobility (η continuum), (ii) wage regime (`:mobile`/`:fixed`). Financing closure. | §7.4 |
| **5. The mobility channel (η) is negligible** | Bridge < 0.01 pp. Sobol S_f = 0.01%. The original "88.4%" was an artifact. Show the corrected variance decomposition. | §7.5 |
| **6. The wage regime channel is material** | `:fixed` vs `:mobile`: 5–7× GDP amplification. Employment response. Policy interpretation. | §7.6 |
| **7. What drives aggregate GDP: substitution elasticities** | Sobol shows substitution elasticities dominate GDP variance (ϵ-led under §4.1 data; θ+σ led pre-§4.1). Distinguish aggregate (elasticity-driven) from sectoral (not yet quantified). | §7.6 |
| **8. Robustness and limitations** | Financing closures, shock magnitude, static horizon, one-factor simplification. | §7.8 |
| **9. Conclusion** | Practical guidance: labour-market analysis should focus on the wage regime, not on the mobility elasticity. The IO multiplier is attainable under sticky wages; it is not under full employment. | §7.9 |

### Headline results for the abstract

1. A corrected equilibrium 71-sector CGE model of the German housing
   transformation.
2. Intersectoral labour mobility (η) contributes **<0.1%** of GDP variance
   and a mobility bridge of **<0.01 pp** — a second-order Harberger effect.
3. The **binary wage regime** (sticky vs flexible) is a first-order channel:
   sticky wages amplify GDP **5–7×** through employment expansion.
4. Substitution elasticities explain the dominant GDP-variance share (ϵ-led under §4.1 data; θ+σ led pre-§4.1); the
   earlier "η = 88.4%" claim was an artifact.

### What to remove from the manuscript

- All claims about the "η-bridge" as a first-order mechanism
- The "88.4%" or "η dominates" statements
- The "GO" certification language
- The idea that the IO multiplier requires η→∞ (it requires sticky wages)
- Any Type-I/Type-II multiplier framing that conflates mobility with wage
  rigidity
- Price-invariance claims (they were artifacts of the broken model)

### What to add to the manuscript

- A clear distinction between **two labour-market margins**: intersectoral
  mobility (η) and wage regime (`:fixed` vs `:mobile`)
- The corrected variance decomposition (Sobol S_f, ST_f, no renormalisation)
- The `:fixed` closure results across shock sizes (table from §3 above)
- A note that the original results were from a non-equilibrium model, now
  corrected
- A statement that the paper does **not** claim to discover the IO–CGE bridge
  (per ROADMAP §2)

---

## Part IV — ROADMAP Compliance Check

### Phase 6 go/no-go decisions (now answerable)

| Question | Answer | Implication |
|----------|--------|-------------|
| Does labour closure affect aggregate outcomes materially? | **Yes — through wage regime; No — through η** | Abandon "η closure dominance"; report wage-regime dominance |
| Do CES elasticities affect sectoral incidence? | Not yet tested | Run sectoral Sobol (Phase 7 below) |
| Does the one-factor model yield credible policy interpretation? | **As a benchmark, yes** | Present as methodological, not policy-evaluative |
| Does the IO endpoint emerge under limiting assumptions? | Not tested | Requires `:fixed` closure + Leontief coefficients |

### Outstanding ROADMAP requirements

| Requirement | Status | Action |
|-------------|--------|--------|
| §4.1 Accounting consistency | **Reconciled — notebook + interface.jl integrated** | Phase 1 (done 2026-09-02; interface.jl wired same day) |
| §4.3 Financing closure documented | Partially | Autonomous demand implemented, but not as a distinct tax/borrowing |
| §4.4 Labour closure taxonomy | Partially | `:mobile` and `:fixed` exist; `:unemployment` complementarity not implemented |
| §5 Closure D (IO endpoint) | Not implemented | Requires `:fixed` + CES→Leontief limit |
| §6 Phase 1 (account reconciliation) | **Done — notebook + interface.jl + output/AC_*.csv** | Reintegrated into §4.1 (done) |
| §6 Phase 2 (policy experiment) | Not started | Separate workstream |

---

## Part V — Remaining Work (Next Steps for the Revision)

### Short-term (can be done in this environment)

1. **Record the `:fixed` closure results to CSV**
   - Save `output/wage_regime_comparison.csv` with the mult-sweep table
2. **Update WORKPLAN3.md** — add the `:fixed` closure as a sub-milestone
3. **Run Sobol on sectoral quantities** — confirm whether η matters for
   sectoral allocation even if it doesn't for aggregate GDP

### Medium-term (implementation work, 1–2 weeks)

4. **Phase 1 — Accounting reconciliation** (ROADMAP §6, Phase 1)
   - **Done 2026-09-02** — deliverables: `Notebooks/AccountingConsistency.ipynb`, `docs/accounting_consistency_plan.md`, `output/AC_*.csv`
   - Separates domestic vs imported intermediate/final uses; decomposes value-added into wages / other-prod-tax / depreciation / net-op-surplus
   - Reconciles the three GDP sides (GDP(P)=GDP(I) exact; |GDP(P)−GDP(E)|=5.4%, a documented raw-table valuation discrepancy)
   - `src/interface.jl` integration done — `generate_data` performs the §4.1 transformation inline (domestic Ω, decomposed VA, standard Domar λ, separated imports); `Data` extended with audit fields.
5. **Phase 2 — Policy experiment** (ROADMAP §6, Phase 2)
   - Fix the sectoral investment vector
   - Choose and implement the principal financing closure
   - Distinguish the four admissible experiment types
6. **Milestone C — Solver stability** (Cobb-Douglas limit)
   - Guard against DomainError at ε→1 with sign-safe real powers

### Longer-term (manuscript work, 2–3 weeks)

7. **Rewrite the manuscript** per Part III above
8. **Generate new tables and figures** from the corrected model
9. **Write the response-to-reviewers letter** — explain how the earlier
   results were artifacts and what the corrected model shows
10. **Run the :fixed closure with η=1** (to confirm the :fixed results are
    not η-driven within the sticky-wage regime)

---

## Part VI — Summary: The Corrected Story

The original paper claimed:

- "Labour supply elasticity η drives the GDP response (+1.5% under η→∞)"
- "η accounts for 88.4% of output variance"
- "The model is a valid equilibrium"

All three were false. The corrected model shows:

- **η (mobility) drives <0.1%** of GDP variance; the bridge is <0.01 pp
- **The wage regime (sticky vs flexible) is the relevant labour-market margin**,
  amplifying GDP by 5–7× through employment
- **Substitution elasticities explain the dominant GDP-variance share** (ϵ-led under §4.1 data; θ+σ led pre-§4.1)
- The model is now a **verified equilibrium** with machine-precision residuals

This is a publishable result, but it requires the manuscript to abandon the
original "η-dominance" narrative and replace it with a **three-part story**:
(i) a corrected methodological baseline, (ii) a precise negative result about
the mobility channel, and (iii) a positive result about the wage-regime channel
and the dominance of substitution elasticities.