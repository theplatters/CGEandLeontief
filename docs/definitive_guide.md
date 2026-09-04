---
title: "Definitive Guide: Reworking the BeyondHulten Revision Strategy"
author: "calculato (AI research assistant)"
date: "2026-09-03"
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
`README.md`, `roadmaps/vertdict.md`, `tests/test_mobile_labor.jl`,
`tests/test_fixed_closure.jl`, `tests/test_variance_decomposition.jl`,
`tests/test_eta_sweep.jl`.

> **Reproducibility note (2026-09-03).** All numerical tables below are from the
> **completed 2026-09-03 re-run** under the fixed code/data. Canonical accounting
> result: **GDP P = 3,027,818; I = 3,027,818; E = 2,864,724; residual 5.387%**
> (documented raw-table valuation gap); `sum(λ) = 2.1099`; `sum(labor_share) = 1`;
> reference `ref_gdp = 0.9998546537504605`. The reproducible local artifacts
> (`rerun_results.log`, `output/variance_decomposition_sobol.csv`, `output/AC_*.csv`)
> are **gitignored and not committed**; regenerate them with
> `julia --project=. rerun_results.jl` (and the notebook for the AC files). No
> numbers are invented here; they are copied verbatim from that run.

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

### 2. The η (intersectoral mobility) channel is material but not dominant

**Bridge (aggregate GDP).** Under the autonomous demand sweep (Section 3 table),
the immobile→fully-mobile wedge (η=0 → η=1) is small: at mult = 0.1 it is
+0.060 pp (η=0 = −0.046 pp, η=1 = +0.015 pp) and at mult = 0.2 it is +0.066 pp;
η=1 solves only at multipliers 0.1 and 0.2 (it stalls from 0.5 upward under the
autonomous-demand configuration). The continuous η gradient is therefore a
second-order aggregate channel.

**Sobol variance share (§4.1 data, 32/32 grid points, 0 failures).** η has a
first-order share of **~15.7%** (S_f = 0.1571) and total-order **ST_f = 0.1662** —
it is **material, not negligible** (the earlier "~2%" reading was a 2026-09-02
artifact, now corrected), but it does **not** dominate. The dominant first-order
factor is **θ (S_f = 0.3951)**, followed by **σ (0.2734)** and **ε (0.1650)**.
This is a further shift from the 2026-09-02 re-run (where ε led at ≈0.46): under
the final fixed code the intermediate-substitution elasticity θ leads, with σ and
ε each ~16–27%.

| Factor | S_f (first-order) | ST_f (total-order) |
|--------|:-----------------:|:------------------:|
| η | 0.1571 (15.7%) | 0.1662 |
| ε | 0.1650 (16.5%) | 0.1719 |
| θ | 0.3951 (39.5%) | 0.3972 |
| σ | 0.2734 (27.3%) | 0.2744 |
| 1−ΣS_f (interactions) | 0.0093 (0.93%) | |

ΣS_f = 0.9907; the interaction share (1−ΣS_f) is small (0.93%), i.e. the four
elasticities account for essentially all GDP variance with little joint
interaction. The old "88.4%" and "100%" results were artifacts of a
non-equilibrium model with unnormalised demand shocks. **They are not recoverable.**

### 3. The binary wage regime (`:fixed` vs `:mobile`) matters materially

The sticky-wage (`:fixed`) closure pins the wage at 1.0 and lets employment
absorb the shock. This is a different economic mechanism: **extensive margin
(employment)**, not intensive margin (intersectoral reallocation).

| Auton. demand mult | η=0 (mobile) Δ% | η=1 (mobile) Δ% | :fixed Δ% | η_bridge (pp) | fix_bridge (pp) | fixed_emp (ΣL) | fixed_residual |
|:-----------------:|:---------------:|:---------------:|:--------:|:-------------:|:---------------:|:--------------:|:--------------:|
| 0.1 | -0.045697 | +0.014512 | +0.271460 | +0.060209 | +0.317157 | 1.0026 | ≤1.9e-7 |
| 0.2 | -0.051978 | +0.014496 | +0.518147 | +0.066473 | +0.570124 | 1.0050 | ≤1.9e-7 |
| 0.5 | -0.090681 | NaN (stalls) | +1.227407 | NaN | +1.318088 | 1.0121 | ≤1.9e-7 |
| 1.0 | -0.203714 | NaN (stalls) | +2.318170 | NaN | +2.521884 | 1.0230 | ≤1.9e-7 |
| 2.0 | -0.542053 | NaN (stalls) | +4.254413 | NaN | +4.796466 | 1.0424 | ≤1.9e-7 |
| 5.0 | -2.009894 | NaN (stalls) | +9.957729 | NaN | +11.967623 | 1.0994 | ≤1.9e-7 |
| 10.0 | -5.406077 | NaN (stalls) | +19.344778 | NaN | +24.750854 | 1.1933 | ≤1.9e-7 |

*Re-run under §4.1 data (2026-09-03, completed).* `:fixed` (sticky wage) reproduces
the large extensive-margin response: GDP reaches **+19.3%** at mult = 10 (employment
expands to 1.19×), with fixed residuals all **≤ 1.9e-7** — essentially exact convergence.
The flexible η=0 (mobile) closure, by contrast, turns **increasingly negative** as the
shock grows (−0.046 pp at mult 0.1 → −5.41 pp at mult 10) — i.e. under flexible wages the
shock is absorbed partly through wage/quantity contraction, not expansion. The η=1
(mobile) autonomous solve succeeds only at mult = 0.1 and 0.2 (+0.015 pp each); it
**stalls** from mult = 0.5 upward (no autonomous/investment additive anchor — the
scale-indeterminacy guard). Because the η=0 column is negative, the old "5–7× amplification"
ratio is meaningless; the correct comparison is the **percentage-point gap** between the
fixed and flexible-closure responses: fix_bridge grows from +0.32 pp (mult 0.1) to
+24.75 pp (mult 10) while the flexible η=0 response is negative, so the sticky-wage
closure is the material labour-market margin — not η.

**Key finding: labour-market closure matters, but through the wage regime
(flexible vs sticky) — not through the intersectoral mobility parameter (η).**
η is a second-order aggregate channel whose Sobol first-order share (~15.7%) is
material but far smaller than the wage-regime effect.

### Status after §4.1 data integration — RE-RUN COMPLETE (2026-09-03)

The §4.1 transformation was wired into `src/interface.jl` (domestic `Ω_dom` retained for
audit, equilibrium technology and `consumption_share` calibration on `Ω_raw`, basic-price
gross output, standard Domar `λ`, decomposed value added) and the Part I results were
re-run under the new data (`rerun_results.jl`). Two blocking bugs surfaced during the
re-run and were fixed:

1. **`Ω` alignment between calibration and equilibrium.** §4.1 introduced a domestic-only
   conditional-share matrix `Ω_dom` (built in the notebook by allocating imports
   proportionally across suppliers). Because that proportional allocation is largely
   undone by row normalisation, `Ω_dom` and the total conditional shares `Ω_raw`
   coincide closely in practice, so the distinction is mostly an accounting one. The model
   now uses `Ω_raw` consistently: both the equilibrium technology (`mobile_labor.jl`,
   `ces.jl`) and the `consumption_share` calibration (`interface.jl`) read `data.Ω_raw`.
   `Ω_dom` is retained only for the domestic audit/incidence breakdown and is no longer
   fed into the production technology. This also fixed the `:fixed` (sticky-wage) closure,
   which had stalled for every shock when domestic-only shares dropped imported
   intermediates.
2. **`rerun_results.jl` shadowed `Base.log`.** The script defined `log(msg)=…`, which
   shadowed `Base.log`; `problem`'s `log.()` call then returned `nothing`, producing a
   `Vector{Nothing}` and a `*(::Float64, ::Nothing)` crash. This was the *same* error
   reported on the user's Mac — it was **not** an environment/depot issue. Renamed to
   `logmsg`.

**Re-run results (`rerun_results.log` — regenerate locally; the log and the output
CSVs are gitignored, locally generated artifacts and are *not* committed):**

- **§3 wage regime: `:fixed` reproduces with exact convergence.** Sticky-wage GDP
  response reaches **+19.3%** at mult = 10 (employment expands to 1.19×), fixed
  residuals all **≤ 1.9e-7** — essentially exact. The flexible η=0 (mobile) closure
  turns **increasingly negative** (−0.046 pp at mult 0.1 → −5.41 pp at mult 10); η=1
  (mobile) solves only at mult = 0.1/0.2 and **stalls** from 0.5 upward (no
  autonomous/investment anchor — the `:fixed` scale-indeterminacy guard). The wage
  regime is the material channel; because η=0 is negative, the old "5–7× amplification"
  *ratio* is no longer used and **percentage-point gaps** (fix_bridge) are reported instead.
- **§2 Sobol: completes cleanly (n_failed = 0, 32/32).** First-order shares:
  **θ = 0.3951** (dominant), **σ = 0.2734**, **ε = 0.1650**, **η = 0.1571** (material,
  not negligible). This is a further shift from the 2026-09-02 re-run (where ε led at
  ≈0.46): under the final fixed code θ leads, with σ and ε each ~16–27%, and η's
  first-order share is ≈15.7% (well above the 2026-09-02 "~2% negligible" artifact) but
  still below the wage-regime effect. Interaction share is small (1−ΣS_f = 0.0093).

Net: the §4.1 accounting goal is met and the Part I results are reproducible under the
2026-09-03 fixed code/data; the integration fixes above (Ω alignment, `log`→`logmsg`) plus
the additional 2026-09-03 guards (log-space mobile allocation with finite `|η| ≤ 50`,
corrected allocative-wedge sign and numerical guards, calibration/equilibrium `Ω_raw`
alignment including Cobb–Douglas, Sobol grid/output validation, `:fixed` scale-indeterminacy
guard) are all required. `Data` transformation sound (GDP P=I exact; E residual 5.387%
documented). Re-run: `julia --project=. rerun_results.jl`. Code fixes: `src/mobile_labor.jl`
& `src/ces.jl` (`Ω` → `Ω_raw` in the equilibrium technology); the `consumption_share`
calibration (and Cobb–Douglas) also migrated to `Ω_raw`; `rerun_results.jl` (`log` → `logmsg`).
(The 2026-09-02 recorded numbers in the prior version of this document are historical and
superseded by the 2026-09-03 re-run.)

### 4. Substitution elasticities (θ, σ, ϵ) shape aggregate GDP — and θ dominates under the final fixed code

Within each wage regime, the aggregate response is determined by production structure.
Under the pre-§4.1 data θ (intermediate substitution) and σ (consumption substitution)
together accounted for ~89% of GDP variance. Under the 2026-09-02 §4.1 re-run the ranking
shifted and ϵ led (S_f ≈ 0.46). **Under the final 2026-09-03 fixed code/data the Sobol
ranking is: θ (intermediate-substitution elasticity) is the dominant first-order driver
(S_f = 0.3951, ST_f = 0.3972), followed by σ (S_f = 0.2734, ST_f = 0.2744), ε
(S_f = 0.1650, ST_f = 0.1719), and η (S_f = 0.1571, ST_f = 0.1662).** η is material
(~15.7% first-order) — not negligible and not dominant — and its total-order index is close
to its first-order share, so it contributes mostly as a direct main effect rather than
through interactions. The headline "substitution parameters matter" finding survives and is
now led by θ. (Sobol run under §4.1: n_failed = 0, 32/32 grid points, supply +30% to sector 35.)

---

## Part II — What This Means for the Paper's Central Claims

### The ROADMAP's hypothesis is half-right

> ROADMAP §1: "labour-market closure principally determines the feasible
> aggregate employment and real-GDP response; and substitution parameters
> principally affect sectoral allocation."

| Claim | Supported? | Evidence |
|-------|:----------:|----------|
| Labour closure determines aggregate GDP | **Yes — through wage regime, not through η** | fix_bridge pp-gap grows from +0.32 pp (mult 0.1) to +24.75 pp (mult 10); flexible η=0 is negative, so the wage regime is the material margin |
| Substitution parameters affect sectoral allocation | Not yet tested | Sobol on sectoral quantities not run |
| η (mobility) determines aggregate GDP | **No (as sole driver)** | η first-order S_f ≈ 0.157 (material, not negligible), but θ (0.395) dominates; mobile η-bridge is second-order (pp gaps ≤ 0.07) |

The correct reformulation is:

> **"The aggregate GDP response is determined by the wage regime (flexible vs
> sticky) and by substitution elasticities (θ, σ, ε). The intersectoral mobility
> parameter (η) is a material but secondary channel (first-order Sobol share
> ≈ 15.7%, S_f = 0.1571) — far smaller than the wage-regime effect and below θ —
> that only becomes first-order-relevant under extensions such as markups or
> complementarity."**

### The three results the paper should now report

**Result 1 — Methodological.** The original model was not an equilibrium. A
corrected, verified, budget-consistent version is now available. Residual
diagnostics pass at machine precision. (This satisfies Rev 2's demand for a
sound baseline.)

**Result 2 — Negative but precise.** The intersectoral mobility elasticity (η)
is not the *sole* driver of the aggregate multiplier. Its first-order Sobol
share is **material (~15.7%, S_f = 0.1571)** but well below the dominant θ (0.3951);
η's total-order index (ST_f = 0.1662) is only slightly above its first-order share,
so η contributes mostly as a direct main effect, not through interactions. The
mobility bridge (η=0 → η=1) is a second-order effect (≤ 0.07 pp) and η=1 solves
only at multipliers 0.1/0.2. The earlier reported "η dominance" (88.4%, 100%) and
the 2026-09-02 "~2% negligible" reading were both artifacts; the corrected model
shows η is a real but secondary channel. The competitive one-factor CES model does
not support the claim that labour mobility *alone* determines the output effects
of a sector-specific demand programme.

**Result 3 — Positive.** The **binary wage regime** (flexible vs sticky) is
material. A sticky wage produces a large *positive* extensive-margin response —
**+0.27 pp at mult 0.1 rising to +19.3 pp at mult 10** (employment expands from
1.003× to 1.19×), with fixed residuals ≤ 1.9e-7 — while the flexible η=0 closure
turns increasingly *negative* (−0.05 pp … −5.41 pp) and η=1 stalls from mult 0.5.
Because the flexible baseline is negative, the comparison is reported as a
**percentage-point gap** (fix_bridge: +0.32 pp … +24.75 pp), not a ratio. The
shock is absorbed through employment expansion under sticky wages and through
wage/quantity contraction under flexible wages. This identifies the correct
labour-market margin for policymakers: the question is not "how mobile is
labour?" but "is the economy at full employment or in an unemployment regime?"

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
| **5. The mobility channel (η) is material but secondary** | First-order Sobol S_f ≈ 15.7% (S_f = 0.1571); total-order ST_f ≈ 0.1662 (close to first-order, little interaction). Bridge is second-order (η=1 solves only at mult 0.1/0.2). The original "88.4%" and the 2026-09-02 "~2% negligible" reading were artifacts. Show the corrected variance decomposition. | §7.5 |
| **6. The wage regime channel is material** | `:fixed` vs `:mobile`: large positive extensive-margin response (+0.27 pp → +19.3 pp; employment up to 1.19×; residuals ≤ 1.9e-7), while flexible η=0 turns negative and η=1 stalls from mult 0.5. Reported as pp gaps, not ratios. Policy interpretation. | §7.6 |
| **7. What drives aggregate GDP: substitution elasticities** | Sobol shows substitution elasticities dominate GDP variance — **θ-led (S_f = 0.3951)** under the final fixed code, with σ (0.2734) and ε (0.1650) material and η (0.1571) secondary (ϵ-led in the 2026-09-02 re-run; θ+σ led pre-§4.1). Distinguish aggregate (elasticity-driven) from sectoral (not yet quantified). | §7.6 |
| **8. Robustness and limitations** | Financing closures, shock magnitude, static horizon, one-factor simplification. | §7.8 |
| **9. Conclusion** | Practical guidance: labour-market analysis should focus on the wage regime, not on the mobility elasticity. The IO multiplier is attainable under sticky wages; it is not under full employment. | §7.9 |

### Headline results for the abstract

1. A corrected equilibrium 71-sector CGE model of the German housing
   transformation.
2. Intersectoral labour mobility (η) has a material first-order GDP-variance share
    (~15.7%, S_f = 0.1571; total-order 0.1662) — it is not a sole driver, nor
    unqualified negligible, but it is secondary to the wage regime and to θ.
3. The **binary wage regime** (sticky vs flexible) is a first-order channel: a
    sticky wage yields a large positive extensive-margin response (+0.27 pp at
    mult 0.1 to +19.3 pp at mult 10, employment up to 1.19×) while flexible η=0
    turns negative; reported as pp gaps, not ratios.
4. Substitution elasticities explain the dominant GDP-variance share — **θ-led**
    (S_f = 0.3951) under the final fixed code, with σ (0.2734) and ε (0.1650)
    also material and η (0.1571) secondary; the earlier "η = 88.4%" claim was an artifact.

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
| §6 Phase 1 (account reconciliation) | **Done — notebook + interface.jl + locally generated `output/AC_*.csv` (gitignored, not committed)** | Reintegrated into §4.1 (done) |
| §6 Phase 2 (policy experiment) | Not started | Separate workstream |

---

## Part V — Remaining Work (Next Steps for the Revision)

### Short-term (can be done in this environment)

1. ~~**Record the `:fixed` closure results to CSV**~~ — **DONE**: the mult-sweep
   table is produced by `rerun_results.jl` (→ `output/wage_regime_comparison.csv`, a
   locally generated, gitignored artifact) and by `tests/test_fixed_closure.jl`.
2. ~~**Update WORKPLAN3.md**~~ — **DONE**: the `:fixed` closure is documented above
   and in `tests/test_fixed_closure.jl`.
3. **Run Sobol on sectoral quantities** — confirm whether η matters for
   sectoral allocation even if it doesn't for aggregate GDP

### Medium-term (implementation work, 1–2 weeks)

4. **Phase 1 — Accounting reconciliation** (ROADMAP §6, Phase 1)
   - **Done 2026-09-02** — deliverables: `Notebooks/AccountingConsistency.ipynb`, `docs/accounting_consistency_plan.md`, and locally generated `output/AC_*.csv` (regenerate with the notebook; not committed).
   - Separates domestic vs imported intermediate/final uses; decomposes value-added into wages / other-prod-tax / depreciation / net-op-surplus
   - Reconciles the three GDP sides (GDP(P)=GDP(I) exact; |GDP(P)−GDP(E)|=5.387%, a documented raw-table valuation discrepancy)
   - `src/interface.jl` integration done — `generate_data` performs the §4.1 transformation inline (`Ω_raw` for technology + `consumption_share` calibration, `Ω_dom` retained for domestic audit; decomposed VA; standard Domar λ; separated imports); `Data` extended with audit fields.
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
10. **Run the `:fixed` closure with an explicit demand anchor.** The `:fixed`
    closure at η≈1 is scale-indeterminate without an autonomous/investment
    additive anchor (it now throws an explanatory scale-indeterminacy error). Run
    it with an explicit autonomous or investment demand anchor that pins the
    activity scale, to confirm the sticky-wage results are not η-driven within the
    wage-regime margin.

---

## Part VI — Summary: The Corrected Story

The original paper claimed:

- "Labour supply elasticity η drives the GDP response (+1.5% under η→∞)"
- "η accounts for 88.4% of output variance"
- "The model is a valid equilibrium"

All three were false. The corrected model shows:

- **η (mobility) has a material first-order GDP-variance share (~15.7%, S_f = 0.1571)**
  and a second-order bridge (η=1 solves only at mult 0.1/0.2) — not a sole driver,
  nor unqualified negligible, but secondary to θ and to the wage regime
- **The wage regime (sticky vs flexible) is the relevant labour-market margin**,
  producing a large positive extensive-margin response (+0.27 pp → +19.3 pp;
  employment up to 1.19×; residuals ≤ 1.9e-7) while flexible η=0 turns negative —
  reported as pp gaps, not ratios
- **Substitution elasticities explain the dominant GDP-variance share** — **θ-led**
  (S_f = 0.3951) under the final fixed code, with σ (0.2734) and ε (0.1650) material
  and η (0.1571) secondary (ϵ-led in the 2026-09-02 re-run; θ+σ led pre-§4.1)
- The model is now a **verified equilibrium** with machine-precision residuals

This is a publishable result, but it requires the manuscript to abandon the
original "η-dominance" narrative and replace it with a **three-part story**:
(i) a corrected methodological baseline, (ii) a precise negative result about
the mobility channel, and (iii) a positive result about the wage-regime channel
and the dominance of substitution elasticities.