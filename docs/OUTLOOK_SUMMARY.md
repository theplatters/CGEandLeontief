---
title: "Outlook Summary -- BFrep Manuscript Revision"
project: Kapeller & Scharnreitner, submitted to Metroeconomica
date: August 2026
tags: [outlook, revision-strategy, response-to-reviewers, priority-matrix, code-status]
---

# Response Strategy: From "Searching for Bridges" to "Measuring the Bridge"

The original manuscript framed the exercise as a search for commensurability between IO and CGE -- a "paradigmatic cleavage." Reviewer 2 called this superfluous: the bridge already exists (Robinson 2006), and a CGE model with fixed prices and unconstrained factor supply simply is an IO multiplier model. Reviewer 1 agreed, more gently: the paper was hard to follow and the central message was unclear.

Both reviewers were right on the merits. The response strategy is not to defend the old framing but to retire it and replace it with something the reviewers did not anticipate: actual measurement.

## The Pivot

The new abstract retires the entire "searching for commensurability" language. The paper now opens with the fact that the bridge is known (Robinson, Rose, Dervis et al.) and immediately asks the follow-up question no one has answered:

> How much does each side of the bridge actually matter?

The labor supply elasticity $\eta \in [0, \infty)$ is the parameter that governs the bridge. At $\eta = 0$ you get full-employment CGE; at $\eta \to \infty$ you get the IO multiplier. The paper's contribution is not discovering this -- it is measuring it across a 71-sector German IO table and quantifying that $\eta$ explains 88.4% of the variation in real GDP, an order of magnitude more than all production technology parameters ($\epsilon, \theta, \sigma$) combined.

This turns the reviewers' strongest criticism -- "the bridge already exists" -- into a strength. Yes, it exists. And now we know how much it matters.

## The New Abstract

It is well understood that the IO multiplier model is a special case of a CGE model with fixed prices and unconstrained factor supply [@Robinson.2006]. The labor supply elasticity~$\eta$ --- running from zero (full-employment CGE) to infinity (IO multiplier) --- is thus the key parameter governing the relationship between these modeling traditions. While this conceptual bridge has been recognized for decades [@Rose.1995; @Dervis.1982], its quantitative implications have not been systematically measured.

We replicate and extend the Baqaee--Farhi framework to quantify the relative importance of factor market closures versus production technology parameters. Using the German housing transformation toward carbon neutrality as a policy-relevant application, we conduct a systematic sensitivity analysis across the labor supply elasticity continuum and a range of CES substitution elasticities. We find that the factor market closure explains roughly an order of magnitude more variation in both aggregate and sectoral outcomes than all production technology parameters combined. The estimated GDP effect of the housing transformation ranges from virtually zero to approximately 1.5% depending on the closure, while CES elasticities shift this estimate by only a fraction of a percentage point. Moving from exogenous labor slack to endogenous labor supply produces the "bridge" the literature already knows exists --- and here we show how much it matters.

**Keywords:** input-output modeling, general equilibrium, labor supply elasticity, factor market closures, CES production functions, socio-ecological transformation

**JEL-Codes:** C67, C68, D57, E24

---

# How the Code Answers Both Reviewers

**Reviewer 2's core demand:** endogenize labor supply instead of exogenous ad-hoc injections. Done. The `mobile_labor.jl` module implements $L = \bar L \cdot w^\eta$ with an economy-wide wage, exactly as Robinson (2006) describes. No more sector-specific immobile labor, no more arbitrary allocation rules.

**Reviewer 2's second demand:** stop claiming CES is novel. Done -- the new abstract and revised manuscript frame CES as standard (Shoven/Whalley 1984, Mansur/Whalley 1984). The Cobb-Douglas = "archetypical neoclassical" language is dropped.

**Reviewer 1's demand:** clarify model structure, show calibration table, consider intersectoral mobility. Done -- mobile labor replaces immobile, the variance decomposition quantifies the aggregate-vs-sectoral distinction that Reviewer 1 flagged.

**Reviewer 1's key reference:** McGregor et al. (1996) -- long-run IO-type results under involuntary unemployment. The $\eta$ continuum subsumes this: $\eta \to \infty$ gives the McGregor et al. result, $\eta = 0$ gives the standard full-employment CGE.

---

# Priority Matrix (from WORKPLAN2.md)

## Completed -- Critical Path Cleared

| Task | Stream | Effort | Verdict |
|------|--------|--------|---------|
| $\eta$ parameter (continuous labor supply elasticity) | B | 3-5 days | Done |
| Mobile labor (economy-wide wage) | B | 2-3 days | Done |
| Variance decomposition ($\eta \times \epsilon \times \theta \times \sigma$) | B | 2-3 days | Done |
| Pilot $\eta$ sweep -- Go/No-Go decision | B | 1 day | GO (88.4% $\eta$, 26/71 sectors) |
| B&F replication data pipeline (76 sectors, BFdata.csv) | C | 3-5 days | Done |
| B&F core equilibrium solver (Newton-Raphson) | C | 3-5 days | Done |
| B&F inflation dynamics analysis | C | 2 days | Done |
| B&F MATLAB comparison (welfare GDP, oil shock) | C | 2 days | Done |

## High Priority -- Pending

| Task | Stream | Effort | Notes |
|------|--------|--------|-------|
| New introduction (drop Kuhnian preamble) | A | 2 days | Pending |
| Rewrite Section 2 (sign error fix, remove supply-vs-demand) | A | 2 days | Pending |
| Write Section 5 (endogenous labor supply with $\eta$) | A | 2-3 days | Results ready |
| Write Section 6 (variance decomposition results) | A | 2 days | Results ready |
| Literature integration (Robinson, Rose, Dervis, McGregor, Shoven/Whalley) | A | 2-3 days | Pending |
| Generate final figures ($\eta$ sweep, variance decomposition bar chart) | B | 1-2 days | Needs plotting code |
| Response letter to reviewers | A | 2-3 days | Needs final manuscript |

## Medium Priority

| Task | Stream | Effort | Notes |
|------|--------|--------|-------|
| Calibration documentation (summary table of $\alpha, \beta, \lambda$, Domar weights) | B | 1 day | Data ready |
| Open economy (imports, exports) or transparent mapping documentation | B | 3-5 days | Fallback: document closed-economy mapping |
| B&F R3 Monte Carlo simulation (50k draws from stfp.csv) | C | 5-7 days | Needs solver robustness fix |

## Low Priority -- Nice to Have

| Task | Stream | Effort | Notes |
|------|--------|--------|-------|
| B&F R4 Elasticity gradient (sweep $\epsilon, \theta, \sigma$) | C | 2-3 days | Not blocking manuscript |
| B&F R5 Mobile labor variant (reallocation) | C | 2-3 days | Already implemented in Stream B |
| B&F R6 Second-order effects | C | 4-6 days | Not blocking manuscript |
| Update notebooks to match current API | B | 1-2 days | Maintenance |

---

# Code Additions Completed (Detailed)

## BeyondHulten Core (`src/`)

**`mobile_labor.jl`** (239 lines)
- New model type `MobileLaborCES` with `MobileLaborCESElasticities($\theta, \epsilon, \sigma, \eta$)`
- Economy-wide wage w as scalar unknown (2N+1 system where N = 71 sectors)
- Sectoral labor demand $L_i = [p_i \cdot A_i^{(\epsilon-1)/\epsilon} \cdot \alpha_i^{1/\epsilon} \cdot y_i^{1/\epsilon} / w]^\epsilon$
- Labor market clearing: $\Sigma L_i = \bar L \cdot w^\eta$
- Numeraire: CPI = 1 (replaces one price equation to break CRTS price-level indeterminacy)
- Verified: price CV = 0.0000 across $\eta$ sweep, quantity CV > 0.01 in 26/71 sectors

**`variance_decomposition.jl`** (419 lines)
- `eta_sweep()` -- run model across $\eta$ grid with warm-start initialization
- `variance_decomposition()` -- full factorial ($5 \times 3 \times 3 \times 3 = 135$ evaluations)
- ANOVA-style $R^2$ for each elasticity factor
- Results: $\eta = 88.4\%$, $\epsilon = 7.7\%$, $\theta = 3.5\%$, $\sigma = 0.4\%$
- `pilot_eta_sweep()` -- Go/No-Go decision function with automatic diagnostics

## B&F Replication (`bf_replication/`)

**`src/data_loader.jl`** (210 lines)
- Reads `BFdata.csv` (Jorgenson 88-sector, 1960-2005)
- Removes sectors 60, 80-88 (government), 8 (uranium ores), 62 (machinery rental)
- Computes $\Omega, \alpha, \beta, \lambda$, L for base year 1980
- All consistency checks pass: Domar weights, factor shares, consumption shares

**`src/model.jl`** (346 lines)
- 2N equilibrium system: N price equations + N market clearing equations
- Custom Newton-Raphson solver with numerical Jacobian
- Converges to machine precision (~$10^{-15}$) for baseline and moderate TFP shocks
- Verified line-by-line against MATLAB `Simulation.m`

**`src/inflation_analysis.jl`** (423 lines)
- Price vs quantity decomposition under TFP shocks
- Network price propagation analysis (direct vs upstream vs downstream vs GE)
- Eta-insensitivity confirmation (CV = 0.000)

**`run_matlab_comparison.jl`** (220 lines)
- Side-by-side comparison: Julia output vs MATLAB specification
- Oil shock (sector 7, A = 0.7): welfare GDP $\Delta \log$ = -0.0355, Hulten = -0.0277, ratio = 1.28x
- Cross-sector comparison: Oil (1.28x amplified), Retail (0.87x dampened), Construction (0.88x dampened)
- Matches paper's Figure S2: oil-is-special finding

### REPLICATION_REPORT.md (224 lines)
- Full replication status report documenting what is verified vs. not yet verified against published paper results
- Qualitative confirmations: asymmetric nonlinearity, sign-dependent ranking, invariant CPI
- Limitation: Monte Carlo (50k draws) not yet run; $\theta = 1$ solver failure; no MATLAB runtime in container

---

# Timeline (Revised Post-Pilot)

| Phase | Duration | Focus |
|-------|----------|-------|
| Writing sprint | Weeks 1-2 | New intro, Sections 5-6, figures, literature |
| Integration | Weeks 3-4 | Full manuscript draft, calibration table, response letter |
| Finalization | Weeks 5-6 | Copy-edit, response letter, final figures |
| Submission | Week 7-8 | All three deliverables: manuscript + response + replication code |

The critical path has shifted from model implementation to manuscript writing. Four of the four critical coding milestones are complete and verified. The B&F replication (Stream C) is no longer on the critical path -- it validates but does not block the manuscript revision.