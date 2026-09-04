---
title: "Top-Level Workplan -- BFrep Manuscript Revision (Post-Pilot)"
project: Kapeller & Scharnreitner, submitted to Metroeconomica
date: July 2026
tags: [workplan, revision, io-modeling, cge, labor-supply, beyond-hulten, overview]
---

# Top-Level Workplan -- BFrep Manuscript Revision (Post-Pilot)

## Overview

This document coordinates three parallel work streams for salvaging the manuscript "A pluralist perspective on input-output modeling: searching for commensurability" (Kapeller & Scharnreitner, submitted to *Metroeconomica*).

The revision is based on the finding that the existing code base (`(3)BeyondHulten`) is approximately 70--80% functional -- far more advanced than earlier assessments assumed. This reduces the total timeline from 12--16 weeks to 6--10 weeks.

**Update (2026-07-29):** The four critical coding tasks from the Priority Matrix are **complete and verified** on Julia 1.12.6. The Go/No-Go decision is **GO**. The critical path has shifted from model implementation to manuscript writing and integration. Remaining timeline is estimated at 4--6 weeks.

> **Superseded (per `roadmaps/vertdict.md`).** The "complete and verified" status, the **GO** decision, the reported **price-invariance finding**, and the **88.4%** variance share from this round are **not supportable** and must not be treated as evidentially admissible. An independent two-sector diagnostic shows the supplied `mobile_labor.jl` / `variance_decomposition.jl` are not an equilibrium (household expenditure ≈ 25.5% exceeds income; a zero-profit equation was dropped instead of a redundant market equation; nominal rather than real wages were used; the high-η limit was misinterpreted as unlimited labor; and the 88.4% figure was a renormalized main-effect share, not a first-order sensitivity index). All mobile-labor sweep results and the GO decision derived from them are invalid until regenerated against the corrected equilibrium core in `ROADMAP.md`.

> **Historical note (2026-09-03).** The `tests/minimal_test/` harness and
> diagnostic scripts (`run_test.jl`, `test_full_pilot.jl`, `diag_*.jl`) referenced
> below have since been **deleted** and replaced by the organised test suite
> (`tests/test_*.jl`). The claims in this document are retained as a historical
> record of the rejected round (see the supersession note above), not as
> evidentially admissible results.

---

## Completed Milestones (This Round)

### ✅ Mobile Labor + Economy-Wide Wage — ⚠️ REJECTED (see `roadmaps/vertdict.md`)

Replaced sector-specific (immobile) wages with a single economy-wide wage `w`. Sectoral labor allocation is now endogenous via the marginal product condition. Unknown vector grew from 2N to 2N+1, with a numeraire constraint (CPI = 1) replacing one price equation to break the price-level indeterminacy inherent in CRTS models.

**Key structural insight (reported, NOT VALIDATED):** Under CRTS + mobile labor, the zero-profit condition determines `(p, w)` as a function of technology parameters alone. The labor supply elasticity η only affects the **scale** of the economy (total employment, GDP, sectoral quantities) -- **not** relative prices or the wage. Verified empirically: price CV = 0.0000 across the full η sweep, while quantity CV > 0.01 in 26/71 sectors. **This price invariance is an artifact of the rejected equilibrium system (dropped zero-profit equation, nominal rather than real wage, unnormalised demand shock), not a validated CGE result.** The "labor supply elasticity continuum" finding must be regenerated after the equilibrium and labor-supply corrections.

### ✅ η Parameter (Labor Supply Elasticity) — ⚠️ REJECTED SPECIFICATION

Implemented `L = L̄ · w^η` as the labor supply function. **Correction required:** labor supply must be a function of the *real* wage, `L = L̄[(w/P)/(w₀/P₀)]ᵝ`, not the nominal wage. The labor market clearing equation pins down the economy's scale.

- **η = 0**: perfectly inelastic (fixed total supply, Leontief-like)
- **η → ∞**: converges to the **fixed-real-wage closure** `w/P = ω̄` (not "unlimited labor"); with `w/P < 1`, `L = L̄(w/P)ᵝ → 0`, so the high-η limit must be implemented as an explicit complementarity/fixed-real-wage closure, not as `η = 10⁶`.
- **0 < η < ∞**: intermediate continuum

### ✅ Variance Decomposition — ⚠️ REJECTED (see `roadmaps/vertdict.md`)

Full factorial design over η × ε × θ × σ (5×3×3×3 = 135 evaluations) with ANOVA-style partial R²:

| Factor | Partial R² | Share |
|--------|-----------|-------|
| **η**  | **0.859** | **88.4%** |
| ε      | 0.075     | 7.7%  |
| θ      | 0.034     | 3.5%  |
| σ      | 0.004     | 0.4%  |
| Total  | 0.971     |       |
| Residual (interactions) | 0.029 | |

**Reported claim:** η dominates — it accounts for ~88% of explained variance in real GDP. **Audit: invalid.** The reported 88.4% is η's share of the *sum of reported main effects*; the percentages were renormalized over main effects rather than expressed as first-order sensitivity indices `S_f = Var(E[Y|f])/Var(Y)` under a documented product measure. The decomposition must be rebuilt per `ROADMAP.md` Phase 5 (explicit probability weights, first-order + total-effect indices, no main-effect renormalization, failures not silently discarded, complete designs required).

### ✅ Pilot η Sweep — **GO** 🟢 ⚠️ SUPERSEDED

| Criterion | Result | Threshold | Pass? |
|-----------|--------|-----------|-------|
| Sectors with CV > 0.01 | 26/71 | ≥5 | ✅* |
| Max sectoral CV | 0.057 | >0.01 | ✅* |
| η partial R² | 0.859 | >0.05 | ✅* |

\* Derived from the rejected equilibrium system; the GO decision resting on these numbers is **not supportable** until regenerated.

Top-5 most affected sectors (by quantity CV):
1. Grundstücks- u. Wohnungswesens (CV=0.057)
2. Kraftwagen und Kraftwagenteile (CV=0.049)
3. Vorb. Baustellen-, Bauinstallations-, Ausbauarbeiten (CV=0.040) — the shocked sector
4. Großhandelsleistungen (CV=0.034)
5. Öffentliche Verwaltung (CV=0.031)

### ✅ Julia Environment Aligned to Host (1.12.6)

- `Project.toml`: added `[compat] julia = "~1.12"`
- `Manifest.toml`: targets `julia_version = "1.12.6"`
- Julia 1.12.6 Linux ARM64 binary persists at `/root/.julia/linux-julia-1.12.6/` (mounted host filesystem)
- Minimal test project (`tests/minimal_test/`) runs all tests successfully
- Code is now portable between container (Linux) and host (macOS)

### Bug Fix

Fixed Julia destructuring bug in `mobile_labor.jl`: `(; supply_shock, factor_share) = (shocks, data)` was creating a plain tuple instead of accessing struct fields. Corrected to separate `(; supply_shock) = shocks` and `(; factor_share) = data`.

---

## Work Streams

### Stream A: Manuscript Revision

**Location:** `(1)Submission/revised/`

**Status:** Scaffold created. The manuscript has been reorganized into a modular structure:

```text
revised/
---- main.tex                      -- new main file with revised abstract
---- BandF.bib                    -- bibliography (copied)
---- pictures/                    -- figures (copied)
---- chapters/
    ---- ch01_introduction.tex    -- % OLD TEXT FROM HERE ON
    ---- ch02_io_leontief.tex     -- % OLD TEXT FROM HERE ON
    ---- ch03_cge_model.tex       -- % OLD TEXT FROM HERE ON
    ---- ch04_data.tex            -- % OLD TEXT FROM HERE ON
    ---- ch05_ces_functions.tex   -- % OLD TEXT FROM HERE ON
    ---- ch06_labor_slack.tex     -- % OLD TEXT FROM HERE ON
    ---- ch07_conclusion.tex      -- % OLD TEXT FROM HERE ON
    ---- ch08_backmatter.tex      -- % OLD TEXT FROM HERE ON
```

The new abstract replaces the "paradigmatic cleavage" framing with a "labor supply elasticity continuum" framing, centered on Robinson (2006) and the quantitative decomposition of closure choice vs. technology choice.

**Remaining work:**

- Rewrite introduction (drop Kuhnian preamble to one paragraph, open with policy question, state upfront that IO is a special case of CGE)
- Rewrite Section 2 (streamline IO exposition, fix sign error in eq. 19, remove "supply-side vs demand-side" distinction, be transparent about open economy)
- Rewrite Section 4 (drop CES novelty claim, drop Cobb-Douglas as "archetypical neoclassical")
- **Replace Section 5** (exogenous labor slack → endogenous labor supply with $\eta$ parameter) — **model code is ready, results verified**
- **Rewrite Section 6** (drop "deus ex machina" hand-wringing, present variance decomposition results: η = 88.4% of explained variance) — **results are ready**
- Integrate new literature: Robinson (2006), Rose (1995), Dervis et al. (1982), de Melo & Tarr (1992), McGregor et al. (1996), Shoven & Whalley (1984)
- Fix minor issues (p.5 sentence, p.16 paragraph, Figure 3 axis)

**Dependencies:** ~~Requires results from Stream B ($\eta$ sweep, variance decomposition) before Sections 5--6 can be finalized.~~ **RESOLVED**: Stream B results (η sweep, variance decomposition, Go/No-Go) are complete. Sections 5--6 can proceed immediately.

### Stream B: Model Extension (BeyondHulten code base)

**Location:** `(3)BeyondHulten/src/`

**Status:** Core code is 70--80% functional. Three model types (CES, Leontief, Cobb-Douglas) work. German data is connected. Labor slack is implemented as a function parameter. Elasticity gradient and plotting infrastructure exist.

**Completed this round (claimed — see supersession note above):**
- ✅ Implement continuous labor supply elasticity $\eta$ in [0, $\infty$) — done in `mobile_labor.jl` ⚠️ used nominal wage, must use real wage
- ✅ Switch from sector-specific to mobile labor with economy-wide wage — done in `mobile_labor.jl` ⚠️ dropped a zero-profit eqn; rejected per `vertdict.md`
- ✅ Implement variance decomposition ($\eta$ × $\varepsilon$ × $\theta$ × $\sigma$ factorial with partial R²) — done in `variance_decomposition.jl` ⚠️ renormalized main-effect share, not first-order index; rejected
- ✅ Go/no-go decision: pilot $\eta$ sweep — **GO** (η = 88.4% of variance, 26/71 sectors with CV > 0.01) ⚠️ SUPERSEDED; not supportable

**Remaining work:**

- Open the economy (imports, exports, domestic vs. imported intermediates) or document the closed-economy mapping transparently
- Fix calibration documentation (summary table of $\alpha$, $\beta$, $\lambda$, Domar weights, elasticities)
- Update notebooks to match current API (some notebooks use older function signatures)
- Generate final figures for the manuscript (η sweep plots, variance decomposition bar chart, sectoral variation heatmap)
- ~~Full $\eta$ sweep + variance decomposition~~ — **DONE**
- ~~Implement $\eta$ parameter~~ — **DONE**
- ~~Switch to mobile labor~~ — **DONE**

**Dependencies:** ~~Can proceed immediately on existing code.~~ **Critical path cleared.** The B&F replication (Stream C) provides validation but is not a prerequisite.

### Stream C: B&F Replication

**Location:** `(3)BeyondHulten/bf_replication/`

**Status:** Workplan created (`REPLICATION_WORKPLAN.md`). No code yet.

**Purpose:** A clean-room replication of Baqaee & Farhi (2019) on US data. Serves two functions:

- Validates the model implementation against published results
- Provides the computational infrastructure for the "replicate and contextualize" strategy from the salvage plan

**Remaining work:** See `bf_replication/REPLICATION_WORKPLAN.md` for the 6-phase plan (R1--R6).

**Dependencies:** Can reuse code patterns from `BeyondHulten/src/`. Independent of Streams A and B.

---

> **Cross-cutting caveat.** All "✅ DONE / GO / 88.4% / results ready" markers elsewhere in this document (Stream B "Completed this round", Priority Matrix, Execution Order, Risk Register, Reviewer Response Strategy) refer to the **rejected** round audited in `roadmaps/vertdict.md`. They are retained as a historical record of what was claimed, **not** as evidentially admissible results. No mobile-labor figure, GO decision, or sensitivity share from this round may be reported or written into the manuscript until regenerated against the corrected equilibrium core in `ROADMAP.md`.

## Cross-References

| This workplan | References |
|---------------|-----------|
| Full salvage strategy | `(4)plans/REVISED_SALVAGE_PLAN.md` |
| Earlier plans (superseded) | `(4)plans/SALVAGE_PLAN.md`, `salvage_plan2.md`, `salvage1_comment.md`, `BF_CONTRIBUTION_RECEPTION.md` |
| Replication details | `(3)BeyondHulten/bf_replication/REPLICATION_WORKPLAN.md` |
| Manuscript source | `(1)Submission/revised/main.tex` + `chapters/` |
| Original manuscript | `(1)Submission/foundations_of_pdf/submission(1).tex` |
| Reviews | `(2)Reviews/metro-rev1.docx`, `metro-rev2.docx` |
| Core code | `(3)BeyondHulten/src/` |
| Original workplan | `(3)BeyondHulten/WORKPLAN.md` |
| Model source files | `src/mobile_labor.jl`, `src/variance_decomposition.jl` |
| Test suite | `tests/minimal_test/` (run_test.jl, test_full_pilot.jl) |
| Julia modeling skill | Hermes skill `beyond-hulten-julia-modeling` |

---

## Execution Order (Revised)

### ~~Weeks 1--2: Foundation~~ — COMPLETED

- ~~**Stream A:** Draft new introduction, rewrite Section 2~~ (pending)
- ~~**Stream B:** Implement $\eta$ parameter, switch to mobile labor, run pilot sweep~~ — **DONE**
- **Stream C:** Phase R1 (data pipeline), start Phase R2 (core solver)

### ~~End of Week 2: Go/No-Go~~ — DECISION: GO ✅

- ~~Check pilot $\eta$ sweep results: are sectoral dynamics interesting?~~ — **YES** (26/71 sectors, CV up to 0.057)
- ~~Check variance decomposition: does $\eta$ dominate ($\varepsilon$, $\theta$, $\sigma$)?~~ — **YES** (η = 88.4% of variance)
- ~~If yes to either: proceed with full plan~~ — **PROCEEDING**
- ~~If no to both: pivot to pure replication strategy~~ — **NOT NEEDED**

### Weeks 1--2: Writing Sprint (Current Phase)

- **Stream A:** Draft new introduction, rewrite Section 2 (fix sign error, remove "supply vs demand", state model closure clearly)
- **Stream A:** Write Section 5 (endogenous labor supply with η) — results ready
- **Stream A:** Write Section 6 (variance decomposition results) — results ready
- **Stream B:** Generate final figures (η sweep, variance decomposition, sectoral variation)
- **Stream B:** Calibration documentation
- **Stream C:** Phase R1--R2 (data pipeline, core solver)

### Weeks 3--4: Integration

- **Stream A:** Rewrite Sections 4--5 (CES, labor supply), integrate literature
- **Stream B:** Open economy (or transparent mapping), update notebooks
- **Stream C:** Phase R2 (solver), Phase R3 (shock simulation)

### Weeks 5--6: Finalization

- **Stream A:** Write results section (aggregate, sectoral, variance decomposition), discussion, conclusion
- **Stream A:** Copy-edit, fix minor issues, finalize manuscript
- **Stream A:** Write response letter to both reviewers
- **Stream C:** Phase R4 (elasticity gradient), Phase R5 (reallocation) — if time permits

### Week 7--8: Submission

- Final review of all three deliverables
- Submit revised manuscript + response letter + replication code

---

## Priority Matrix (Revised)

| Task | Stream | Priority | Status | Effort |
|------|--------|----------|--------|--------|
| ~~Implement $\eta$ parameter~~ | B | ~~Critical~~ | ✅ DONE | ~~3--5 days~~ |
| ~~Mobile labor (economy-wide wage)~~ | B | ~~Critical~~ | ✅ DONE | ~~2--3 days~~ |
| ~~Variance decomposition~~ | B | ~~Critical~~ | ✅ DONE | ~~2--3 days~~ |
| ~~Pilot $\eta$ sweep~~ | B | ~~Critical~~ | ✅ GO | ~~1 day~~ |
| New introduction | A | High | Pending | 2 days |
| Rewrite Section 2 (fixes) | A | High | Pending | 2 days |
| Write Section 5 (labor supply with η) | A | High | Pending — **results ready** | 2--3 days |
| Write Section 6 (variance decomposition) | A | High | Pending — **results ready** | 2 days |
| Literature integration | A | High | Pending | 2--3 days |
| Generate final figures | B | High | Pending | 1--2 days |
| Calibration documentation | B | Medium | Pending | 1 day |
| Open economy | B | Medium | Can be transparent mapping | 3--5 days |
| Rewrite Sections 4 | A | High | Pending | 2--3 days |
| Results section | A | High | Needs figures | 2--3 days |
| Response letter | A | High | Needs final manuscript | 2--3 days |
| Replication R1--R2 | C | Medium | Validates model | 5--10 days |
| Replication R3--R6 | C | Low | Nice to have | 10--15 days |
| Update notebooks | B | Low | Maintenance | 1--2 days |

---

## Risk Register (Revised)

| Risk | Likelihood | Impact | Mitigation | Status |
|------|-----------|--------|-----------|--------|
| ~~$\eta$ sweep results are smooth/monotonic~~ | ~~Medium~~ | ~~High~~ | ~~Size claim is insurance~~ | **RESOLVED**: η dominates (88.4%), 26/71 sectors show variation |
| $\eta$ affects only scale, not prices | — | Medium | **Embrace as finding**: standard CGE result, strengthens "closure choice" framing | NEW — key insight for manuscript |
| Reviewer 2 irreconcilable | Medium | High (desk reject) | Maximally gracious response; size claim is genuine contribution | Unchanged |
| ~~Mobile labor implementation difficult~~ | ~~Low~~ | ~~Medium~~ | ~~Well-understood CGE modification~~ | **RESOLVED**: implemented and verified |
| Open economy too hard | Medium | Medium | Fallback: transparent documentation of closed-economy mapping | Unchanged |
| Replication takes too long | Medium | Low | Replication is not blocking; model extension is complete | Unchanged |
| Full project (GLMakie etc.) OOM in container | High | Low | Run on host Mac (more RAM); use minimal_test for container work | NEW |

---

## Reviewer Response Strategy

### Reviewer 2 (devastating, recommends rejection)

- Acknowledge core criticism fully: "the bridge already exists" (Robinson 2006)
- Show endogenous labor supply replaces exogenous injection — **now implemented and verified**
- Present variance decomposition as genuinely new quantitative finding — **η = 88.4% of explained variance**
- Address every specific point (sign error, closed economy, immobile labor, CES novelty, etc.)
- Be maximally gracious -- do not defend the original framing
- **New ammunition**: The CRTS result (η affects scale, not prices) is a clean theoretical finding that demonstrates understanding of CGE closure properties

### Reviewer 1 (moderate, constructive)

- Address each point systematically
- Show restructured paper with clear model description
- Include calibration summary table
- Engage with McGregor et al. (1996)
- Clarify aggregate vs sectoral distinction (now quantified via variance decomposition: η = 88.4% aggregate, 26/71 sectors show sectoral variation)
