---
title: "Top-Level Workplan -- BFrep Manuscript Revision"
project: Kapeller & Scharnreitner, submitted to Metroeconomica
date: July 2026
tags: [workplan, revision, io-modeling, cge, labor-supply, beyond-hulten, overview]
---

# Top-Level Workplan -- BFrep Manuscript Revision

## Overview

This document coordinates three parallel work streams for salvaging the manuscript "A pluralist perspective on input-output modeling: searching for commensurability" (Kapeller & Scharnreitner, submitted to *Metroeconomica*).

The revision is based on the finding that the existing code base (`(3)BeyondHulten`) is approximately 70--80% functional -- far more advanced than earlier assessments assumed. This reduces the total timeline from 12--16 weeks to 6--10 weeks.

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
- Replace Section 5 (exogenous labor slack -> endogenous labor supply with $\eta$ parameter)
- Rewrite Section 6 (drop "deus ex machina" hand-wringing, present variance decomposition results)
- Integrate new literature: Robinson (2006), Rose (1995), Dervis et al. (1982), de Melo & Tarr (1992), McGregor et al. (1996), Shoven & Whalley (1984)
- Fix minor issues (p.5 sentence, p.16 paragraph, Figure 3 axis)

**Dependencies:** Requires results from Stream B ($\eta$ sweep, variance decomposition) before Sections 5--6 can be finalized.

### Stream B: Model Extension (BeyondHulten code base)

**Location:** `(3)BeyondHulten/src/`

**Status:** Core code is 70--80% functional. Three model types (CES, Leontief, Cobb-Douglas) work. German data is connected. Labor slack is implemented as a function parameter. Elasticity gradient and plotting infrastructure exist.

**Remaining work:**

- Implement continuous labor supply elasticity $\eta$ in [0, $\infty$) as a new labor slack function
- Switch from sector-specific to mobile labor with economy-wide wage
- Open the economy (imports, exports, domestic vs. imported intermediates) or document the closed-economy mapping transparently
- Implement variance decomposition ($\eta$ x $\varepsilon$ x $\theta$ x $\sigma$ factorial with partial R^2 or Sobol indices)
- Fix calibration documentation (summary table of $\alpha$, $\beta$, $\lambda$, Domar weights, elasticities)
- Update notebooks to match current API (some notebooks use older function signatures)
- Go/no-go decision: run pilot $\eta$ sweep, check whether sectoral results show interesting variation

**Dependencies:** Can proceed immediately on existing code. The B&F replication (Stream C) provides validation but is not a prerequisite.

### Stream C: B&F Replication

**Location:** `(3)BeyondHulten/bf_replication/`

**Status:** Workplan created (`REPLICATION_WORKPLAN.md`). No code yet.

**Purpose:** A clean-room replication of Baqaee & Farhi (2019) on US data. Serves two functions:

- Validates the model implementation against published results
- Provides the computational infrastructure for the "replicate and contextualize" strategy from the salvage plan

**Remaining work:** See `bf_replication/REPLICATION_WORKPLAN.md` for the 6-phase plan (R1--R6).

**Dependencies:** Can reuse code patterns from `BeyondHulten/src/`. Independent of Streams A and B.

---

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

---

## Execution Order

### Weeks 1--2: Foundation

- **Stream A:** Draft new introduction, rewrite Section 2 (fix sign error, remove "supply vs demand", state model closure clearly)
- **Stream B:** Implement $\eta$ parameter, switch to mobile labor, run pilot sweep
- **Stream C:** Phase R1 (data pipeline), start Phase R2 (core solver)

### End of Week 2: Go/No-Go

- Check pilot $\eta$ sweep results: are sectoral dynamics interesting?
- Check variance decomposition: does $\eta$ dominate ($\varepsilon$, $\theta$, $\sigma$)?
- If yes to either: proceed with full plan
- If no to both: pivot to pure replication strategy (Stream C + Stream A)

### Weeks 3--4: Core Results

- **Stream A:** Rewrite Sections 4--5 (CES, labor slack -> endogenous labor supply), integrate literature
- **Stream B:** Full $\eta$ sweep + variance decomposition, open economy (or transparent mapping), calibration documentation
- **Stream C:** Phase R2 (solver), Phase R3 (shock simulation)

### Weeks 5--6: Integration

- **Stream A:** Write results section (aggregate, sectoral, variance decomposition), discussion, conclusion
- **Stream B:** Generate final figures, update plotting for new results
- **Stream C:** Phase R4 (elasticity gradient), Phase R5 (reallocation) -- if time permits

### Weeks 7--8: Finalization

- **Stream A:** Copy-edit, fix minor issues, finalize manuscript
- **Stream A:** Write response letter to both reviewers
- **Stream C:** Phase R6 (second-order effects) -- if time permits

### Week 9--10: Submission

- Final review of all three deliverables
- Submit revised manuscript + response letter + replication code

---

## Priority Matrix

| Task | Stream | Priority | Blocking? | Effort |
|------|--------|----------|-----------|--------|
| Implement $\eta$ parameter | B | Critical | Blocks A (Sec 5--6) | 3--5 days |
| Mobile labor (economy-wide wage) | B | Critical | Blocks $\eta$ implementation | 2--3 days |
| Variance decomposition | B | Critical | Blocks A (results section) | 2--3 days |
| Pilot $\eta$ sweep | B | Critical | Go/no-go decision | 1 day |
| New introduction | A | High | -- | 2 days |
| Rewrite Section 2 (fixes) | A | High | -- | 2 days |
| Literature integration | A | High | -- | 2--3 days |
| Open economy | B | Medium | Can be transparent mapping | 3--5 days |
| Calibration documentation | B | Medium | Blocks A (data section) | 1 day |
| Rewrite Sections 4--5 | A | High | Needs B results | 3--4 days |
| Results section | A | High | Needs B results | 2--3 days |
| Response letter | A | High | Needs final manuscript | 2--3 days |
| Replication R1--R2 | C | Medium | Validates model | 5--10 days |
| Replication R3--R6 | C | Low | Nice to have | 10--15 days |

---

## Risk Register

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|-----------|
| $\eta$ sweep results are smooth/monotonic | Medium | High (weakens novelty) | Size claim (variance decomposition) is insurance policy |
| Reviewer 2 irreconcilable | Medium | High (desk reject) | Maximally gracious response; size claim is genuine contribution |
| Mobile labor implementation difficult | Low | Medium | Well-understood CGE modification; existing code is modular |
| Open economy too hard | Medium | Medium | Fallback: transparent documentation of closed-economy mapping |
| Replication takes too long | Medium | Low | Replication is not blocking; model extension (Stream B) is the critical path |

---

## Reviewer Response Strategy

### Reviewer 2 (devastating, recommends rejection)

- Acknowledge core criticism fully: "the bridge already exists" (Robinson 2006)
- Show endogenous labor supply replaces exogenous injection
- Present variance decomposition as genuinely new quantitative finding
- Address every specific point (sign error, closed economy, immobile labor, CES novelty, etc.)
- Be maximally gracious -- do not defend the original framing

### Reviewer 1 (moderate, constructive)

- Address each point systematically
- Show restructured paper with clear model description
- Include calibration summary table
- Engage with McGregor et al. (1996)
- Clarify aggregate vs sectoral distinction (now quantified via variance decomposition)
