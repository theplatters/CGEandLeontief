## `WORKPLAN2.md` (Manuscript-Level, Three Streams)

### Stream A --- Manuscript Revision

| Task | Priority | Status |
|------|----------|--------|
| New introduction (drop Kuhnian preamble) | High | ❌ Pending |
| Rewrite Section 2 (sign error, remove supply vs demand) | High | ❌ Pending |
| Write Section 5 (endogenous labor supply with $\eta$) | High | ❌ Results ready, needs writing |
| Write Section 6 (variance decomposition) | High | ❌ Results ready, needs writing |
| Rewrite Section 4 (CES framing) | High | ❌ Pending |
| Literature integration (Robinson, Rose, Dervis, McGregor, Shoven/Whalley) | High | ❌ Pending |
| Results section | High | ❌ Needs figures |
| Response letter to reviewers | High | ❌ Needs final manuscript |

### Stream B --- Model Extension (BeyondHulten Core)

| Task | Priority | Status |
|------|----------|--------|
| $\eta$ parameter (continuous labor supply elasticity) | Critical | ✅ Complete |
| Mobile labor (economy-wide wage) | Critical | ✅ Complete |
| Variance decomposition ($\eta \times \epsilon \times \theta \times \sigma$) | Critical | ✅ Complete ($\eta$ = 88.4%) |
| Pilot $\eta$ sweep --- Go/No-Go | Critical | ✅ **GO** |
| Generate final figures ($\eta$ sweep, variance decomposition bar chart) | High | ❌ Pending |
| Calibration documentation (summary table) | Medium | ❌ Pending |
| Open economy (imports/exports) or transparent mapping | Medium | ❌ Pending |
| Update notebooks to match current API | Low | ❌ Pending |

### Stream C --- B&F Replication (bf_replication, 2019 Paper)

| Phase | Description | Status |
|-------|-------------|--------|
| R1 | Data pipeline (BFdata.csv loading) | ✅ Complete |
| R2 | Core solver (Newton-Raphson, verified against MATLAB) | ✅ Complete |
| R3 | Shock simulation (Monte Carlo, 50k draws from stfp.csv) | ❌ Pending |
| R4 | Elasticity gradient (sweep $\epsilon$, $\theta$, $\sigma$) | ❌ Pending |
| R5 | Mobile labor variant (reallocation) | ❌ Already implemented in Stream B |
| R6 | Second-order effects | ❌ Pending |
