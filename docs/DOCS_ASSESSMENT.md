---
title: "Document Assessment: What Remains Relevant in docs/"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-09-14"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [assessment, docs-audit, revision, beyondhulten, metroeconomica]
last-updated: September 2026
---

**Version 7** \textcolor{revisionV6}{(September 2026)}
**Version 6** \textcolor{revisionV5}{(September 2026)}
**Version 5** \textcolor{revisionV4}{(September 2026)}
**Version 4** \textcolor{revisionV3}{(September 2026)}
**Version 3** \textcolor{revisionV2}{(September 2026)}
**Version 2** \textcolor{revisionV1}{(September 2026)}
**Version 1** (September 2026)

\textcolor{revisionV2}{This document is the working plan for the
Metroeconomica revision: Foundation I (labour-market closures),
Foundation II (financial closures), the 5 x 3 evaluation matrix, and the
staged workplan with integrated reviewer-response maps. It evolved from
the document assessment of Versions 1--2, whose chronological review now
lives in \texttt{docs/archive/document\_review\_leftovers.md}.
Governing documents: \texttt{ROADMAP.md} and
\texttt{roadmaps/vertdict.md}; active references:
\texttt{definitive\_guide.md}, \texttt{labor\_closures.md}, and
\texttt{docs/reviews/}.}

\textcolor{revisionV4}{Version 5 reports the executed 5 x 3 matrix on the re-anchored A-bill calibration, the accounting correction that closed the external-identity canary, and the two closure fixes (real-wage supply, verified scale determinacy) that made the fixed-wage rows executable. Section 4.2 is new; where the earlier intermediate-results section disagrees with it, section 4.2 governs.}

\textcolor{revisionV5}{Version 6 separates real GDP from welfare in the measurement layer (ADR-0018). Real GDP is now the income-side index w*sum(L) deflated by a seven-component Tornqvist GDP deflator, with the expenditure-side Divisia index as its dual and the external wedge reported as a diagnostic; the legacy household-consumption index is relabelled as welfare and is no longer reported as GDP. Because the kernel measurement layer changed, the matrix is re-minted as matrix\_5x3\_v4 (the same 15 cells and parameters as v3; v1-v3 stay as history, ADR-0004). Under the demand-only design the deflator is exactly one, so real GDP is flat (0.000 percent) in every full-employment mobile row: the F2 column shows the tax-financed programme as GDP-neutral with a 1.69 percent welfare loss, not as a GDP contraction. The v5 reading of the F2 column is corrected in section 4.2; the manifest-generated v4 flow tables are in \texttt{paper/tables/matrix\_5x3\_v4\_flows.md}.}

\textcolor{revisionV6}{Version 7 harmonizes the external account (ADR-0019). All N goods-market clearings are now enforced in every regime: the mobile eta = 1 system carries the endogenous net external transfer F with E = (1-tau) w sum(L) + F, the eta = 0 endpoint pins F = 0, and programme financing under F3 is booked as B\_gov = sum(p g), so the booked identity S + T\_int + M - (I+X) = F + B\_gov holds to about 1e-16 at the mobile cells. Because the kernel changed, the matrix is re-minted as matrix\_5x3\_v5 (15 of 15 cells executed; v1-v4 stay as history, ADR-0004); the manifest-generated flow tables are in \texttt{paper/tables/matrix\_5x3\_v5\_flows.md}. Financing neutrality is the generation's finding: F2 and F3 have identical real allocations in the mobile regime with F\_F3 = F\_F2 - B\_gov and an identical net external position, so the F2/F3 differences reported in v1-v4 were artifacts of the omitted-market shortcut, not economics. Measured movement against v4: mobile F1 welfare +0.043 percent becomes +0.143 percent, mobile F2/F3 welfare -1.694 percent becomes -1.823 percent, and the fixed rows and the BF rows are unchanged. Correction of the same date, within this version: the BF eta = 0 external-position entries are booked financing rather than positions -- the F = 0 pin reduces them to the programme cost B\_gov = +1.331 percent of GDP, BF-F3's account is open by -0.7936 percent of GDP, and the eta = 0 position is not identified (every pin F = c is an exact root); see section 4.2 and ADR-0020 (accepted, option C, and promoted the same day: the matrix is re-minted as matrix\_5x3\_v6, whose BF row now has an identified external position and a closed account while the twelve non-BF cells reproduce v5).}

\textcolor{revisionV3}{Version 4 reports the intermediate implementation
results of the notebook pipeline \texttt{cbase2/} (v3 open-economy
recalibration): the Stage 0 chapter is promoted to its own section with an
updated status matrix, the remaining stages form the new Workplan section,
and Stage 1 carries its updated state.}

The document now opens with the two foundation sections whose decisions
gate everything else -- the labour-market closure (Foundation I) and the
financing closure (Foundation II) -- followed by the chronological
workplan. \textcolor{revisionV1}{All of this is part of the Version 2
reworking; the option tables in the foundation sections carry, per row,
formulation, status, referee point, cost, and recommended role, plus the
provenance tags (reviewer, old docs, Baqaee--Farhi usage).}

# Foundation I: Labour-Market Closure \textcolor{revisionV1}{\normalsize [added v2]}

\textcolor{revisionV1}{The referee demand (R2) is an explicit
labour-supply elasticity; the code currently implements a geometric
reallocation parameter and two wage regimes. The option tables
below consolidate the full design space, grouped into the selected
closures (BF plus the four standard alternatives ALPHA--DELTA), notable
unselected candidates, and the dropped exogenous-shift closure. Provenance tags: "reviewer"
attributes R1/R2; "old docs" names documents (archived under
\texttt{docs/archive/}, plus the governing \texttt{ROADMAP.md} and
\texttt{roadmaps/vertdict.md}); "used by BF" points to the Baqaee--Farhi
papers (both source-verified in \texttt{bf\_replication/} and
\texttt{bf\_replication2/}):}

- BF 2019 (Econometrica): `Dokumente/Baqaee-2019-The Macroeconomic Impact of Microeconomic Shocks- Beyond Hulten's Theorem-Econometrica.pdf`
- BF 2022 (AER): `Dokumente/Baqaee-2022-Supply and Demand in Disaggregated Keynesian Economies with an Application to the COVID-19 Crisis-American Economic Review_4.pdf`

## Selected closures \textcolor{revisionV1}{\normalsize [added v2]}

\textcolor{revisionV1}{Suggested narrative: the paper measures how far the estimated effects of a
sector-specific green-investment programme travel between the CGE and IO
endpoints as the labour-market closure varies. Closure BF operationalizes
the neoclassical short-run argument as a market friction: BF 2019's own
framework contrasts full reallocation with none, and the $\eta$ parameter
makes that friction continuous -- limits on intersectoral reallocation,
not a change of paradigm, separate the short run from the benchmark. On
the shared production-network core, four standard closures span the
space: ALPHA fixes aggregate labour at full employment (the neoclassical
standard); BETA relaxes the endowment through a labour--leisure
elasticity (the voluntary response the referees ask for); GAMMA fixes the
real wage and lets employment absorb the shock (production-agnostic
Keynesian: the mechanism operates on the labour market rather than
through the production function, so the extensive-margin result survives
any production structure); DELTA fixes the real factor price with
unconstrained quantity and fixed coefficients, emulating the Leontief
inverse (Robinson 2006). The measurement claim: the wage regime
(ALPHA/BETA vs.\ GAMMA) and the substitution elasticities govern the
distance to DELTA; the reallocation friction BF is material but
secondary; and the financing closure keeps every experiment admissible.}

| **Closure BF: Geometric reallocation parameter ($\eta$, market friction)** | *Endogenous sectoral allocation between immobile and cost-minimizing; the BF-comparable short-run friction* |
|---|---|
| **Formulation** | $L_i = L_{i,\mathrm{fixed}}^{1-\eta}\, L_{i,\mathrm{costmin}}^{\eta}$ with $\eta \in [0,1]$ (extrapolation outside); total employment fixed at $\bar L$; **not** a labour-supply elasticity |
| **Status** | Implemented and tested (`MobileLaborCES`, `:mobile`); Sobol first-order share about 15.7% |
| **Referee point addressed** | The mobility margin only; total supply stays fixed |
| **Cost** | None |
| **Recommended role** | **Selected**: the BF-comparable market-friction closure; the interpolation parameter is reported under an honest name ("reallocation parameter") |
| **Mentioned by reviewer** | both (R1.6 mobility/wage-rule question; R2.5 immobile-labour criticism -- indirectly) |
| **Mentioned in old docs** | `labor_closures.md`; `definitive_guide.md`; `WORKPLAN3.md` (2026-09-03 addendum) |
| **Used by BF** | Yes at its endpoints -- the original BF 2019 replication ships both variants: **immobile labour** ("no reallocation"; sector-specific wage vector, `Simulation_Derivs.m`) and **mobile labour** ("full reallocation"; single wage pinned at the numeraire, `Simulation_Derivs_realloc.m`; `Read_Me.txt`: "the first replicates results with immobile labor and the second replicates results with mobile labor"). The geometric interpolation parameter *between* the endpoints is project-specific |

| **Closure ALPHA: Full-employment mobile labour (neoclassical standard)** | *The standard full-employment CGE endpoint* |
|---|---|
| **Formulation** | $\sum_i L_i = \bar L$; one economy-wide wage; cost-minimizing allocation; the $\eta = 1$ limit of the reallocation parameter (as `:mobile`) |
| **Status** | Implemented (`:mobile`; ROADMAP Closure A); baseline endpoint of the wage-regime sweep |
| **Referee point addressed** | R2.5 (economy-wide mobile labour); R1.6 (unified wage rule enabling mobility) |
| **Cost** | None |
| **Recommended role** | Baseline endpoint of the closure continuum |
| **Mentioned by reviewer** | both (R2.5 explicitly recommends it; R1.6 asks for it) |
| **Mentioned in old docs** | `REVISED_SALVAGE_PLAN.md`; `WORKPLAN.md`; `WORKPLAN2.md`; `vertdict.md`; `definitive_guide.md` |
| **Used by BF** | Yes -- BF 2019 baseline: competitive core with inelastically supplied (fixed) composite-factor endowment and uniform factor prices (see filename above) |

| **Closure BETA: Elastic total labour supply (neoclassical, labour--leisure)** | *Voluntary labour-supply response along the real wage; the continuum from vertical to horizontal supply* |
|---|---|
| **Formulation** | $L^s = \bar L\, [(w/P)/(w_0/P_0)]^{\eta_s}$ with a labour--leisure interpretation; $\eta_s = 0$ fixed supply; $\eta_s \to \infty$ approaches the fixed-real-wage limit |
| **Status** | Not implemented; specified in ROADMAP sections 4.2 and 5 (Closure B) |
| **Referee point addressed** | R2's core demand ("the one obvious set of elasticities that really matters") |
| **Cost** | Moderate: extends the $2N{+}1$ system by the supply equation and its anchor |
| **Recommended role** | **Selected** as the referee-facing principal addition |
| **Mentioned by reviewer** | both (R2 explicitly and centrally; R1 implicitly via the slack-scenario discussion) |
| **Mentioned in old docs** | `REVISED_SALVAGE_PLAN.md` (original `elastic_labor_slack` proposal); `WORKPLAN.md`; `WORKPLAN2.md`; `PRELIMINARY-ASSESSMENT.md` (baseline-anchoring correction); `vertdict.md` (ElasticLaborSupply); `ROADMAP.md` |
| **Used by BF** | No (2019: inelastic endowments; 2022: sticky wages -- vertdict names the 2019 base as the right foundation for adding it) |

| **Closure GAMMA: Sticky real wage, uncapped (production-agnostic Keynesian)** | *Fixed real wage; employment absorbs the shock -- the extensive margin* |
|---|---|
| **Formulation** | $w/P = \bar\omega$ (wage pinned at baseline); employment endogenous and uncapped; the gap $\bar L - \sum_i L_i$ is computed post-solve, not an equilibrium constraint |
| **Status** | Implemented; large positive GDP responses (+19.3 pp at mult 10, employment to 1.19x), residuals at most $1.9 \times 10^{-7}$ |
| **Referee point addressed** | R1.5 (wage-setting rule: fixed real/nominal wage); R2's "fix the wage rate at baseline and drop the labour supply constraint" recipe |
| **Cost** | None |
| **Recommended role** | Report as-is, with the uncapped-employment caveat stated |
| **Mentioned by reviewer** | both (R2 gives the recipe explicitly; R1.5 demands the wage-setting rule) |
| **Mentioned in old docs** | `WORKPLAN3.md` (Milestone F + 2026-09-03 addendum); `definitive_guide.md` (Result 3); `labor_closures.md` |
| **Used by BF** | Yes (adapted) -- BF 2022's sticky-wage mechanism (downward wage rigidity, demand-determined employment) in simplified one-wage form (see filename above) |

| **Closure DELTA: IO endpoint (emulating the Leontief inverse)** | *Fixed real factor price, unconstrained factor quantity; the exact conditions under which the CGE reproduces the IO multiplier* |
|---|---|
| **Formulation** | Fixed real primary-factor price; unconstrained factor quantity within the experiment; Leontief input coefficients (or a demonstrated CES special case with unchanged relative prices) |
| **Status** | Not implemented as an explicit closure (ROADMAP Closure D); GAMMA approximates its spirit |
| **Referee point addressed** | R2 explicitly (Robinson 2006: fix primary input prices, unconstrained supplies -- the multiplier model); R1 explicitly (McGregor--Swales--Yin long-run zero-price-change IO-type result) |
| **Cost** | Small--moderate |
| **Recommended role** | Required: R2's entire critique turns on the exact conditions of equivalence |
| **Mentioned by reviewer** | both (R2 with the Robinson 2006 quote; R1 with McGregor et al. 1996) |
| **Mentioned in old docs** | `vertdict.md` (FixedRealWage endpoint); `ROADMAP.md` (Closure D); `definitive_guide.md` (Part IV); `REVISED_SALVAGE_PLAN.md` (the superseded $\eta \to \infty$ endpoint idea) |
| **Used by BF** | No (origin is the IO tradition and Robinson 2006, not the BF papers; BF 2019 computes efficient reallocations, not fixed-price multipliers) |

## Notable candidates, but unselected \textcolor{revisionV1}{\normalsize [added v2]}

\textcolor{revisionV1}{Structural alternatives that did not make the
selection: the partial-mobility 3N system (allocation-axis alternative)
and the unemployment complementarity (optional, policy-evaluation only).}

| **Standard partial-mobility formulation (3N system)** | *Structural alternative: sectoral wages with wage-responsive sectoral supply* |
|---|---|
| **Formulation** | $L_i^s = \bar L\, s_i^0 (w_i/\bar w_i)^{\nu} \,/\, \sum_j s_j^0 (w_j/\bar w_j)^{\nu}$; sectoral labour demand from Shephard's lemma; $N$ zero-profit + $(N{-}1)$ market-clearing + $N$ sectoral labour-market + one numeraire = $3N$ equations |
| **Status** | Open design proposal in `labor_closures.md`; not implemented |
| **Referee point addressed** | None directly; referee-adjacent (R1.6, R2.5) |
| **Cost** | Largest: a new equilibrium system |
| **Recommended role** | Fallback if the one-wage supply formulation of Closure B is judged too reduced-form |
| **Mentioned by reviewer** | no |
| **Mentioned in old docs** | `labor_closures.md` (the proposal, with references: Hsieh--Klenow 2009; Baqaee--Farhi 2020 QJE; Artuc--Chaudhuri--McLaren 2010) |
| **Used by BF** | No (closest analogue: BF 2022's $N$ sector-specific labour markets, but those are sticky, not supply-responsive) |

| **Closure C: Unemployment complementarity (capped)** | *Employment capped at $\bar L$; full-employment or wage-floor regime switches endogenously* |
|---|---|
| **Formulation** | $0 \le \bar L - L \perp w/P - \bar\omega \ge 0$: complementarity between the employment gap and the real-wage gap |
| **Status** | Not implemented (ROADMAP section 5, Closure C; definitive guide Part IV records "partially") |
| **Referee point addressed** | R1.5 (wage-setting rule, wage curve, involuntary unemployment); R2 (DDR/Robinson unemployment closures; skill-class split) |
| **Cost** | Larger: complementarity solver, new results |
| **Recommended role** | Optional; mandatory only if the paper claims policy evaluation (ROADMAP Phase 6) |
| **Mentioned by reviewer** | both (R1.5 wage curve; R2's closure-literature demands) |
| **Mentioned in old docs** | `vertdict.md` (UnemploymentComplementarity); `ROADMAP.md`; `definitive_guide.md` (Part IV status) |
| **Used by BF** | Yes (adapted) -- BF 2022's complementarity machinery for binding wage floors / labour capacity, to be stripped down; its $N$-sector sticky-labour structure must not be copied unchanged (verdict; see filename above) |


## Dropped \textcolor{revisionV1}{\normalsize [added v2]}

\textcolor{revisionV1}{Regretfully retired: the exogenous labour-endowment
shift in its legacy-callback implementation -- the original paper's own
mechanism. The lineage note below records its direct resonance with
BF 2022 and with the response letter's methodological story.}

| **Dropped: Exogenous labour-endowment shift (the legacy labour-slack callback)** | *The original paper's own mechanism: exogenous sectoral slack vectors -- the same closure BF 2022 calibrates as labour-supply shifts* |
|---|---|
| **Formulation** | Sectoral labour endowments shifted exogenously: $\bar L_i \to \bar L_i (1+g_i)$ with calibrated or assumed $g_i$; no wage response required. Code implementations: `full_labor_slack` / `full_labor_slack_alt` (Leontief-projection sectoral vectors, manuscript Section 5.1) and `empirical_labor_slack` (uniform $\mathbf{l}_{\mathrm{new}} = (1-\mu)^{-1}\,\mathbf{l}$, Section 5.2) |
| **Status** | Implemented (legacy); the original submission's Sections 5.1--5.2; ROADMAP section 4.4 lists it as a distinct closure concept |
| **Referee point addressed** | R1.5 describes precisely this mechanism ("is it simply a matter of increasing labour supply by the number of unemployed persons?"); R2 criticizes its use as the central bridge (deus ex machina; R2.10 uniform allocation) |
| **Cost** | None |
| **Recommended role** | **Dropped** from the model; its economics survive in the lineage note below and in the response letter's methodological story |
| **Mentioned by reviewer** | both (R1 describes the mechanism; R2 criticizes its use as deus ex machina) |
| **Mentioned in old docs** | `REVISED_SALVAGE_PLAN.md`; `WORKPLAN.md`; `WORKPLAN2.md`; `labor_closures.md`; `definitive_guide.md`; `WORKPLAN3.md` |
| **Used by BF** | Yes -- BF 2022 calibrates exogenous sectoral labour-supply shifts to BLS May-2020 sectoral hours (contractions; `bf_replication2/src/model.jl` applies them on the labour block); run in reverse, the same closure is the slack injection. BF 2019 shocks factor endowments in its general framework |

\textcolor{revisionV1}{Lineage note (original manuscript). The resonance is
direct, not analogical. The original submission's Section 5 already framed
labour slack as exactly this closure: "the introduction of such a labor
slack by itself is to be conceived and interpreted as the advent of an
exogenous shock ... a spontaneous expansion of the labor force." Its
legitimizing narrative was borrowed from BF 2022's channel: "the context
of Corona-related lockdowns that partially required people to abstain from
work due to social distancing regulations" (with BF 2022's minimum-wage
unemployment explicitly set aside because it "does not relax capacity
constraints"). The code implements the manuscript's two variants:
\texttt{full\_labor\_slack} and \texttt{full\_labor\_slack\_alt}
(Leontief-projection-calibrated sectoral vectors; Section 5.1, aggregate
employment +1.7\%) and \texttt{empirical\_labor\_slack} (the uniform rule
$\mathbf{l}_{\mathrm{new}} = (1-\mu)^{-1}\,\mathbf{l}$ with $\mu$ the ILO
unemployment rate of 3.2\%; Section 5.2, +3.1\%). The referees' treatments
map onto these two variants: R1 questions 5.1's rationale ("unsure about
the rationale behind this scenario") and 5.2's missing wage-setting rule;
R2's "deus ex machina" targets the injection as the bridge mechanism and
the uniform transfer as arbitrary (R2.10). BF 2022 runs the same machinery
as calibrated contractions; the original manuscript ran it as calibrated
expansions. The revision therefore changes not the closure's existence but
its role: from the central bridge mechanism to a named scenario closure,
with the bridge carried by endogenous labour supply (Closure B), the wage
regime, and an explicit financing closure.}


## Coverage check: BF closures outside our lists \textcolor{revisionV1}{\normalsize [added v2]}

\textcolor{revisionV1}{Verified against the two papers
(\texttt{Dokumente/}) and the replications; it spans both foundation
lists. BF 2019 itself uses nothing
beyond Closure A: its factor supply is inelastic at the endowment ("the
supply of the factor is inelastic"), factors are non-reproducible goods,
and the paper explicitly makes "no allowance for amplification of shocks
through endogenous labor supply and capital accumulation" -- elastic
labour supply (Closure B) is named there as a channel they do not use.
BF 2022, however, uses four closure elements that our two foundation
lists do not cover:}

| BF 2022 element | Coverage status in our lists |
|---|---|
| Downward **nominal** wage rigidity, sector-specific ($N$ sticky labour factors) | Partially covered by Closure C (adapted, stripped down). Two differences: nominal floor vs.\ our real-wage floor (moot only because the static model has CPI = 1; it bites under deflation), and $N$ sector-specific floors vs.\ one economy-wide floor. The replication implements sticky sectoral labour but not the nominal/deflation side |
| Credit constraints / hand-to-mouth vs.\ Ricardian households | **Not covered.** Demand-side income-distribution closure (who spends the shock); our model has one representative household; ROADMAP leaves it out of scope |
| Zero lower bound / monetary policy | **Not covered.** Nominal closure; absent from the local replication (which keeps only nominal GDP accounting) and outside our static real model; relevant only if the paper's question shifts, per the verdict's 2019-vs-2022 criterion |
| Multiple primary factors (capital alongside labour; flexible capital vs.\ sticky labour in BF 2022) | **Not covered by a foundation table.** Our model and the 2019 replication are one-factor (value added as composite); ROADMAP section 5 carries it as the "optional extension" (mandatory only for policy-evaluation claims) |
| Firm capacity constraints (social distancing) | **Not covered.** BF 2022 models a supply shock as a labour-supply shift **or** a firm capacity constraint ("firms could be forced to operate at lower capacity ... due to social distancing") -- a rationing element outside our price-based static system |
| Reallocation frictions (adjustment costs, kappa) | **Not covered by a foundation table.** BF 2019 Online Appendix 2 ("Stuck Intermediates and Adjustment Costs", the kappa variants); ours is a static model -- `labor_closures.md` names convex employment-adjustment costs for a dynamic extension |

\textcolor{revisionV1}{None of these gaps blocks the current plan: the
first is exactly the "stripped-down 2022 closure" the verdict prescribes,
and the remaining three enter only via the optional extensions already
flagged in Stage 0 (scope decisions) and ROADMAP section 5. They are
listed here so the response letter can state precisely which BF elements
are and are not used. One replication-code nuance to verify before citing
the BF 2019 mobile variant as Closure A exactly: in the MATLAB code that
variant pins $w = 1$ (numeraire) and solves only the $N$ price equations,
so the aggregate labour constraint is implicit rather than imposed; the
paper itself states inelastic aggregate factor supplies. The exogenous
labour-endowment shift, though present in both BF papers, is a different
closure from the sticky wage: it moves supply without a wage response,
whereas \texttt{:fixed} moves the wage floor and lets employment respond.}

# Foundation II: Financial Closures \textcolor{revisionV1}{\normalsize [added v2]}

\textcolor{revisionV1}{The four admissible experiment types of ROADMAP
Phase 2. The current \texttt{demand\_shock} is an unfinanced add-on and
inadmissible by the ROADMAP's own rule until one principal and one
robustness variant are chosen; the explicit demand anchor introduced with
the financing closure is also expected to resolve the $\eta = 1$ stall
(scale indeterminacy) observed from multiplier 0.5 upward. Neither BF
paper supplies a fiscal financing closure (verdict); only the
compositional-shifter normalization has a BF 2022 analogue.}

## Selected \textcolor{revisionV1}{\normalsize [added v2]}

\textcolor{revisionV1}{The three financing variants carried into the
evaluation matrix: the preference-reallocation/tax-financed bracket as
principal, external debt as the third variant.}

| **Budget-neutral preference reallocation** | *Compositional household demand shift with the budget closed by construction* |
|---|---|
| **Budget identity** | $\tilde\beta_i = \beta_i d_i \,/\, \sum_j \beta_j d_j$; $Z(p) = \sum_j \tilde\beta_j p_j^{1-\sigma}$; $c_i^h = E^h \tilde\beta_i p_i^{-\sigma} / Z(p)$, so $\sum_i p_i c_i^h = E^h$ holds exactly |
| **Economic interpretation** | Compositional demand shift; the deduction runs through the **household's** consumption composition (other categories shrink within the fixed $E^h$); government spending untouched; must **not** be described as an autonomous investment multiplier (verdict warning) |
| **Cost / status** | Low: normalization is half-done already (Milestone A normalized household demand) |
| **Recommended role** | Robustness / diagnostic benchmark |
| **Mentioned by reviewer** | no (neither referee demands it; R1's demand-shock rationale question is adjacent) |
| **Mentioned in old docs** | `vertdict.md` (redesign section 1); `ROADMAP.md` (Phase 2); `definitive_guide.md` |
| **Used by BF** | Yes (analogous treatment) -- BF 2022 renormalizes compositional demand shifters so the budget holds exactly (\texttt{bf\_replication2/src/model.jl:442}); but as a demand shifter, not a financing device (see filename above) |

| **Tax-financed public investment** | *Autonomous green investment with an explicit tax counterparty* |
|---|---|
| **Budget identity** | $\sum_i p_i g_i = T$ with a lump-sum or specified tax rule; separate government/investment vector $g_i$ |
| **Economic interpretation** | Autonomous green investment financed by taxation |
| **Cost / status** | Medium: new $g_i$ vector plus budget equation; not started |
| **Recommended role** | **Principal** (recommendation) |
| **Mentioned by reviewer** | no (the financing gap was identified internally; the referees demand shock transparency only) |
| **Mentioned in old docs** | `vertdict.md`; `ROADMAP.md` (Phase 2); `PRELIMINARY-ASSESSMENT.md` (financing admissibility); `definitive_guide.md` (Part V) |
| **Used by BF** | No (verdict: neither model supplies the fiscal financing closure) |

| **Debt/foreign-financed net expenditure** | *Demand expansion with an intertemporal or external counterpart* |
|---|---|
| **Budget identity** | $\sum_i p_i g_i = B + F$ with a defined saving--investment or current-account rule |
| **Economic interpretation** | Net new expenditure against borrowing or the external balance. Static-horizon nuance: debt finance ($B$) is resource-equivalent to the lump-sum tax (the household's current consumption gives up the same resources; the distinction becomes real only in a dynamic setting), while foreign finance ($F$) is genuinely distinct -- resources arrive via the external balance (imports serve the investment) |
| **Cost / status** | Medium--high: additional institutional balance |
| **Recommended role** | Third variant (external debt); $B$ collapses into the tax variant at this horizon |
| **Mentioned by reviewer** | no |
| **Mentioned in old docs** | `vertdict.md`; `ROADMAP.md` (Phase 2); `definitive_guide.md` |
| **Used by BF** | No |

\textcolor{revisionV1}{Financing design (user decision): taking the
housing stimulus as vantage point and ruling out unfinanced manna, the
realistic financing lies between the preference-reallocation and
tax-financed variants, so the principal experiments are a bracket of the
two. A third variant with external debt (foreign financing, $F$) is
added: within this static real model it is the closest representation of
the endogenous-money intuition that a borrowed euro need not be withheld
from domestic consumption ex ante -- the purchasing power enters, and
the ex-post identity $I = S$ (here: the external deficit) records where
the resources came from. The model itself cannot represent endogenous
money: with the nominal side collapsed it displays only the ex-post
accounting, so the $F$ variant is an accounting open-economy closure,
not a monetary mechanism. Imports feedback: at the current pipeline the
imported content of the $g_i$ vector is fixed by sector import shares,
so the external variant absorbs the demand through net imports
mechanically; a behavioral feedback (relative prices shifting import
shares) would require an Armington-type extension, for which the
existing $\Omega_{\mathrm{dom}}/\Omega_{\mathrm{raw}}$ split is the
natural hook. Deduction routes recap: preference reallocation deducts
from the household's consumption composition, expenditure switching from
old government expenditure, tax financing from household income.
Recommendation: run the preference-reallocation/tax-financed bracket as
principal, foreign financing as the third variant; debt ($B$) collapses
into the tax variant at this horizon; expenditure switching remains
available as a secondary option. The full 5 x 3 evaluation matrix is
assembled in the next section.}

\textcolor{revisionV1}{Financing design (user decision): taking the
housing stimulus as vantage point and ruling out unfinanced manna, the
realistic financing lies between the preference-reallocation and
tax-financed variants, so the principal experiments are a bracket of the
two; external debt is added as the third selected variant. Deduction
routes recap: preference reallocation deducts from the household's
consumption composition, expenditure switching from old government
expenditure, tax financing from household income via the lump-sum tax.
Imports feedback: at the current pipeline the imported content of the
$g_i$ vector is fixed by sector import shares, so the external variant
absorbs the demand through net imports mechanically; a behavioral
feedback (relative prices shifting import shares) would require an
Armington-type extension, for which the existing
$\Omega_{\mathrm{dom}}/\Omega_{\mathrm{raw}}$ split is the natural
hook. The full 5 x 3 evaluation matrix is assembled in the next
section.}

## Dropped \textcolor{revisionV1}{\normalsize [added v2]}

\textcolor{revisionV1}{Regretfully retired as a separate variant:
expenditure switching. Pure debt financing ($B$) is not given its own
variant at all -- within the static horizon it is resource-equivalent to
the tax variant (the distinction becomes real only in a dynamic
setting), so it would only re-label the principal tax experiment.}

| **Expenditure-switching investment** | *Reallocation within a fixed public budget* |
|---|---|
| **Budget identity** | $\sum_i p_i g_i = T$ with $T$ raised by cutting other public spending |
| **Economic interpretation** | Public budget reallocated across sectors; the deduction runs through **old government expenditure**; household budgets untouched; no new purchasing power |
| **Cost / status** | Low--medium |
| **Recommended role** | Robustness alternative |
| **Mentioned by reviewer** | no |
| **Mentioned in old docs** | `vertdict.md`; `ROADMAP.md` (Phase 2); `definitive_guide.md` |
| **Used by BF** | No |

\textcolor{revisionV1}{The endogenous-money issue. The motivation for
the external-debt variant is the endogenous-money intuition: in a world
where credit money is created alongside the loan, a borrowed euro need
not be withheld from domestic consumption ex ante; the ex-post identity
$I = S$ holds by accounting, but it does not describe an ex-ante fund
that had to be taken from anyone. The static real model cannot represent
this mechanism itself: with the nominal side collapsed there is no
credit creation, no banking sector, no interest-rate channel -- the
model displays only the ex-post accounting. What the $F$ variant does
capture is the shadow of the intuition: purchasing power enters without
an ex-ante sacrifice of domestic consumption, and the external deficit
records where the resources came from. It should therefore be sold as an
accounting open-economy closure, not as a monetary mechanism. The
contrast between the tax-financed and external-debt variants is
nonetheless paradigmatically productive: it operationalizes, inside one
framework, the difference between a loanable-funds reading (domestic
saving must finance investment ex ante) and an endogenous-money reading
(finance precedes saving; the balance closes ex post) -- exactly the
kind of commensurability the revised paper set out to measure.}

# The Evaluation Matrix (5 x 3) \textcolor{revisionV1}{\normalsize [added v2]}

\textcolor{revisionV1}{Decision (user): three financial closures -- the
preference-reallocation and tax-financed bracket plus the external-debt
variant -- crossed with the five selected labour closures gives a 5 x 3
matrix of 15 model variants to evaluate. This is not overblown: each cell
is one solve of a 71-sector static system, computationally trivial; the
interpretive burden is managed by reporting one headline set per cell
(Tornqvist real GDP response, employment, external balance) with
sectoral detail in an appendix. Two design amendments: BF's row uses the
pre-registered endpoints $\eta \in \{0,1\}$ (the interpolated $\eta^{*}$
is retired, ADR-0010), and DELTA is a corner of the matrix rather than an
independent equilibrium row -- its cells are reached as
GAMMA plus Leontief technology. The matrix serves the paradigmatic
main point directly: it turns the search for commensurability into a
measured map, with financing closures setting the scale of the stimulus
and labour closures its transmission margin.}

\textcolor{revisionV1}{Expected signatures, pre-registered before any
Stage 2 run (headline quantities per cell; sectoral detail deferred):}

| Labour closure | Financing: preference reallocation | Financing: tax-financed | Financing: external debt $F$ |
|---|---|---|---|
| BF (friction, $\eta^{*}$) | Composition shift under friction: aggregate ~ flat; sectoral mix most distorted | Tax withdrawal with frictional reallocation; GDP response smaller than ALPHA | Externally financed demand, frictional allocation; import leakage plus misallocation |
| ALPHA (full employment) | Composition only; aggregate ~ 0 (the $\eta$-invariance diagnostic); sectoral reallocation | Crowding-out via wages/prices; aggregate ~ 0 (classic full-employment CGE) | Demand leaks to imports; wage rise moderated; external deficit = $F$ |
| BETA (elastic supply) | Composition only; employment ~ flat along the supply elasticity | Positive but modest; employment and wage share the adjustment ($\eta_s$-dependent) | Positive; import leakage dampens the domestic employment response |
| GAMMA (sticky real wage) | Composition shift; small positive employment response | The extensive-margin result: employment absorbs the shock (current: +19.3 pp at mult 10) | Employment absorbs the domestically-produced share; imports take the rest |
| DELTA (IO endpoint, corner) | IO multiplier on the compositional shift (value-added effects) | IO multiplier net of the tax withdrawal | Full IO multiplier with import-adjusted inverse; external deficit = $F$ |

\textcolor{revisionV4}{Measured values for every cell are in section 4.2. Three readings of the table above need correction: the BF row's label carries the retired interpolated $\eta$-star (the endpoints $\eta \in \{0,1\}$ are what runs); the parenthetical "current: +19.3 pp at mult 10" is a retired UNFINANCED result; and \"external deficit = F\" describes the gross external inflow, not the net position the model reports (section 4.2; the decomposition is in \texttt{paper/tables/matrix\_5x3\_v3\_flows.md}).}

\textcolor{revisionV1}{Limit relations across the matrix (reported as
results, not hidden as redundancy): BF at $\eta = 1$ coincides with
ALPHA, and BF at $\eta = 0$ is the immobile benchmark; BETA's high
$\eta_s$ limit approaches GAMMA's fixed-real-wage regime; and GAMMA plus
Leontief technology is DELTA. DELTA is therefore a corner, not an
independent row; BF's endpoints duplicate ALPHA and the immobile
benchmark. The coincidences are the bridge result: the same endpoint is
reached from different directions, which is precisely the
paradigmatic-geography claim of the revised paper.}

![Figure: Paradigmatic directions towards the short run. The horizontal
axis is the labour-supply-elasticity continuum (ALPHA at zero, BETA
intermediate, DELTA at the IO limit); BF is the orthogonal
reallocation-friction dimension anchored at the fixed-supply end; GAMMA
approaches the IO endpoint from the wage-rigidity direction (GAMMA plus
Leontief technology). Drafted after the user's sketch.](pictures/paradigmatic_directions.png)

# Stage 0: Conceptual decisions \textcolor{revisionV2}{\normalsize [reworked v3]} \textcolor{revisionV3}{\normalsize [section 4 since v4]}

\textcolor{revisionV2}{Resolved in this round: the labour-closure
selection (Foundation I: BF, ALPHA, BETA, GAMMA, DELTA), the financing
design (Foundation II: the preference-reallocation/tax-financed bracket
as principal, external debt as the third selected variant; expenditure
switching dropped, pure debt subsumed by the tax variant), and the
consolidation of the exogenous-shift/legacy closure as dropped. The
status matrix below completes the Stage 0 answer by locating every cell
of the evaluation matrix in the current implementation. Two items remain
open: the audit-trail disclosure decision (Stage 4) and the scope locks
(skill classes, open-economy depth, Stream C off the critical path).}

\textcolor{revisionV2}{The two paradigmatic wedges. The matrix has two
axes, and both are paradigmatic wedges. The first -- the labour-market
closure -- is the referees' axis: it sets the transmission margin of a
demand stimulus, from allocation friction (BF) through full employment
(ALPHA) and voluntary supply (BETA) to wage rigidity (GAMMA) and the IO
endpoint (DELTA). The second -- the financing closure -- runs through
the financing assumptions: it sets where the purchasing power of the
stimulus comes from and through which budget it closes (household
composition vs.\ household income vs.\ external balance), i.e.\ the
loanable-funds versus endogenous-money contrast. Recommended strength of
the second dimension: first-class as structure, bounded as claim. It is
a full matrix axis and appears in every results table; but in the static
real model it operates as accounting closures, so the paper claims the
financing contrast as measured sensitivity, not as a theory of money.}

\textcolor{revisionV2}{Status matrix (Stage 0 answer).} Symbols: ![](pictures/emoji/2705.png){width=9pt}
implemented and validated; ![](pictures/emoji/1f7e0.png){width=9pt} partially available -- the machinery
exists but needs admissibility re-runs, explicit framing, or hardening
(the missing piece is stated); ![](pictures/emoji/1f534.png){width=9pt} not implemented; ![](pictures/emoji/1f525.png){width=9pt} exists in the code
but is inadmissible under the ROADMAP gates (the unfinanced autonomous
shock). No cell is ![](pictures/emoji/2705.png){width=9pt} yet -- under the corrected admissibility standard that is
precisely the finding, and the matrix defines the Stage 1 work.

| Labour closure | F1: preference reallocation | F2: tax-financed | F3: external debt |
|---|---|---|---|
| BF (friction, $\eta^{*}$) | ![](pictures/emoji/1f7e0.png){width=9pt} geometric $\eta$ implemented; missing $\eta^{*}$ pre-registration and the explicit $\tilde\beta$ experiment framing | ![](pictures/emoji/1f534.png){width=9pt} needs the $g_i$ vector and the $T$ rule (the existing unfinanced autonomous shock is ![](pictures/emoji/1f525.png){width=9pt} as input) | ![](pictures/emoji/1f534.png){width=9pt} needs the external-account closure |
| ALPHA (full employment) | ![](pictures/emoji/1f7e0.png){width=9pt} normalized composition shock runs (the +0.067% $\eta$-invariance diagnostic); needs explicit framing and reporting | ![](pictures/emoji/1f525.png){width=9pt} the autonomous-demand version is unfinanced; the tax rule is missing | ![](pictures/emoji/1f534.png){width=9pt} |
| BETA (elastic supply) | ![](pictures/emoji/1f534.png){width=9pt} closure not implemented | ![](pictures/emoji/1f534.png){width=9pt} | ![](pictures/emoji/1f534.png){width=9pt} |
| GAMMA (sticky real wage) | ![](pictures/emoji/1f7e0.png){width=9pt} :fixed implemented; composition run needs a validation pass | ![](pictures/emoji/1f7e0.png){width=9pt} the +19.3 pp extensive-margin sweep exists but rests on the ![](pictures/emoji/1f525.png){width=9pt} unfinanced shock; must be regenerated under F2 | ![](pictures/emoji/1f534.png){width=9pt} |
| DELTA (IO endpoint) | ![](pictures/emoji/1f534.png){width=9pt} corner not implemented (fixed factor price + Leontief limit) | ![](pictures/emoji/1f534.png){width=9pt} | ![](pictures/emoji/1f534.png){width=9pt} |

\textcolor{revisionV2}{Stage 0 remaining items: (i) audit-trail
disclosure (brief mention vs.\ full appendix vs.\ none; recommended:
brief, per the supersession-note culture); (ii) scope locks: skill-class
disaggregation (mandatory only for policy-evaluation claims),
open-economy depth beyond the section 4.1 mapping (Armington extension
optional), B\&F replication off the critical path (confirmed by all
surviving plans).}

## Intermediate results: the cbase2 pipeline (v3 open-economy recalibration) \textcolor{revisionV3}{\normalsize [added v4]}

\textcolor{revisionV3}{The notebook pipeline \texttt{cbase2/} (see
\texttt{cbase2/documentation.md} for the neutral description and
\texttt{cbase2/process\_comments.md} for the dated observations) has
implemented and verified the v3 open-economy recalibration. Three
structural findings shaped it. First, the DELTA investigation showed that
the fixed-wage rows are indeterminate without marginal leakages: with the
import margin alone, the only solution is the degenerate corner E = 0,
because the import leak has no offsetting injection. Second, the
national-accounts identity S = I + X - M forces the saving rate: the
household's non-consumption share is the accounting partner of the
exogenous investment and export injections, so the v3 calibration derives
it from the data (s = 0.398; government share tau0 = 0.214, export share
0.422, investment share 0.162) instead of assuming it. Third, the
equilibrium formulation was corrected throughout: the mobile system keeps
N-1 clearing equations plus the CPI = 1 numeraire, and the omitted N-th
market is the residual external account (it is NOT Walras-redundant once
imports leak, so it is exposed and asserted rather than silently dropped;
ADR-0010); the F1/F3 tax is the price-indexed baseline budget T = p'gG
(F3 reads exactly as "baseline tax unchanged, programme externally
financed"), and F2 keeps the genuinely balanced-budget rule.}

\textcolor{revisionV3}{Verified results (acceptance test
\texttt{v2\_verify.jl}): the injection continuation runs in 5.4 s with no
stalls; the price-explosion branch found by the direct solve is
eliminated; the F1/F2/F3 mobile rows solve to residuals of at most
3.7e-7 with EXACT budget identities (sum p c = (1-s)E); the F3 external
balance is recorded (F = 0.00147 at the m = 1 programme); and the DELTA
analytic equivalence is EXACT for both F2 and F3 (rel y error 0.0;
L = 1.1656 and 1.1767 at the m = 1 programme). The updated status matrix,
same layout as above:}

| Labour closure | F1: preference reallocation | F2: tax-financed | F3: external debt |
|---|---|---|---|
| BF (endpoint selector, $\eta \in \{0,1\}$) | ![](pictures/emoji/1f7e0.png){width=9pt} machinery verified; the explicit $\tilde\beta$ tilt ported to the canonical pipeline ($d = 1 + G0\cdot\psi/c0$, \texttt{f1\_shift = tilt\_g0\_over\_c0}); matrix-cell runs pending | ![](pictures/emoji/2705.png){width=9pt} implemented and verified (balanced-budget rule; resid 1.1e-7) | ![](pictures/emoji/2705.png){width=9pt} implemented and verified (external balance recorded) |
| ALPHA (full employment) | ![](pictures/emoji/1f7e0.png){width=9pt} v3 re-run pending ($\eta = 1$) | ![](pictures/emoji/1f7e0.png){width=9pt} v3 re-run pending ($\eta = 1$) | ![](pictures/emoji/1f7e0.png){width=9pt} v3 re-run pending ($\eta = 1$) |
| BETA (elastic supply) | ![](pictures/emoji/1f7e0.png){width=9pt} implemented ($\eta_s$, continuation solver; elasticity identification verified at the v2 stage); v3 verification pending | ![](pictures/emoji/1f534.png){width=9pt} verification pending | ![](pictures/emoji/1f534.png){width=9pt} verification pending |
| GAMMA (sticky real wage) | ![](pictures/emoji/1f534.png){width=9pt} v3 run pending | ![](pictures/emoji/1f7e0.png){width=9pt} the DELTA corner (its Leontief limit) verified exact; the CES-elasticity GAMMA run pending | ![](pictures/emoji/1f534.png){width=9pt} v3 run pending |
| DELTA (IO endpoint, corner) | ![](pictures/emoji/1f534.png){width=9pt} F1 cell pending | ![](pictures/emoji/2705.png){width=9pt} implemented and verified -- EXACT analytic equivalence (rel y error 0.0) | ![](pictures/emoji/2705.png){width=9pt} implemented and verified -- EXACT analytic equivalence (rel y error 0.0) |

\textcolor{revisionV3}{Open items carried into Stage 1: the BETA
verification (the $\eta$-continuation runtime under v3), the ALPHA and GAMMA
v3 re-runs at $\eta = 1$ and the CES elasticities, the F1 tilted-cell runs
(the explicit tilt is ported to the canonical pipeline as
\texttt{f1\_shift = tilt\_g0\_over\_c0}; the matrix cells are pending),
and the interpretive pass on the v3 baseline units (the no-numeraire
units experiment is retired with ADR-0010: the CPI = 1 numeraire is
restored, so the CPI-normalized wage is the real wage; the 0.54 reading
belonged to the retired no-numeraire units and must be re-checked in
notebook 04).}

\textcolor{revisionV3}{Post-v4 findings (folded into this version, no
version change). Three further results qualify the intermediate state.
(i) The branch mystery was resolved: the "depressed" equilibrium found by
the solver is the same real allocation at a lower price scale; the fixed
nominal tax made its real burden price-level-dependent, so the government
budget was made homogeneous (real purchases $gG$ financed at current
prices, $T = p'gG$) and the CPI numeraire restored as the scale selector.
(ii) Sector 71 was dropped on documented grounds (residual catch-all,
1.8\% of gross output, 37.5\% self-loop), but the instability then moved
to other high-self-loop service sectors -- the (f) provenance diagnostic
showed the diagonals are genuine data (self-shares identical before and
after the import split; 23 sectors above a 0.3 self-share), so no data
correction is available and the trade-off is structural. Three dataset
variants are now defined in \texttt{cbase2/src/calibration.jl}
(\texttt{full} = 71 sectors, \texttt{70s} = drop 71, \texttt{reduced} =
additionally drop the self-share $> 0.45$ class 13, 18, 19, 48, 53, 58,
68; coverage: 100\%/99.2\%/87.4\% of gross output).
(iii) The solver question is reframed: the theta-homotopy showed that
even at $\theta = 2$ (gross substitutes, no self-referencing spiral
possible) the FD-Newton stalls at the same $3.6\times10^{-4}$ floor --
the difficulty is solver-setup reliability (kinked residuals under
ForwardDiff), not model structure; the identical system converged to
$2\times10^{-10}$ under a different init. This matches the Baqaee--Farhi
practice of solving equilibria as general NLPs (KNITRO/fmincon in
MATLAB): the recommended cbase2 route is a constrained-NLP formulation
(IPOPT via JuMP, the free KNITRO analogue) if the quick
finite-difference/least-squares experiments do not suffice.}


\textcolor{revisionV4}{The v4 numbers in this section are the CLAMPED cbase2 calibration and are superseded. The A-bill calibration (ADR-0012, ADR-0013) re-anchors the saving rate to s = 0.1199 (not 0.398), deletes the clamp entirely (the household residual is the table's own domestic household column, with zero negative sectors), books the product taxes on intermediate use (raw-table row 75) as the third external leak, and makes the omitted-market canary identity exact at 1.7e-16 on full-71 (9.6e-17 on the 70s variant). The fixed-wage indeterminacy reported above is superseded too: ADR-0014 replaced the heuristic admissibility guard with the verified round-gain criterion and the open calibration is determinate (full-rank Jacobian; three inits converge to one root). The raw table's own production-vs-expenditure residual (5.387 \%) is unrelated and remains open.}

## Matrix results: the 5 x 3 evaluation on the A-bill calibration \textcolor{revisionV4}{\normalsize [added v5]}

\textcolor{revisionV4}{Design \texttt{matrix\_5x3\_v3} (the ADR-0014/ADR-0015 generation): 15 of 15 cells executed, no locked and no dropped scenarios. Below is the pre-registered headline set -- Tornqvist real GDP against the baseline, employment, and the model's own identity-consistent external position (the omitted-market canary, S - (I+X-M) + T). The per-cell accounting decomposition (saving, tax, investment plus exports, and the three import components) is in \texttt{paper/tables/matrix\_5x3\_v3\_flows.md}; the paper cites the \texttt{matrix\_5x3-v3-*} run ids.}

| Cell | Real GDP rel. | Employment | Net external position |
| --- | ---: | ---: | ---: |
| `BF-F1` | +0.000433 | 1.000000 | +0.000453 |
| `BF-F2` | -0.016936 | 1.000000 | -0.000548 |
| `BF-F3` | -0.000000 | 1.000000 | -0.007947 |
| `ALPHA-F1` | +0.000433 | 1.000000 | +0.000465 |
| `ALPHA-F2` | -0.016936 | 1.000000 | -0.000601 |
| `ALPHA-F3` | +0.000000 | 1.000000 | -0.008495 |
| `BETA-F1` | +0.000433 | 1.000000 | +0.000465 |
| `BETA-F2` | -0.016936 | 1.000000 | -0.000601 |
| `BETA-F3` | +0.000000 | 1.000000 | -0.008495 |
| `GAMMA-F1` | -0.000806 | 0.999027 | -0.000000 |
| `GAMMA-F2` | -0.015333 | 1.001260 | -0.000000 |
| `GAMMA-F3` | +0.022644 | 1.017796 | -0.000000 |
| `DELTA-F1` | -0.000806 | 0.999027 | -0.000000 |
| `DELTA-F2` | -0.015333 | 1.001260 | -0.000000 |
| `DELTA-F3` | +0.022644 | 1.017796 | -0.000000 |

\textcolor{revisionV4}{Findings against the pre-registered signatures -- deviations are findings, not failures (ADR-0004):}

1. \textcolor{revisionV4}{The F1 column is admissible in the fixed-wage rows. \texttt{GAMMA-F1} and \texttt{DELTA-F1} were blocked by a heuristic guard whose premise is false on the open calibration: the Jacobian at the F1 point is full rank (sigma-min/sigma-max = 0.1586) and three different inits converge to the same root. The verified criterion (round-gain column sums below 1) admits them and still rejects closed fixtures, where the sums are exactly 1. The same criterion exposed a latent bug: the retired heuristic had admitted closed-fixture solves whenever manna was present, although manna is a constant and cannot remove a unit root.}
2. \textcolor{revisionV4}{BETA reproduces ALPHA to machine precision. Under a demand-only programme the price block is demand-invariant (p = 1), so the equilibrium real wage sits at its anchor and L = Lbar * $(w/w_0)^{\eta_s}$ = Lbar for any $\eta_s$: the supply elasticity is unidentified by THIS experiment, not mis-specified. It bites under a supply-side scenario (a +20 \% sector-1 shock moves the real wage to 1.0049 and separates $\eta_s$ = 0.5 from $\eta_s$ = 2: L = 1.0024 versus 1.0098).}
3. \textcolor{revisionV4}{The F2 column is not "aggregate ~ 0" for the mobile rows: real GDP falls 1.69 \% with employment pinned at exactly 1.0 -- a composition effect at fixed aggregate labour input, for a tax withdrawal of 1.33 \% of GDP.}
4. \textcolor{revisionV4}{The limit relations hold: BF at $\eta$ = 0 coincides with ALPHA in F1 and F2 at the aggregate, and GAMMA coincides with DELTA to 5e-12 in F2 and F3. F3 is the extensive-margin column -- employment rises 1.78 \% in the fixed-wage rows.}
5. \textcolor{revisionV4}{The external account: the F3 programme's gross inflow is F = 0.0133 of GDP, of which 0.0029 is its own import content; the NET position is -0.0085 in the mobile rows (the endogenous response offsets 36 \% of the inflow) and zero to machine precision in the fixed-wage rows, where all N markets clear and the inflow is absorbed by imports.}

\textcolor{revisionV4}{\textbf{Recombination option (left out).} The matrix varies two labour margins that the design parametrises independently: BF's allocation selector $\eta \in \{0,1\}$ (who gets the labour) and BETA's aggregate supply elasticity $\eta_s$ (how much labour in total). The 5 x 3 samples three of the four corners -- ALPHA ($\eta$ = 1, $\eta_s$ = 0), BETA (1, 0.5) and BF (0, 0) -- and leaves out their recombination, the corner ($\eta$ = 0, $\eta_s$ = 0.5): immobile allocation WITH elastic supply. It is left out deliberately because it would be a null result here: with p = 1 the real wage sits at its anchor, so the recombined cell coincides with BF to about 1e-15, exactly as BETA coincides with ALPHA. It becomes informative only in the supply-side scenario of finding 2, where a ($\eta$, $\eta_s$) 2 x 2 would separate the two margins -- the allocation friction moving composition and the supply elasticity moving the level. Recorded as an available extension, not as a gap in the current design.}

\textcolor{revisionV5}{Version 6 measurement update (ADR-0018). The v5 table above reports the household-consumption Tornqvist index in its "Real GDP rel." column; under ADR-0018 that index is the welfare aggregator, not GDP. Real GDP is the income-side index w*sum(L) deflated by the seven-component GDP deflator, which is exactly one here because the design pins all prices at one. The matrix\_5x3\_v4 generation records both concepts; the table below is taken from the 15 v4 manifests. Finding 3 above is corrected: the F2 mobile rows are GDP-neutral (0.000 percent), and the -1.69 percent is the tax-induced welfare loss.}

| Cell | Real GDP rel. | Consumption rel. | Employment | Net external position |
| --- | ---: | ---: | ---: | ---: |
| `BF-F1` | +0.000000 | +0.000433 | 1.000000 | +0.000452 |
| `BF-F2` | +0.000000 | -0.016936 | 1.000000 | -0.000548 |
| `BF-F3` | +0.000000 | +0.000000 | 1.000000 | -0.007947 |
| `ALPHA-F1` | +0.000000 | +0.000433 | 1.000000 | +0.000465 |
| `ALPHA-F2` | +0.000000 | -0.016936 | 1.000000 | -0.000601 |
| `ALPHA-F3` | +0.000000 | +0.000000 | 1.000000 | -0.008495 |
| `BETA-F1` | +0.000000 | +0.000433 | 1.000000 | +0.000465 |
| `BETA-F2` | +0.000000 | -0.016936 | 1.000000 | -0.000601 |
| `BETA-F3` | +0.000000 | +0.000000 | 1.000000 | -0.008495 |
| `GAMMA-F1` | -0.000973 | -0.000806 | 0.999027 | -0.000000 |
| `GAMMA-F2` | +0.001260 | -0.015333 | 1.001260 | -0.000000 |
| `GAMMA-F3` | +0.017796 | +0.022644 | 1.017796 | -0.000000 |
| `DELTA-F1` | -0.000973 | -0.000806 | 0.999027 | -0.000000 |
| `DELTA-F2` | +0.001260 | -0.015333 | 1.001260 | -0.000000 |
| `DELTA-F3` | +0.017796 | +0.022644 | 1.017796 | -0.000000 |

\textcolor{revisionV5}{The full per-cell accounting decomposition -- the seven GDP components (C, G plus programme, I, X, M\_final, M\_int, T\_int), the external wedge Sigma(V) - w*sum(L), and the deflator -- is in \texttt{paper/tables/matrix\_5x3\_v4\_flows.md}. Findings 1, 2, 4 and 5 of the v5 list stand: the limit relations, the unidentified supply elasticity under demand-only shocks, and the external-account reading are unchanged, and the v4 runs reproduce the v3 equilibria to 1.2e-11, so only the measurement changed. The recombination option remains a deliberate omission.}

\textcolor{revisionV6}{\normalsize [added v7: the ADR-0019 generation]}

\textcolor{revisionV6}{Design \texttt{matrix\_5x3\_v5} (the ADR-0019 generation): 15 of 15 cells executed under the all-N closure with the explicit external transfer, no locked and no dropped scenarios. The headline set below adds the booked external account (F, B\_gov, net position F + B\_gov) to the v6 columns; the per-cell accounting decomposition is in \texttt{paper/tables/matrix\_5x3\_v5\_flows.md}, and the paper cites the \texttt{matrix\_5x3-v5-*} run ids for the twelve non-BF rows and the \texttt{matrix\_5x3-v6-*} run ids for the BF row, whose eta = 0 closure changed after this table was written (marker below).}

| Cell | Real GDP rel. | Consumption rel. | Employment | F | B_gov | Net ext. pos. |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `BF-F1` | +0.000000 | +0.000433 | 1.000000 | +0.00000000 | 0.00000000 | +0.00000000 |
| `BF-F2` | +0.000000 | -0.016936 | 1.000000 | +0.00000000 | 0.00000000 | +0.00000000 |
| `BF-F3` | +0.000000 | +0.000000 | 1.000000 | +0.00000000 | 0.01330991 | +0.01330991 |
| `ALPHA-F1` | +0.000000 | +0.001432 | 1.000000 | +0.00078481 | 0.00000000 | +0.00078481 |
| `ALPHA-F2` | +0.000000 | -0.018226 | 1.000000 | -0.00101403 | 0.00000000 | -0.00101403 |
| `ALPHA-F3` | +0.000000 | -0.018226 | 1.000000 | -0.01432394 | 0.01330991 | -0.00101403 |
| `BETA-F1` | +0.000000 | +0.001432 | 1.000000 | +0.00078481 | 0.00000000 | +0.00078481 |
| `BETA-F2` | +0.000000 | -0.018226 | 1.000000 | -0.00101403 | 0.00000000 | -0.00101403 |
| `BETA-F3` | +0.000000 | -0.018226 | 1.000000 | -0.01432394 | 0.01330991 | -0.00101403 |
| `GAMMA-F1` | -0.000973 | -0.000806 | 0.999027 | +0.00000000 | 0.00000000 | +0.00000000 |
| `GAMMA-F2` | +0.001260 | -0.015333 | 1.001260 | +0.00000000 | 0.00000000 | +0.00000000 |
| `GAMMA-F3` | +0.017796 | +0.022644 | 1.017796 | +0.00000000 | 0.01330991 | +0.01330991 |
| `DELTA-F1` | -0.000973 | -0.000806 | 0.999027 | +0.00000000 | 0.00000000 | +0.00000000 |
| `DELTA-F2` | +0.001260 | -0.015333 | 1.001260 | +0.00000000 | 0.00000000 | +0.00000000 |
| `DELTA-F3` | +0.017796 | +0.022644 | 1.017796 | +0.00000000 | 0.01330991 | +0.01330991 |

\textcolor{revisionV6}{Marker (same date, within this version): the three BF rows of the table above are the SUPERSEDED pin values -- the F = 0 pin of ADR-0019 D2 leaves the external position unidentified, so the entry collapses to the programme cost B\_gov = 1.331 percent of GDP and the account stands open by -0.794 percent at BF-F3. The current BF results (sectoral wages at eta = 0, F solved, identified external position, closed account) are in the promotion block below and in \texttt{paper/tables/matrix\_5x3\_v6\_flows.md}; the twelve non-BF rows of this table, GAMMA and DELTA included, are unchanged in v6.}

\textcolor{revisionV6}{Three earlier readings are superseded by ADR-0019. First, the pre-registered external deficit = F signature: F3's gross programme inflow is B\_gov = 1.331 percent of GDP while the net position is F + B\_gov (-0.101 percent in the mobile rows, +1.331 percent in the fixed-wage rows, and not identified at the BF eta = 0 rows; see the correction below). Second, the wedge-as-external-position reading: gdp\_wedge is now minus the identity gap (about 1e-16 at the mobile eta = 1 cells, about 1e-12 at the fixed-wage eta = 1 cells at residual level, and the documented factor-market gap of +4.34e-04 / -5.62e-04 / -7.94e-03 at the BF eta = 0 cells), and the external position is the booked F + B\_gov. Third, the F2/F3 separation: the two columns are real-neutral in the mobile regime (identical consumption at -1.823 percent and identical net position to about 2.6e-17, with F\_F3 = F\_F2 - B\_gov), so their v1-v4 differences were artifacts of the omitted-market shortcut. Findings 1, 2 and 4 of the v5 list stand unchanged, finding 3 stands in its v6 corrected form, and the recombination option remains a deliberate omission.}

\textcolor{revisionV6}{Correction (same date, within this version): the BF eta = 0 external-position entries are booked financing, not positions. The eta = 0 endpoint pins F = 0 (ADR-0019 D2), so the reported entry F + B\_gov collapses to the programme cost B\_gov = sum(p g) = 1.331 percent of GDP -- identical to the last digit in BF-F3, GAMMA-F3 and DELTA-F3, i.e. a design constant that no labour closure can move (GAMMA/DELTA-F3 close their account to about 1e-12, BF-F3 does not). The same cell's resource side S + T\_int + M - (I+X) reads +0.5374 percent of GDP at BF-F3 against the booked +1.3310 percent, the account open by -0.7936 percent of GDP; that gap is exactly the labour-market residual the pin leaves, -w (sum L\_cm - sum L), with sum L\_cm = 1.0079361 against the frozen bar 1 at BF-F3 (verified at all 15 cells; at the twelve eta = 1 cells the two quantities agree to machine precision). The position is therefore not identified at eta = 0: replacing the pin by F = c leaves every c an exact root (residual 2.2e-15) and moves the entry one-for-one, from -0.669 to +3.331 percent of GDP for c in -0.02 to +0.02, with the resource side moving too. Solving the eta = 0 system with the pin at ALPHA's solved F reproduces the ALPHA cell to 2.8e-17 in the canonical vector and in the consumption block, and replacing the pin by the labour equation sum L\_cm = Lbar closes the account to 4e-16 and returns BF-F1/F2/F3 = ALPHA-F1/F2/F3 exactly (with financing neutrality then holding at eta = 0 as well), so in this demand-only design the BF/ALPHA difference -- including the BF-F3 welfare entry of 0.000 percent, which becomes -1.8226 percent at the labour-clearing F -- is carried by the pin, not by the frozen allocation. Read the BF rows from the resource side and the gap in \texttt{paper/tables/matrix\_5x3\_v5\_flows.md}; the eta = 0 closure options and their scope for a v6 generation are in ADR-0020 (accepted: option C, with the specification and probe measurement at the end of this section).}

\textcolor{revisionV6}{Sources of the recorded deviations (classification). None of the recorded deviations was an accounting error: the identities (zero-profit, budgets, the booked external account) hold to machine precision throughout. They fall into four classes. (i) Measurement--definition: the legacy "Real GDP" column is the household-consumption welfare index; corrected by the income-side measure (ADR-0018, version 6). (ii) Closure artifact: the F2/F3 separation and the earlier external-position readings came from the omitted N-th market; corrected by all-N clearing with the booked external account (ADR-0019, this version; the approaches considered are catalogued in its 2026-09-18 amendment), under which the mobile F2 and F3 cells are real-neutral. (iii) Model-class property: BETA coincides with ALPHA under demand-only shocks because the one-factor CRS price block is demand-invariant (Samuelson nonsubstitution), and DELTA's exactness rests on the same p = 1; eta\_s is identifiable only under a supply-side scenario or with a second production factor -- a property of the closure design, not a mis-specification. (iv) Corrected admissibility mathematics: the GAMMA/DELTA-F1 rejections were a retired heuristic guard, replaced by the measured actual-matrix criterion (ADR-0014, ADR-0017).}

\textcolor{revisionV6}{Decision (same date, within this version): ADR-0020 is accepted with option C -- sectoral wages at the eta = 0 endpoint -- so the BF row's external position is repaired in the closure rather than patched by convention. The frozen sectoral allocation stays in place and each sector's wage is set by that sector's own marginal product at the frozen allocation: unknowns [p; y; w(1:N); F], 3N + 1 in all, with the N zero-profit conditions p\_i = cost\_i(p, w\_i), the N sectoral first-order conditions at the frozen allocation, the N clearing equations with household wage income the sum of w\_i L\_i, and the CPI numeraire. F is kept because the block of zero-profit, FOC and clearing conditions is homogeneous of degree 1 in (p, w, F): N sectoral conditions replace the mobile system's single aggregate labour equation, so the block carries one equation more than it has directions to pin and the demand side needs one free scalar; dropping a clearing equation is the ADR-0010 shortcut that ADR-0019 retired. Two consequences, both measured: the external account closes at eta = 0, because the identity gap is exactly the value-weighted difference between the cost-minimizing and the frozen labour allocation and that difference is zero when the frozen allocation is cost-minimizing at these wages; and financing neutrality holds at eta = 0 as well (F\_F3 = F\_F2 - B\_gov, with an identical booked position). Option B (labour equation instead of the pin, common wage) is recorded in ADR-0020 as a future venue: it also closes the account, but it collapses the BF row onto ALPHA at 2.8e-17, so it buys accounting coherence by giving up the immobile benchmark as a distinct economy; option A is superseded by this choice.}

\textcolor{revisionV6}{Probe measurement (full-71 A-bill calibration, before any src change; \texttt{experiments/probes/probe7\_sectoral\_wages\_eta0.jl}):}

| Cell | max abs residual | identity gap | F | B_gov | booked F + B_gov |
| --- | ---: | ---: | ---: | ---: | ---: |
| `BF-F1` | 1.4e-13 | +1.9e-14 | -0.00581497 | 0.00000000 | -0.00581497 |
| `BF-F2` | 4.2e-13 | -4.2e-13 | -0.00852841 | 0.00000000 | -0.00852841 |
| `BF-F3` | 1.7e-12 | -5.0e-13 | -0.02337469 | +0.01484628 | -0.00852841 |

\textcolor{revisionV6}{The eta = 0 row is therefore no longer a normalisation of ALPHA: it carries its own wages (0.970 to 1.672 around the CPI-pinned level, wage bill 1.00702), its own prices and an identified external position (F-column norm 0.163; Jacobian condition number 9.4e7 against 52.8 for the mobile all-N system, so the v6 acceptance gate and the polish target must be set from measurement rather than copied). Because the F3 booking is priced at the eta = 0 prices, which the sectoral wages move, it reads +1.4846 percent of GDP there against the baseline +1.3310 percent. Scope of the change: it is the eta = 0 branch only. ALPHA and BETA are the eta = 1 mobile system and GAMMA and DELTA are the fixed-wage regime (w = 1 pinned, no labour equation, no F) -- a different system from \texttt{problem} -- so those rows' v5 results are unaffected by design; and the one shared helper the change touches (\texttt{\_cost\_minimizing\_labor}, whose scalar-wage form is log(w)) becomes generic through log.(w), a no-op at a scalar wage, which the promotion must verify by re-running the v5 cells and asserting the twelve non-BF rows reproduce to machine precision. Promotion into \texttt{src/} and the \texttt{matrix\_5x3\_v6} generation is the next step and is not yet done.}

\textcolor{revisionV6}{Promotion (same date, within this version): option C is in the kernel and executed as the \texttt{matrix\_5x3\_v6} generation (15 of 15 cells executed, commit \texttt{6db1da5}). The change is dispatch-only: \texttt{\_cost\_minimizing\_labor} evaluates \texttt{log.(w)} and so accepts a wage vector, \texttt{\_wage\_bill} keeps the scalar arithmetic bit-identical while summing w\_i L\_i for a vector, \texttt{problem\_sectoral} is the new eta = 0 system (3N + 1 unknowns), \texttt{problem} rejects eta = 0, and \texttt{solve} dispatches the new system with a tighter primary tolerance (1e-8) and a longer polish ladder (6 steps, target 1e-13) because it is stiff. Verification: the twelve non-BF cells reproduce their v5 manifests -- ALPHA/BETA to 6.7e-16 and GAMMA/DELTA to 1.6e-11, where the same differences appear on the pristine kernel (stash test, \texttt{probe9\_nonbf\_reproduction.jl}), so they are the fixed-wage sensitivity to the warm start (the system is well conditioned: Jacobian cond about 6.3 measured on the GAMMA-F2 cell, so the amplification is in the measurement layer, ~100x the residual level) rather than an effect of the change. The eta = 0 identity gate is set from the measured floor (harness 1e-9, test 1e-10; worst 1.04e-11 at BF-F1) instead of the 1e-12 written in ADR-0020's plan text. The BF row now reads:}

| Cell | F | B\_gov | Net ext. pos. | Gap | Consumption rel. |
| --- | ---: | ---: | ---: | ---: | ---: |
| \texttt{BF-F1} | -0.005815 | 0 | -0.005815 | -1.0e-11 | +0.096 percent |
| \texttt{BF-F2} | -0.008528 | 0 | -0.008528 | -1.8e-13 | -1.980 percent |
| \texttt{BF-F3} | -0.023375 | +0.014846 | -0.008528 | -1.9e-13 | -1.980 percent |

\textcolor{revisionV6}{against the v5 pin generation (net external position 0 / 0 / +1.331 percent of GDP, gap +4.3e-04 / -5.6e-04 / -7.9e-03, consumption rel. +0.043 / -1.694 / 0.000 percent). Financing neutrality now holds at eta = 0 as well (\texttt{F\_F3 = F\_F2 - B\_gov} to 3.5e-14, with identical positions), the BF/ALPHA coincidence is gone, and the paper-facing table is \texttt{paper/tables/matrix\_5x3\_v6\_flows.md}. The correction above describes the superseded v5 generation: cite the \texttt{matrix\_5x3-v6-*} run ids. A harness episode belongs to the record: the first v6 batch recorded the three BF cells as failed on \texttt{MethodError: isfinite(::Vector\{Float64\})}, because the \texttt{wage} metric carried the sectoral wage vector at eta = 0; the metric now reports the wage-bill-weighted aggregate (dispersion in the new \texttt{wage\_min} / \texttt{wage\_max} diagnostics), the three crash directories were removed and the whole generation re-run so that all fifteen cells come from one harness version. The eta = 0 FOC row also had a real defect before that: using the frozen allocation in the FOC made those N equations identically zero and the system degenerate, which the external-account gate caught immediately (gap 3.3e-04) -- the acceptance test earned its keep.}

\textcolor{revisionV6}{The same generation's fixed-wage rows are unchanged, so the current GAMMA and DELTA results remain those of the v7 headline table above: the promotion was verified row by row (ALPHA/BETA reproduce v5 to 6.7e-16, GAMMA/DELTA to 1.6e-11, the fixed-wage warm-start sensitivity recorded below). For reference, the current GAMMA row is:}

| Cell | Real GDP rel. | Consumption rel. | Employment | F | B_gov | Net ext. pos. |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `GAMMA-F1` | -0.000973 | -0.000806 | 0.999027 | 0.00000000 | 0.00000000 | 0.00000000 |
| `GAMMA-F2` | +0.001260 | -0.015333 | 1.001260 | 0.00000000 | 0.00000000 | 0.00000000 |
| `GAMMA-F3` | +0.017796 | +0.022644 | 1.017796 | 0.00000000 | 0.01330991 | +0.01330991 |

\textcolor{revisionV6}{These are the values to cite for the fixed-wage rows. DELTA reproduces them to six digits (its Leontief corner), the external account closes to about 1e-12, i.e. at the level of the solver residual, and F3 is the extensive-margin column: employment rises 1.78 percent and no external transfer is involved (the fixed-wage regime has no external unknown).}

\textcolor{revisionV6}{\normalsize [added v7, same date: variation in GAMMA]}

\textcolor{revisionV6}{The fixed-wage rows carry one degree of freedom in the matrix. F2 and F3 are real-neutral by theorem (financing neutrality) and DELTA reproduces GAMMA to six digits, so the three GAMMA cells are effectively two economies; and the flatness is a property of the design rather than of the model class, because the closure's parameter space already contains a wage vector that the design sets to its degenerate value (all sectors equal) and the kernel hard-codes as the scalar pin. ADR-0021 (proposed, same date) would make that vector a parameter of the existing closure, with the wage structure as a preregistered design dimension instead of a new closure; the options, their costs and the claim-by-claim provenance are in \texttt{docs/VariationinGamma.md}.}

\textcolor{revisionV6}{Measured on the GAMMA-F2 cell (full-71 A-bill, \texttt{experiments/probes/probe11\_gamma\_wage\_structure.jl}, 2026-09-18; a prototype that mirrors the fixed-wage residual, so the magnitudes await a kernel run): the degenerate pin reproduces the cell (residual 8.9e-16, employment 1.00125981). A UNIFORM rescale of the pin is a numeraire change and moves nothing real -- wbar = 0.9 and 1.1 leave the allocation, employment and real income invariant to 8.9e-16, with the CPI scaling exactly -- so the wage LEVEL is not an instrument. A RELATIVE change is a genuine shock: a ten percent wage push in the seven programme sectors moves employment by about half a percent and real income by about one percent (1.00686452 and 1.01186464 against 1.00125981), a single-sector tilt moves the allocation by about two percent in log terms, and a broad tilt of sectors 1 to 5 by about 3.5 percent -- larger than the closure's entire existing F1 to F2 spread. The external account still closes at the solver-residual level under every tilt (price-weighted clearing residual 2.9e-16 at the degenerate pin, at most 7.0e-12 tilted), validated against the kernel canary along quantity perturbations (ratio 1.0000 at one and ten per mille). The Jacobian condition number stays between 6.32 and 6.34, so the fixed-wage system is well conditioned at this calibration and its metric sensitivity to the warm start is a measurement-layer amplification of the residual level rather than a singular solve.}

\textcolor{revisionV6}{Three limits belong with these numbers. First, the pinned vector is an exogenous scenario parameter, never an equilibrium object: the variation is an instrument response, like the programme incidence or the shock size, and the scenarios must be preregistered before any run that reports them. Second, the level invariance is closure-conditional -- it holds because the model has no nominal anchor besides the wage, and it would fail with a nominal government bill, nominal debt or money. Third, the magnitudes come from a mirrored residual: the ADR-0020 promotion showed that translation into the kernel is where defects appear, so they should be reproduced through the kernel before citation, which also requires extending \texttt{external\_balance\_canary} and \texttt{gdp\_components} so that the acceptance gate evaluates a tilted pin on kernel arithmetic rather than on a proxy. The zero-cost alternative, needing no kernel change, is a sectoral supply-shock design, which also makes the allocation margin bite. The two routes it leaves open are recorded the same date: the thematic supply-side version (a climate-change scenario with a sectoral incidence vector) as \texttt{docs/ideas/IDEA-0001-climate-productivity-shocks.md} in the new ideas register, and the demand-side route (a free wage in a segmented wage structure, which is what makes prices demand-sensitive) as \texttt{docs/WORKPLAN\_SENSITIVE\_PRICES.md}. Measured while scoping the latter: the eta = 0 row is already the only row with demand-sensitive prices (max abs(p - 1) of 0.221 to 0.279 against 2.7e-15 in the ALPHA and BETA rows), because its sectoral wages are solved given demand, while a pinned wage vector produces price heterogeneity without price sensitivity.}

# Workplan \textcolor{revisionV3}{\normalsize [section 5 since v4; reworked v3]}

## Stage 1: Model completion \textcolor{revisionV2}{\normalsize [reworked v3]}

\textcolor{revisionV2}{Ordered by the status matrix -- the financing
core first (it makes every experiment admissible), then the missing
labour closures, then the guards:}

1. **Financing core**: the shared final-demand and institutional-budget
   core; the government/investment vector $g_i$ with the tax rule
   $\sum_i p_i g_i = T$ (F2); the external-account closure with
   $\sum_i p_i g_i = F$ and import content at fixed sector shares (F3);
   the explicit preference-reallocation experiment with renormalized
   weights $\tilde\beta_i$ (F1); retire the unfinanced autonomous
   shock (![](pictures/emoji/1f525.png){width=9pt}).
   \textcolor{revisionV3}{-- DONE in v3, superseded in scope: the v3
   open-economy calibration implements the government block, investment,
   exports, import margins, the saving rate and all three financings
   (see section 4, intermediate results). The explicit $\tilde\beta$
   tilt is ported to the canonical pipeline (\texttt{experiments/run.jl},
   \texttt{f1\_shift = tilt\_g0\_over\_c0}); the tilted matrix cells are
   pending.}
2. **BETA**: elastic total labour supply on the real wage,
   $L^s = \bar L\,[(w/P)/(w_0/P_0)]^{\eta_s}$, with the labour--leisure
   interpretation and the numerical elasticity test
   $\mathrm{d}\log L / \mathrm{d}\log(w/P) \approx \eta_s$.
   \textcolor{revisionV3}{-- Implemented ($\eta_s$, $\eta$-continuation
   solver, single-point elasticity identification verified at the v2
   stage: 0.5/1.0/2.0 recovered); the v3 verification run is pending
   ($\eta$-continuation runtime).}
3. **DELTA corner**: fixed real factor price, unconstrained factor
   quantity, Leontief limit of the CES core; document the exact
   equivalence conditions (R2's own point, answered).
   \textcolor{revisionV3}{-- DONE: implemented as :fixed + Leontief
   limit; the analytic equivalence is EXACT for F2 and F3 (rel y error
   0.0).}
4. **Pre-register the BF endpoints** before any Stage 2 run: only
   $\eta \in \{0,1\}$ are kept (ADR-0010); no interpolated $\eta^{*}$.
5. **Cobb-Douglas limit guard** (Milestone C): sign-safe real powers at
   $\varepsilon \to 1$.
   \textcolor{revisionV3}{-- Implemented (analytic branch at
   $\varepsilon = 1$); the continuous-epsilon assertion is pending its
   v3 re-run.}
6. **Full residual validation** (ROADMAP Phase 4):
   \textcolor{revisionV6}{all-N clearing
   (ADR-0019; the omitted-market rotation is retired),}
   household-expenditure
   exhaustion for every experiment type, homogeneity, multi-start
   convergence, Tornqvist consistency. Machine-precision residuals are
   the gate.
   \textcolor{revisionV3}{-- Partially done: budget identities exact at
   every solved equilibrium; the fixed-wage regime enforces all N
   clearings, the mobile regime enforces N-1 plus the CPI numeraire and
   asserts the omitted N-th market against S - (I+X-M) (ADR-0010); the
   omitted-equation rotation and multi-start battery are pending.}
   \textcolor{revisionV6}{-- v7 update: every regime now enforces all N
   clearings with the explicit external account (ADR-0019); the mobile
   eta = 1 identity S + $T_{\mathrm{int}}$ + M - (I+X) = F + $B_{\mathrm{gov}}$ and the all-N
   clearing residuals are gated in tests and experiments, with the eta = 0
   factor-market gap reported rather than gated. The equation-ordering
   rotation is retired with the omitted market; a multi-start battery on
   the full grid remains open.}

## Stage 2: Simulation and sensitivity \textcolor{revisionV2}{\normalsize [reworked v3]}

1. **Fill the 5 x 3 matrix**: headline set per cell (Tornqvist real
   GDP, employment, external balance), sectoral detail into the
   appendix; compare each cell against its pre-registered expected
   signature -- deviations are findings, not failures.
   \textcolor{revisionV4}{-- DONE (v5): all 15 cells executed on the A-bill
   calibration (section 4.2); the F1 column's fixed-wage cells became admissible
   once the heuristic guard was replaced (ADR-0014). The recombination option
   (allocation x supply elasticity) is recorded in section 4.2 as an available
   extension, deliberately left out because it is a null result for a
   demand-only programme.}
2. **Sectoral Sobol**: does the reallocation friction (BF) matter for
   sectoral allocation even where it is aggregate-second-order?
   (R1.7's aggregate-vs-sectoral demand.)
3. **Robustness**: shock magnitudes; the financing variants are already
   matrix columns; state the $B$-collapses-into-$T$ equivalence as a
   limitation.
4. **Elasticity table** (R1.3): values, ranges, sources, distributions,
   pre-registered.

\textcolor{revisionV2}{Exit criterion (ROADMAP Phase 5): headline
rankings qualitatively stable across financing closures and shock
sizes; otherwise ranking instability itself is the result.}

## Stage 3: Manuscript and reviewer response \textcolor{revisionV2}{\normalsize [reworked v3]}

\textcolor{revisionV2}{The manuscript is written from the two-wedge
design; every referee point is addressed at a named location. The
nine-section structure (definitive guide Part III, adapted):}

1. **Introduction** -- the two wedges; state upfront that the IO--CGE
   bridge is known (Robinson 2006, Rose 1995); one short Kuhnian
   paragraph at most (R2.3; R1.1).
2. **Related literature** -- the closure taxonomy; Robinson (2006)
   central; Rose (1995); Dervis--de Melo--Robinson (1982); de Melo \&
   Tarr (1992); McGregor--Swales--Yin (1996); Shoven \& Whalley (1984);
   Mansur \& Whalley (1984); Willenbockel (1994) (R2.2; R1.8).
3. **Data and accounting** -- the section 4.1 pipeline, import
   separation, decomposed value added, the documented 5.387\% residual,
   the calibration table (R1.3; R2.4).
4. **Model and closures** -- verbal model overview before equations
   (R1.2); static horizon stated; the two wedges as design dimensions;
   financing closures explicit (R2.4, R2.5).
5. **Aggregate results** -- matrix rows; wage-regime vs elasticity
   channels (R1.7).
6. **Sectoral results** -- sectoral Sobol, bottlenecks (R1.7).
7. **The financing wedge** -- bracket results, the external-debt
   variant, the bounded endogenous-money paragraph.
8. **Robustness and limitations** -- financing variants, shock sizes,
   one-factor aggregation, static horizon, open-economy treatment.
9. **Conclusion** -- practical closure guidance; no bridge-discovery
   claims (R2.1/R2.3 answered by measurement).

\textcolor{revisionV2}{Reviewer 1 point map: R1.1 structure/central
argument -> sections 1 and 3 (the matrix); R1.2 verbal model description
-> section 4; R1.3 calibration summary table -> section 3 and the
elasticity table; R1.4 Type-I/II multipliers -> clarify or drop in
section 2; R1.5 wage-setting rule -> GAMMA/BETA with explicit floors;
R1.6 mobility/unified wage -> BF/ALPHA; R1.7 aggregate vs sectoral ->
sections 5--6 with the sectoral Sobol; R1.8 McGregor--Swales--Yin ->
literature and the DELTA corner; R1.9 minor copy-edits (p. 5 sentence,
p. 16 paragraph, Figure 3 axis).}

\textcolor{revisionV2}{Reviewer 2 point map: R2.1 labour-supply
elasticity as bridging parameter -> BETA and the matrix; R2.2 closure
literature -> section 2; R2.3 cleavage framing -> domain-of-applicability
reframing, section 1; R2.4 closed-economy miracle -> section 3 pipeline
plus the external-debt closure; R2.5 sector-specific immobile
labour -> the mobile core (ALPHA/BF); R2.6 eq.~(19) sign error -> section 4
fix; R2.7 supply-vs-demand characterization -> removed; R2.8 CES novelty
claim -> dropped; R2.9 Cobb-Douglas archetype -> dropped; R2.10 arbitrary
uniform allocation -> eliminated by endogenous allocation (BF/ALPHA);
R2.11 "fundamentally contested" allocation claim -> dropped or
referenced; R2.12 internal mobility contradiction -> resolved by the
mobile core; R2.13 McGregor et al. -> section 2 and DELTA; R2.14 skill
classes -> optional extension with the limitation stated.}

\textcolor{revisionV2}{Methodological correction: the
corrected-equilibrium story (Result 1 of the definitive guide) is the
response letter's opening; the disclosure level remains a Stage 0
decision. Must not appear: all $\eta$-dominance claims, the
88.4\%/100\% figures, GO certification language, price-invariance
claims, and "the IO multiplier requires $\eta \to \infty$" (it
requires sticky wages).}

## Stage 4: Response letter and submission \textcolor{revisionV2}{\normalsize [reworked v3]}

\textcolor{revisionV2}{Reviewer 2: maximally gracious; concede "the
bridge exists" fully; present the corrected equilibrium and the measured
matrix as the contribution. Reviewer 1: point-by-point via the map
above; calibration table; McGregor engagement; aggregate-vs-sectoral
quantification. Then the final gates (ROADMAP section 11): clean-checkout
reproduction, reconciled accounts, explicit financing, identity tests,
sensitivity shares conditional on stated distributions, tag the
submission commit.}

## Immediate next step \textcolor{revisionV2}{\normalsize [reworked v3]}

\textcolor{revisionV2}{Stage 0 is resolved except the disclosure decision
and the scope locks. Begin Stage 1 with the financing core ($g_i$, the
$T$ rule, the external account, the explicit F1 experiment), then BETA,
then the DELTA corner; pre-register the BF endpoints
($\eta \in \{0,1\}$, ADR-0010) before any Stage 2 run.}

\textcolor{revisionV3}{Update (v4): the financing core and the DELTA
corner are done; the remaining Stage 1 items are the BETA verification
run, the ALPHA/GAMMA v3 re-runs ($\eta = 1$, CES elasticities), the F1
tilted-cell runs (tilt ported to the canonical pipeline; cells pending),
and the interpretive pass on the v3 baseline
units in notebook 04 -- then re-preregister the design
($\eta \in \{0,1\}$, ADR-0010) and open Stage 2.}

## Running the pipeline (container and Mac) \textcolor{revisionV3}{\normalsize [added v4, folded]}

\textcolor{revisionV3}{Execution venue. Notebooks 01 (data wrangling) and
02 (accounting consistency) are validated end-to-end and reproduce the
parent \texttt{AC\_*} artifacts exactly; notebook 03 (financing closures)
needs a v3 re-run once the solver route is settled. Notebooks are
executed in the container (Julia 1.12.7, depot
\texttt{/opt/julia-depot}, kernel \texttt{julia-1.12}); the Mac is a
second venue for the same code, useful for solver experiments, not for
speed.}

\textcolor{revisionV3}{Can runtime be ruled out? Yes, in the narrow
sense: a single 140-dimensional solve takes seconds and the failures are
convergence stalls at floors that are identical across algorithms and
run-times -- more CPU time does not fix them. A Mac run is nevertheless
NOT superfluous, but its purpose is solver robustness in a clean
environment: the session evidence is init/path-sensitive (the identical
system converged to $2\times10^{-10}$ under one init and stalled at
$3.6\times10^{-4}$ under another), so a controlled init sweep is the
cheapest decisive experiment.}

\textcolor{revisionV3}{Commands. Container (repo root
\texttt{BFRep/(3)BeyondHulten}): \texttt{julia -{}-project=. -e 'using
Pkg; Pkg.instantiate()'} once; acceptance test
\texttt{julia -{}-threads=4 cbase2/scripts/verify\_v3.jl} (the standing
gate: baseline continuation, F1/F2/F3 rows, budget identities, DELTA
equivalence). Mac (same repo checked out, \texttt{juliaup} 1.12.x):
first \texttt{julia -{}-project=. -e 'using Pkg; Pkg.instantiate(); using
BeyondHulten'}; then (1) the acceptance test as above, (2) an init sweep
on one calibrated economy -- default \texttt{[ones(N); lambda; 1.0]},
warm linear fixed point, and small random perturbations of each -- over
$\theta \in \{2.0, 1.0, 0.5\}$ with \texttt{drops=[71]}, recording
retcode and residual per run (a short script to be added as
\texttt{cbase2/scripts/solver\_sweep.jl}); (3) the same with
\texttt{AutoFiniteDiff()} and with Levenberg--Marquardt as primary. If
any route reaches machine-tolerance residuals, the pipeline proceeds
unchanged; if all stall, implement the IPOPT/JuMP formulation.}

\textcolor{revisionV2}{The chronological document map, the
per-document evaluation of documents 1--10, and the synthesis have been
moved to \texttt{docs/archive/document\_review\_leftovers.md} at
Version 3; this working document now contains the foundations, the
evaluation matrix, and the workplan only.}

# Revision Log \textcolor{revisionV1}{\normalsize [added v2]}

- **Version 1** (September 2026) --- Initial chronological assessment of `docs/`.
- **Version 2** \textcolor{revisionV1}{(September 2026)} --- Fixed a
  LaTeX build error (backticked paths with underscores inside
  \texttt{\textcolor} spans); added the chronological roadmap chapter
  (Stages 0--4) with Decision overview I (labour-market closures) and
  Decision overview II (financing closures) as Stage 0 tables with
  recommendations, folding in all still-relevant details from documents
  1--8; documents 1--8 moved to \texttt{docs/archive/}
  (\texttt{definitive\_guide.md} and \texttt{labor\_closures.md} remain
  active); version history block added; roadmap restructured into
  Foundation I (labour-market closure) and Foundation II (financial
  closures) with per-option two-column provenance tables (formulation,
  status, referee point, cost, recommended role, reviewer, old docs,
  Baqaee--Farhi usage), followed by the workplan; corrected the geometric
  reallocation row ("used by BF" -- the BF 2019 replication ships both
  the immobile and the mobile variant, verified in the original MATLAB
  files); added the exogenous labour-endowment shift as a ninth closure
  table (BF 2022 BLS labour-supply shifts; the original paper's slack
  injection); extended the coverage check with firm capacity constraints
  and adjustment-cost frictions, plus a to-verify nuance on the BF 2019
  mobile variant.; added a lineage note under the exogenous-shift table
  documenting the direct resonance with the original manuscript's
  Section 5 (slack narrative borrowed from BF 2022's social-distancing
  channel; code mapping of the three slack implementations to the
  manuscript's two variants); Foundation I regrouped into the mandatory
  minimum (reallocation margin, Closure A, sticky wage, Closure D) and
  the extended option set (Closure B, partial mobility, Closure C,
  legacy callback, exogenous shift); final selection implemented: 1.1
  selected closures BF (explicitly justified vis-a-vis BF 2019's
  neoclassical short-run market-friction argument) and ALPHA--DELTA
  (neoclassical standard, neoclassical labour--leisure,
  production-agnostic Keynesian, Leontief-inverse emulation) with a
  pinned suggested narrative; 1.2 notable unselected candidates (partial
  mobility, complementarity, exogenous shift as scenario device); 1.3
  dropped legacy callback; financing deduction routes clarified
  (household composition vs.\ government budget vs.\ lump-sum tax),
  static-horizon nuance added (debt collapses into the tax variant;
  foreign financing is the genuinely distinct robustness), recommendation
  updated to tax principal with preference-reallocation and foreign
  robustness; header numbering repaired (subsection level under the
  foundation sections), coverage check moved to Foundation I as its
  fourth subsection, plausibility parenthetical removed; financing
  design set per user decision: principal bracket of preference
  reallocation and tax financing, external-debt (foreign financing)
  third variant as the static analog of the endogenous-money intuition,
  imports absorbed mechanically at fixed sector import shares
  (Armington extension flagged via the existing Omega split); new
  section "The Evaluation Matrix (5 x 3)" added with pre-registered
  expected signatures per cell, the limit relations (BF endpoints, BETA
  to GAMMA, GAMMA plus Leontief = DELTA corner), and the sketch-derived
  figure "Paradigmatic directions towards the short run"
  (\texttt{docs/pictures/paradigmatic\_directions.png}); Foundation I
  merged the exogenous-shift table into the Dropped section (one closure,
  legacy callback = endowment shift); Foundation II restructured into
  Selected (bracket + external debt) and Dropped (expenditure switching;
  pure debt subsumed by the tax variant) with a final endogenous-money
  paragraph; compile-all.zsh repaired (pandoc now called with
  --resource-path so relative image paths resolve against the md
  directory).
- **Version 3** \textcolor{revisionV2}{(September 2026)} --- Stage 0
  answered by the status matrix (5 x 3 with symbols: what exists, what
  is missing per cell); the second paradigmatic wedge named (financing;
  loanable-funds vs endogenous-money) with a bounded-strength
  recommendation; Workplan reworked (Stages 0--4) with inline
  reviewer-response maps (R1.1--R1.9, R2.1--R2.14) and the
  two-wedge manuscript structure; chronological review sections moved
  to \texttt{docs/archive/document\_review\_leftovers.md}.
- **Version 4** \textcolor{revisionV3}{(September 2026)} --- Intermediate
  implementation results of the \texttt{cbase2/} pipeline reported: the
  Stage 0 chapter promoted to its own section (section 4) with an updated
  status matrix; the remaining stages form the Workplan section
  (section 5). New subsection "Intermediate results" records the three
  structural findings of the v3 open-economy recalibration (the
  fixed-wage indeterminacy without leakages and the degenerate E = 0
  corner; the national-accounts identity S = I + X - M forcing the
  saving rate s = 0.398, calibrated from the data together with tau0 =
  0.214, export share 0.422, investment share 0.162; the corrected
  equilibrium formulation with N-1 mobile clearings plus the CPI
  numeraire and the asserted residual external account, the price-indexed
  F1/F3 tax and the balanced-budget F2 rule) and the verified
  results (injection continuation 5.4 s without stalls, price-explosion
  branch eliminated, F1/F2/F3 mobile rows at resid <= 3.7e-7 with exact
  budget identities, F3 external balance recorded, DELTA analytic
  equivalence EXACT for F2 and F3). Stage 1 items annotated with their
  state (financing core and DELTA corner done; BETA implemented with the
  verification run pending; CD guard implemented; residual validation
  partially done); immediate next step updated (BETA verification,
  ALPHA/GAMMA v3 re-runs, F1 tilted-cell runs, notebook 04 interpretive
  pass, then pre-registration).

- **Version 5** \textcolor{revisionV4}{(September 2026)} --- The executed 5 x 3
  matrix and the records that made it runnable. New section 4.2 reports all 15
  cells of the \texttt{matrix\_5x3\_v3} design (headline set plus the findings
  against the pre-registered signatures) and records the recombination option
  (the allocation x supply-elasticity corner) as a deliberate omission. Section
  4.1's clamped-calibration numbers and its fixed-wage indeterminacy finding are
  corrected in place (ADR-0012/0013/0014); section 3's signature table gains the
  same pointer; Stage 2's matrix item is marked done. Accounting: the A-bill
  calibration (row 73 domestic bill, row 74 imports, row 75 product taxes as
  three booked leaks) makes the external-identity canary exact at 1.7e-16; the
  raw table's 5.387 \% production-vs-expenditure residual stays open. Solver: the
  cbase2 ladder is retired (DE-0010) and the polish target is 1e-10 (ADR-0015),
  so cell metrics no longer depend on where the solver stopped. One
  pre-existing unescaped percent sign (Stage 3, the documented production-vs-
  expenditure residual) is escaped so the line no longer truncates its tail; no
  other text changed.

- **Version 6** \textcolor{revisionV5}{(September 2026)} --- Real GDP and
  welfare separated in the measurement layer (ADR-0018): income-side real GDP
  w*sum(L) deflated by a seven-component Tornqvist GDP deflator, the
  expenditure Divisia index as its dual, and the external wedge
  Sigma(V) - w*sum(L) reported as a diagnostic. New kernel API
  (gdp\_components, gdp\_deflator, gdp\_income, gdp\_expenditure, gdp\_wedge,
  real\_consumption); manifest schema v2 records the gdp and consumption
  metric keys plus the seven component diagnostics. Because the measurement
  layer changed, the matrix is re-minted as matrix\_5x3\_v4: 15 of 15 cells
  executed, all gates pass, and the v4 equilibria reproduce v3 to 1.2e-11.
  Section 4.2 carries the v6 table and corrects the v5 F2 reading (mobile
  income GDP flat, -1.69 percent is welfare); the manifest-generated flow
  tables are in \texttt{paper/tables/matrix\_5x3\_v4\_flows.md}. The F2 open
  gate in \texttt{registry/closures.toml} is resolved by the same runs.

- **Version 7** \textcolor{revisionV6}{(September 2026)} --- The external
  account harmonized (ADR-0019): all N clearings enforced in every regime,
  the mobile eta = 1 system carrying the endogenous net external transfer F
  and F3 booking the programme financing as B\_gov, with the booked identity
  S + T\_int + M - (I+X) = F + B\_gov. The matrix is re-minted as
  matrix\_5x3\_v5 (15 of 15 cells executed; v1-v4 stay as history); section
  4.2 carries the v7 headline table with the booked external account and
  corrects the superseded readings (external deficit = F, the
  wedge-as-external-position, the F2/F3 separation); the manifest-generated
  flow tables are in \texttt{paper/tables/matrix\_5x3\_v5\_flows.md}. The
  generation's finding is financing neutrality (F2 and F3 real-neutral in
  the mobile regime; the v1-v4 F2/F3 differences were omitted-market
  shortcut artifacts). Measured movement against v4: mobile F1 welfare
  +0.043 percent to +0.143 percent, mobile F2/F3 welfare -1.694 percent to
  -1.823 percent, fixed and BF rows unchanged in their real allocation; the
  BF eta = 0 external-position entry moves by definition, from -0.7947
  percent (the v4 residual canary) to +1.3310 percent (the booked B\_gov),
  and is corrected in section 4.2 of this version as booked financing over
  an open account rather than an identified position (ADR-0020, accepted the
  same day with option C: the eta = 0 endpoint moves to sectoral wages, the
  account then closes and financing neutrality extends to eta = 0; the
  specification and the probe measurement are at the end of section 4.2).
  \textcolor{revisionV6}{Option C was promoted and executed the same day: the
  matrix is re-minted as \texttt{matrix\_5x3\_v6} (15 of 15 cells, commit
  \texttt{6db1da5}), in which the BF row carries an identified external
  position (-0.581 / -0.853 / -0.853 percent of GDP) with the account closed to
  1.0e-11 or better and financing neutrality extended to eta = 0, while the
  twelve non-BF cells reproduce v5 (ALPHA/BETA to 6.7e-16; GAMMA/DELTA to
  1.6e-11, the fixed-wage sensitivity to the warm start (the system is well conditioned: Jacobian cond about 6.3 measured on the GAMMA-F2 cell, so the amplification is in the measurement layer, ~100x the residual level), verified to be
  independent of the change). \texttt{matrix\_5x3\_v5} and the correction above
  stay as the record of the superseded pin.}
  \textcolor{revisionV6}{The same date added, within this version: the current
  GAMMA row as an explicit reference block (unchanged from the v7 headline
  table, since the promotion reproduced the twelve non-BF cells), and a note on
  variation in GAMMA --- the one-degree-of-freedom diagnosis, the measured wage
  level invariance and wage structure responses, the closed external account
  under a tilted pin, the conditioning result, the three limits that belong with
  the numbers, and the pointer to \texttt{docs/VariationinGamma.md} and
  ADR-0021 (proposed).}
  The F3 definition gates in
  \texttt{registry/closures.toml} are resolved by the same runs; the
  \texttt{revisionV6} colour used here assumes the render template defines
  it (as \texttt{revisionV5} was introduced for Version 6). Section 4.2 also
  classifies the sources of the recorded deviations (measurement definition,
  closure artifact, model-class property, corrected admissibility
  mathematics), and the ADR-0019 amendment of the same date records the
  harmonization approaches considered (status-quo canary, all-N plus transfer
  F, savings-driven investment, Armington expenditure switching, direct
  external-balance imposition).
