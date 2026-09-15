---
title: "Document Assessment: What Remains Relevant in docs/"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-09-14"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [assessment, docs-audit, revision, beyondhulten, metroeconomica]
last-updated: September 2026
---

**Version 2** \textcolor{revisionV1}{(September 2026)}
**Version 1** (September 2026)

This document evaluates every planning document in `docs/` in chronological
order and records, for each, whether anything in it remains relevant for the
remaining work process. The governing documents outside `docs/`
(`ROADMAP.md`, `roadmaps/vertdict.md`) are included as anchors because
several `docs/` files reference them as their arbiter.

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
reallocation parameter and two wage regimes. The nine option tables
below consolidate the full design space, grouped into the mandatory
minimum (endpoints and the BF-comparable reallocation margin) and the
extended option set. Provenance tags: "reviewer"
attributes R1/R2; "old docs" names documents (archived under
\texttt{docs/archive/}, plus the governing \texttt{ROADMAP.md} and
\texttt{roadmaps/vertdict.md}); "used by BF" points to the Baqaee--Farhi
papers (both source-verified in \texttt{bf\_replication/} and
\texttt{bf\_replication2/}):}

- BF 2019 (Econometrica): `Dokumente/Baqaee-2019-The Macroeconomic Impact of Microeconomic Shocks- Beyond Hulten's Theorem-Econometrica.pdf`
- BF 2022 (AER): `Dokumente/Baqaee-2022-Supply and Demand in Disaggregated Keynesian Economies with an Application to the COVID-19 Crisis-American Economic Review_4.pdf`

### Mandatory minimum: endpoints and the BF reallocation margin \textcolor{revisionV1}{\normalsize [added v2]}

\textcolor{revisionV1}{The group the paper cannot do without: the two
endpoints its question lives between (Closure A, Closure D), the
BF-comparable reallocation margin whose endpoints are exactly the two BF
2019 variants, and the sticky wage -- one step from the IO endpoint and
the carrier of the current headline result (the wage-regime table,
Result 3).}

| **Geometric reallocation parameter (current $\eta$)** | *Endogenous sectoral allocation between immobile and cost-minimizing, total labour fixed* |
|---|---|
| **Formulation** | $L_i = L_{i,\mathrm{fixed}}^{1-\eta}\, L_{i,\mathrm{costmin}}^{\eta}$ with $\eta \in [0,1]$ (extrapolation outside); total employment fixed at $\bar L$; **not** a labour-supply elasticity |
| **Status** | Implemented and tested (`MobileLaborCES`, `:mobile`); Sobol first-order share about 15.7% |
| **Referee point addressed** | The mobility margin only; total supply stays fixed |
| **Cost** | None |
| **Recommended role** | Keep as sub-analysis under an honest name ("reallocation parameter") |
| **Mentioned by reviewer** | both (R1.6 mobility/wage-rule question; R2.5 immobile-labour criticism -- indirectly) |
| **Mentioned in old docs** | `labor_closures.md`; `definitive_guide.md`; `WORKPLAN3.md` (2026-09-03 addendum) |
| **Used by BF** | Yes at its endpoints -- the original BF 2019 replication ships both variants: **immobile labour** ("no reallocation"; sector-specific wage vector, `Simulation_Derivs.m`) and **mobile labour** ("full reallocation"; single wage pinned at the numeraire, `Simulation_Derivs_realloc.m`; `Read_Me.txt`: "the first replicates results with immobile labor and the second replicates results with mobile labor"). The geometric interpolation parameter *between* the endpoints is project-specific |

| **Closure A: Full-employment mobile labour** | *The standard full-employment CGE endpoint* |
|---|---|
| **Formulation** | $\sum_i L_i = \bar L$; one economy-wide wage; cost-minimizing allocation; the $\eta = 1$ limit of the reallocation parameter (as `:mobile`) |
| **Status** | Implemented (`:mobile`); baseline endpoint of the wage-regime sweep |
| **Referee point addressed** | R2.5 (economy-wide mobile labour); R1.6 (unified wage rule enabling mobility) |
| **Cost** | None |
| **Recommended role** | Baseline endpoint of the closure continuum |
| **Mentioned by reviewer** | both (R2.5 explicitly recommends it; R1.6 asks for it) |
| **Mentioned in old docs** | `REVISED_SALVAGE_PLAN.md`; `WORKPLAN.md`; `WORKPLAN2.md`; `vertdict.md`; `definitive_guide.md` |
| **Used by BF** | Yes -- BF 2019 baseline: competitive core with inelastically supplied (fixed) composite-factor endowment and uniform factor prices (see filename above) |

| **Sticky wage, uncapped (current \texttt{:fixed})** | *Fixed real wage; employment absorbs the shock -- the extensive margin* |
|---|---|
| **Formulation** | $w/P = \bar\omega$ (wage pinned at baseline); employment endogenous and uncapped; the gap $\bar L - \sum_i L_i$ is computed post-solve, not an equilibrium constraint |
| **Status** | Implemented; large positive GDP responses (+19.3 pp at mult 10, employment to 1.19x), residuals at most $1.9 \times 10^{-7}$ |
| **Referee point addressed** | R1.5 (wage-setting rule: fixed real/nominal wage); R2's "fix the wage rate at baseline and drop the labour supply constraint" recipe |
| **Cost** | None |
| **Recommended role** | Report as-is, with the uncapped-employment caveat stated |
| **Mentioned by reviewer** | both (R2 gives the recipe explicitly; R1.5 demands the wage-setting rule) |
| **Mentioned in old docs** | `WORKPLAN3.md` (Milestone F + 2026-09-03 addendum); `definitive_guide.md` (Result 3); `labor_closures.md` |
| **Used by BF** | Yes (adapted) -- BF 2022's sticky-wage mechanism (downward wage rigidity, demand-determined employment) in simplified one-wage form (see filename above) |

| **Closure D: IO endpoint** | *Fixed real factor price, unconstrained factor quantity; the exact conditions under which the CGE reproduces the IO multiplier* |
|---|---|
| **Formulation** | Fixed real primary-factor price; unconstrained factor quantity within the experiment; Leontief input coefficients (or a demonstrated CES special case with unchanged relative prices) |
| **Status** | Not implemented as an explicit closure; `:fixed` approximates its spirit |
| **Referee point addressed** | R2 explicitly (Robinson 2006: fix primary input prices, unconstrained supplies -- the multiplier model); R1 explicitly (McGregor--Swales--Yin long-run zero-price-change IO-type result) |
| **Cost** | Small--moderate |
| **Recommended role** | Required: R2's entire critique turns on the exact conditions of equivalence |
| **Mentioned by reviewer** | both (R2 with the Robinson 2006 quote; R1 with McGregor et al. 1996) |
| **Mentioned in old docs** | `vertdict.md` (FixedRealWage endpoint); `ROADMAP.md` (Closure D); `definitive_guide.md` (Part IV); `REVISED_SALVAGE_PLAN.md` (the superseded $\eta \to \infty$ endpoint idea) |
| **Used by BF** | No (origin is the IO tradition and Robinson 2006, not the BF papers; BF 2019 computes efficient reallocations, not fixed-price multipliers) |

### Extended option set \textcolor{revisionV1}{\normalsize [added v2]}

\textcolor{revisionV1}{The recommended addition (Closure B, the
referee-facing supply elasticity), the structural partial-mobility
alternative on the allocation axis, the optional complementarity closure,
and the two lineage closures inherited from the original manuscript.}

| **Closure B: Elastic total labour supply (the referee's elasticity)** | *Voluntary labour-supply response along the real wage; the continuum from vertical to horizontal supply* |
|---|---|
| **Formulation** | $L^s = \bar L\, [(w/P)/(w_0/P_0)]^{\eta_s}$ with a labour--leisure interpretation; $\eta_s = 0$ fixed supply; $\eta_s \to \infty$ approaches the fixed-real-wage limit |
| **Status** | Not implemented; specified in ROADMAP sections 4.2 and 5 (Closure B) |
| **Referee point addressed** | R2's core demand ("the one obvious set of elasticities that really matters") |
| **Cost** | Moderate: extends the $2N{+}1$ system by the supply equation and its anchor |
| **Recommended role** | **Principal addition** (recommendation) |
| **Mentioned by reviewer** | both (R2 explicitly and centrally; R1 implicitly via the slack-scenario discussion) |
| **Mentioned in old docs** | `REVISED_SALVAGE_PLAN.md` (original `elastic_labor_slack` proposal); `WORKPLAN.md`; `WORKPLAN2.md`; `PRELIMINARY-ASSESSMENT.md` (baseline-anchoring correction); `vertdict.md` (ElasticLaborSupply); `ROADMAP.md` |
| **Used by BF** | No (2019: inelastic endowments; 2022: sticky wages -- vertdict names the 2019 base as the right foundation for adding it) |

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

| **Legacy exogenous labour-slack callback** | *The original paper's mechanism: exogenous sectoral slack vectors* |
|---|---|
| **Formulation** | Exogenous `labor_slack` functions (`full_labor_slack`, `full_labor_slack_alt`, `empirical_labor_slack`) in the base CES model; sector-specific slack vectors |
| **Status** | Implemented (legacy); the base of the original submission's Sections 5.1--5.2 |
| **Referee point addressed** | Itself the criticised mechanism |
| **Cost** | None |
| **Recommended role** | Retire to benchmark only |
| **Mentioned by reviewer** | both (as the criticised object: R2 "deus ex machina"; R1 "unsure about the rationale behind this scenario") |
| **Mentioned in old docs** | `REVISED_SALVAGE_PLAN.md` (inventory); `WORKPLAN.md`; `WORKPLAN2.md`; `labor_closures.md`; `definitive_guide.md` |
| **Used by BF** | No |

| **Exogenous labour-endowment shift (labour-slack injection)** | *Exogenous sectoral labour-supply changes, independent of wages* |
|---|---|
| **Formulation** | Sectoral labour endowments shifted exogenously: $\bar L_i \to \bar L_i (1+g_i)$ with calibrated or assumed $g_i$; no wage response required |
| **Status** | Implemented (the legacy callback is exactly this mechanism); the original submission's Sections 5.1--5.2; ROADMAP section 4.4 lists it as a distinct closure concept |
| **Referee point addressed** | R1.5 describes precisely this mechanism ("is it simply a matter of increasing labour supply by the number of unemployed persons?"); R2 criticizes its use as the central bridge |
| **Cost** | None |
| **Recommended role** | Keep available as a named, honestly-labelled scenario closure (shock design); not as the paper's central mechanism |
| **Mentioned by reviewer** | both (R1 describes the mechanism; R2 criticizes its use as deus ex machina) |
| **Mentioned in old docs** | `REVISED_SALVAGE_PLAN.md`; `WORKPLAN.md`; `WORKPLAN2.md`; `labor_closures.md`; `definitive_guide.md` |
| **Used by BF** | Yes -- BF 2022 calibrates exogenous sectoral labour-supply shifts to BLS May-2020 sectoral hours (in the COVID application they are **contractions**: lockdowns and social distancing "reduced the supply of labor"; \texttt{bf\_replication2/src/model.jl} applies them as \texttt{A\_shock} on the labour block). Run in reverse, the same closure is an endowment expansion -- the slack injection. BF 2019 shocks factor endowments in its general framework |

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
regime, and an explicit financing closure.}\textcolor{revisionV1}{Recommendation (unchanged): add Closure B as the
referee-facing principal specification; keep the implemented pair
(reallocation parameter, wage regime) as measured sub-analyses; Closure C
only if the policy-evaluation framing survives ROADMAP Phase 6; Closure D
in either case, because R2's critique turns on the exact equivalence
conditions. The sticky-wage row must always be reported with its
uncapped-employment caveat.}

# Foundation II: Financial Closures \textcolor{revisionV1}{\normalsize [added v2]}

\textcolor{revisionV1}{The four admissible experiment types of ROADMAP
Phase 2. The current \texttt{demand\_shock} is an unfinanced add-on and
inadmissible by the ROADMAP's own rule until one principal and one
robustness variant are chosen; the explicit demand anchor introduced with
the financing closure is also expected to resolve the $\eta = 1$ stall
(scale indeterminacy) observed from multiplier 0.5 upward. Neither BF
paper supplies a fiscal financing closure (verdict); only the
compositional-shifter normalization has a BF 2022 analogue.}

| **Budget-neutral preference reallocation** | *Compositional household demand shift with the budget closed by construction* |
|---|---|
| **Budget identity** | $\tilde\beta_i = \beta_i d_i \,/\, \sum_j \beta_j d_j$; $Z(p) = \sum_j \tilde\beta_j p_j^{1-\sigma}$; $c_i^h = E^h \tilde\beta_i p_i^{-\sigma} / Z(p)$, so $\sum_i p_i c_i^h = E^h$ holds exactly |
| **Economic interpretation** | Compositional demand shift; must **not** be described as an autonomous investment multiplier (verdict warning) |
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

| **Expenditure-switching investment** | *Reallocation within a fixed public budget* |
|---|---|
| **Budget identity** | $\sum_i p_i g_i = T$ with $T$ raised by cutting other public spending |
| **Economic interpretation** | Public budget reallocated across sectors; no new purchasing power |
| **Cost / status** | Low--medium |
| **Recommended role** | Robustness alternative |
| **Mentioned by reviewer** | no |
| **Mentioned in old docs** | `vertdict.md`; `ROADMAP.md` (Phase 2); `definitive_guide.md` |
| **Used by BF** | No |

| **Debt/foreign-financed net expenditure** | *Demand expansion with an intertemporal or external counterpart* |
|---|---|
| **Budget identity** | $\sum_i p_i g_i = B + F$ with a defined saving--investment or current-account rule |
| **Economic interpretation** | Net new expenditure against borrowing or the external balance |
| **Cost / status** | Medium--high: additional institutional balance |
| **Recommended role** | Robustness (or headline extension) |
| **Mentioned by reviewer** | no |
| **Mentioned in old docs** | `vertdict.md`; `ROADMAP.md` (Phase 2); `definitive_guide.md` |
| **Used by BF** | No |

\textcolor{revisionV1}{Recommendation (unchanged): tax-financed as
principal, preference reallocation as robustness. Tax financing keeps the
"investment programme" framing, is the easiest variant for a referee to
audit, and the preference-reallocation variant then doubles as the
compositional-shift benchmark that separates the two experiment families
cleanly.}

### Coverage check: BF closures outside our lists \textcolor{revisionV1}{\normalsize [added v2]}

\textcolor{revisionV1}{Verified against the two papers
(\texttt{Dokumente/}) and the replications. BF 2019 itself uses nothing
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

# Workplan \textcolor{revisionV1}{\normalsize [added v2]}

## Stage 0: Conceptual decisions (now -- discussion, no code)

\textcolor{revisionV1}{The first steps are conceptual; each decision
determines what is implemented and written later. None of them requires
running anything. The option tables with full provenance tags live in
Foundation I and Foundation II above; the recommendations stated there
are the proposed agenda.}

Remaining items not covered by the foundation tables:

1. **Disclosure of the internal audit trail.** Decide how much of the
   rejected-round history (verdict, 88.4% artifact, corrected Sobol) is
   disclosed in the response letter; see Stage 4.
2. **Scope decisions.** Skill-class disaggregation (R2 stretch goal,
   ROADMAP optional extension, becomes mandatory only if the one-factor
   model is presented as policy evaluation); open-economy depth beyond the
   section 4.1 documented mapping; whether B\&F replication (Stream C) stays
   off the critical path (all surviving plans agree it validates but does
   not block).

## Stage 1: Model completion (after Stage 0 decisions)

1. **Shared final-demand and institutional-budget core** (verdict redesign
   \S1, ROADMAP Phase 3): household demand
   $c_i^h = E^h \tilde{\beta}_i p_i^{-\sigma}/Z(p)$ with renormalized
   weights, government/investment budget $\sum_i p_i g_i = T + B + F$.
2. **Implement the chosen financing closure** plus the explicit demand
   anchor. This is expected to also fix the $\eta=1$ stall from multiplier
   0.5 upward (scale-indeterminacy guard without an autonomous/investment
   anchor) and the increasingly negative flexible-wage response pattern
   (interpretation depends on the financing choice).
3. **Implement the labor-closure decision from Stage 0.1** (rename, or
   partial-mobility, or Closure B). Keep the three dimensions of
   `labor_closures.md` separate: legacy callback, reallocation parameter,
   wage regime.
4. **Cobb-Douglas limit guard** (Milestone C, still open): sign-safe real
   powers so $\varepsilon \to 1$ solves without `DomainError`.
5. **Closure D (IO endpoint)** (ROADMAP section 5): fixed real factor
   price, unconstrained factor quantity, exact or qualified equivalence to
   the Leontief multiplier; document the equivalence conditions. This is
   the referee-facing "the bridge exists, and here is exactly when it
   appears" result.
6. **Full residual validation** (ROADMAP Phase 4): omitted-equation
   invariance (rotate the omitted commodity market), household-expenditure
   exhaustion, homogeneity, multi-start convergence, real-GDP index
   consistency. Machine-precision residuals are the gate.

## Stage 2: Simulation and sensitivity (gates in ROADMAP Phase 5)

1. **Re-run the canonical sweeps under the final specification**
   (`rerun_results.jl`): wage-regime table and Sobol decomposition must be
   regenerated after any Stage-1 change; the 2026-09-03 numbers are
   canonical only for the current code state.
2. **Sectoral Sobol** (definitive guide Part V.3, still pending): does the
   reallocation parameter matter for sectoral allocation even if not for
   aggregate GDP? R1 explicitly expects the aggregate-vs-sectoral contrast
   to be quantified (R1.7). Build it on the `SobolResult` conventions from
   `milestone_D_plan.md`: first-order and total-order indices, absolute
   shares, no renormalization, complete designs, no silently dropped cells.
3. **Financing and shock-magnitude robustness** (ROADMAP Phase 5): repeat
   headline results across the robustness financing closure and reasonable
   shock sizes; report ranking stability.
4. **Elasticity table** (R1.3, ROADMAP required table): values, ranges,
   sources, distributions for the chosen parameters, pre-registered in the
   repository before final results are inspected.

## Stage 3: Manuscript rewrite (ROADMAP section 7; definitive guide Part III)

\textcolor{revisionV1}{The surviving structural guidance is: nine-section
structure per \texttt{definitive\_guide.md} Part III, written into the
\texttt{revised\_manuscript/chapters/} scaffold (which still contains only
old text). Per-section details carried forward from the archived plans:}

1. **Abstract and introduction** from the corrected three-part story
   (methodological correction; mobility material-but-secondary; wage regime
   first-order). Retire all commensurability language; state upfront that
   the IO--CGE bridge is known (Robinson 2006, Rose 1995); keep one short
   Kuhnian paragraph at most. The `OUTLOOK_SUMMARY.md` abstract draft is
   quantitatively dead (built on the rejected 88.4\%) and must not be used
   even as a stylistic base.
2. **Section 2 (model exposition) fixes:** sign error in eq.~(19)
   (R2.6); remove the "supply-side vs.\ demand-side" characterization
   (R2.7); streamline the IO exposition; clarify the purpose of the
   Type-I/Type-II multiplier discussion or drop it (R1.4); state static
   horizon, closure rules, and calibration replication explicitly (R1.2).
3. **Data and accounting section** from `accounting_consistency_plan.md`:
   section 4.1 pipeline, import separation, decomposed value added, the
   documented 5.387\% expenditure residual, shock-incidence rule (domestic
   final demand), and the calibration/summary tables R1 asked for
   (intermediate shares, VA shares, Domar weights, elasticities).
4. **Results sections** with the corrected numbers only: wage-regime table
   (pp gaps, not ratios), Sobol shares ($\theta$-led, $\eta \approx
   15.7\%$ first-order), aggregate-vs-sectoral contrast; figures with
   clarified axes, units, and GDP measures (R1.9, R2 last bullet).
5. **Terminology discipline** per `labor_closures.md`: never conflate the
   legacy callback, the reallocation parameter, and the wage regime; state
   the `:fixed` closure's economic meaning exactly; drop the "deus ex
   machina" hand-wringing in favor of the positive result.
6. **Literature integration** (R2.2, R1.8): Robinson (2006) central, plus
   Rose (1995), Dervis--de Melo--Robinson (1982), de Melo \& Tarr (1992),
   McGregor--Swales--Yin (1996), Shoven \& Whalley (1984), Mansur \&
   Whalley (1984), Willenbockel (1994). Remove CES-novelty and
   Cobb-Douglas-archetype claims (R2.8, R2.9); use the
   domain-of-applicability framing instead of the Kuhnian cleavage (R2.3).
7. **Minor copy-edits** (R1.9): the p.~5 convoluted sentence, the
   incomplete p.~16 paragraph, the Figure~3 axis question.
8. **What must not appear** (definitive guide "remove" list): all
   $\eta$-dominance claims, 88.4\%/100\% figures, "GO" certification
   language, price-invariance claims, and the claim that the IO multiplier
   requires $\eta \to \infty$ (it requires sticky wages).

## Stage 4: Response letter and submission

\textcolor{revisionV1}{The archived plans contain a fully formed reviewer
response strategy; its notes are folded in here as subitems of the
drafting stage.}

1. **Reviewer 2 (devastating, recommends rejection):** maximally gracious;
   concede the core criticism ("the bridge already exists") fully and
   without defending the original framing; present the corrected
   equilibrium and the wage-regime/Sobol results as the genuinely new
   quantitative contribution; address every specific point (sign error,
   closed economy, immobile labor, CES novelty, labor allocation); engage
   the closure literature. Checked against the R2.1--R2.14 map in
   `REVISED_SALVAGE_PLAN.md` section 2 (archived).
2. **Reviewer 1 (constructive):** systematic point-by-point against
   R1.1--R1.9; show the restructured paper, the calibration summary table,
   the McGregor--Swales--Yin engagement, and the quantified
   aggregate-vs-sectoral distinction.
3. **Methodological-correction section:** decide (Stage 0.4) how much of
   the internal audit trail to disclose; if disclosed,
   `PRELIMINARY-ASSESSMENT.md` and WORKPLAN3 Milestone~A (archived)
   supply the technical content, and the supersession-note style of
   WORKPLAN2 is the disclosure template.
4. **Final gates** (ROADMAP section 11, definition of done): every
   quantitative claim reproducible from a clean checkout; accounts
   reconcile; financing closure explicit; all identity tests pass;
   sensitivity shares conditional on stated distributions; tag the exact
   submission commit; manuscript, response letter, figures, and code all
   agree with one version.

## Immediate next step

\textcolor{revisionV1}{The proposed agenda for the next conceptual
discussion: the Foundation I recommendation (Closure B as referee-facing
principal) and the Foundation II recommendation (tax-financed principal,
preference reallocation as robustness), plus the two remaining conceptual
items in Stage 0.}

# Chronological document map

| # | Document | Date | Status in one line |
|---|----------|------|--------------------|
| 1 | `REVISED_SALVAGE_PLAN.md` | 2026-07-29 | Foundational strategy; superseded in its results, still the best review-point map |
| 2 | `WORKPLAN.md` | 2026-07-29 | Superseded; historical record only |
| 3 | `WORKPLAN2.md` | 2026-07-29 (upd. later) | Superseded; its "GO / 88.4%" round was formally rejected |
| 4 | `OUTLOOK_SUMMARY.md` | 2026-08-02 | Superseded; abstract draft is built on invalid numbers |
| 5 | `PRELIMINARY-ASSESSMENT.md` | 2026-08-18 | The audit that killed the false claims; bug documentation still valuable |
| 6 | `milestone_D_plan.md` | 2026-08-19 | Executed and complete; historical |
| 7 | `WORKPLAN3.md` | 2026-08-18 (addendum 2026-09-03) | Repair record; its 2026-09-03 addendum holds canonical Part I numbers |
| 8 | `accounting_consistency_plan.md` | 2026-09-02 | Executed and integrated; now reference documentation |
| 9 | `definitive_guide.md` | 2026-09-03 | **Current authoritative narrative**; the manuscript should be written from it |
| 10 | `labor_closures.md` | 2026-09-05 | Current code semantics; contains an open design question that needs a decision |

External anchors:

- `ROADMAP.md` (2026-08-18/20): the governing plan, validation gates, and
  definition of done. Not superseded by anything in `docs/`.
- `roadmaps/vertdict.md` (2026-08-18): the independent audit verdict that
  rejected the first mobile-labor round. It is the arbiter WORKPLAN2 defers to.
- `docs/reviews/metro-rev1.docx` and `metro-rev2.docx`: the referee reports
  (R1 = constructive major revision; R2 = devastating reject). Extracted text
  is available at `docs/reviews/metro-rev1.txt` / `metro-rev2.txt` if needed.

\textcolor{revisionV1}{Note (v2): documents 1--8 of the map now live in \texttt{docs/archive/}; only \texttt{definitive\_guide.md}, \texttt{labor\_closures.md}, this assessment, and \texttt{reviews/} remain in \texttt{docs/} directly.}

# Evaluation, oldest first

## 1. REVISED_SALVAGE_PLAN.md (2026-07-29)

The foundational pivot document. It contains three things of unequal value.

**Still relevant:**

- **Section 2, the review-point map** (R1.1-R1.9, R2.1-R2.14) is the most
  complete mapping of referee demands to revision actions anywhere in the
  repository. Every response letter will be checked against it. Keep.
- The strategy insight of section 3.1 ("closure dominance claim") survives in
  modified form: acknowledge that the IO-CGE bridge is known (Robinson 2006,
  Rose 1995) and contribute measurement, not discovery. The corrected model
  qualifies the finding (see doc 9), but the *strategic posture* -- graceful
  concession to R2 plus a quantitative contribution -- is unchanged.

**Superseded / must not resurface:**

- The claim that "the core infrastructure exists and works" (sections 1.1
  and 7) was written before the mobile-labor audit and is wrong for the
  files it praises (`mobile_labor.jl`, `variance_decomposition.jl` at that
  date).
- The proposed labor-supply function `L = L_bar * (w/w_bar)^eta` (eta as
  labor-supply elasticity, eta in [0, inf)) is **not** what was finally
  implemented. Current code uses eta as a geometric intersectoral
  reallocation parameter with eta in [0, 1] (see doc 10). Any sentence in
  the manuscript drafted from this plan must be re-checked against
  `labor_closures.md`.
- The 6-10 week timeline and the phase plan are dead; the corrected critical
  path is in `ROADMAP.md` section 10.

**Dispositions for remaining work:** extract the review-point map into the
response-letter skeleton; otherwise treat as archive.

## 2. WORKPLAN.md (2026-07-29)

The three-stream coordination document (A: manuscript, B: model extension,
C: B&F replication).

**Still relevant:** only the *topology* -- three parallel streams with the
go/no-go gate between model work and writing -- and the Stream C pointer to
`bf_replication/REPLICATION_WORKPLAN.md`. The stream structure is inherited
by every later plan, so the file has organizational value as the origin of
the scheme.

**Superseded:** all statuses, all timelines, the priority matrix, and the
paths (it references `(1)Submission/revised/`, which now lives at
`revised_manuscript/` inside this repo). Its per-task efforts fed WORKPLAN2's
false "complete and verified" certification, which `vertdict.md` rejected.

**Disposition:** archive. Nothing unique except the stream topology.

## 3. WORKPLAN2.md (2026-07-29, updated)

The post-pilot workplan. Its headline claims -- "complete and verified", GO
decision, price invariance, eta = 88.4% -- were **formally superseded** by
the supersession note at the top of the document itself, which cites
`roadmaps/vertdict.md`. The document is internally annotated: it now reads as
a historical record of the rejected round.

**Still relevant:**

- The supersession note itself is a model of how the project handles
  retracted results; it should be cited in the response letter's
  methodological-correction section if the authors choose to disclose the
  internal audit trail.
- The Julia environment alignment record (1.12.6, container/host portability)
  remains operationally true.
- The destructuring-bug fix note (named-tuple destructuring of structs) is a
  real Julia pitfall worth remembering for future code review.

**Superseded:** every checkmark/GO/88.4% marker, the "Stream B results are
ready, Sections 5-6 can proceed immediately" instruction, and the
eta-as-labor-supply semantics throughout.

**Disposition:** archive, clearly marked as the rejected round. Do not use
for the manuscript.

## 4. OUTLOOK_SUMMARY.md (2026-08-02)

The response strategy essay ("From Searching for Bridges to Measuring the
Bridge") plus a draft abstract.

**Still relevant:**

- The rhetorical frame of the pivot is sound and still governs: retire the
  commensurability language, concede R2's point, contribute measurement.
- The priority matrix's pending Stream A tasks (introduction, Section 2
  rewrite, literature integration, figures, response letter) are an accurate
  list of what is *still* pending today -- Stream A has not advanced since.

**Superseded:**

- The draft abstract is quantitatively dead: it reports the eta = 88.4% and
  "order of magnitude" claims and the eta-as-labor-supply-elasticity reading
  of `mobile_labor.jl`. It must not be used even as a stylistic base without
  replacing every number and the eta interpretation.
- The claim "Done. The `mobile_labor.jl` module implements
  `L = L_bar * w^eta`" describes the since-rejected specification.

**Disposition:** keep the pivot rhetoric in mind; discard the abstract draft.
The abstract will be written anew from `definitive_guide.md` Part III.

## 5. PRELIMINARY-ASSESSMENT.md (2026-08-18)

The independent audit that caught the false certification. Verified by
execution that: the eta-sweep response was inverted (eta = 10 collapsed GDP
to 12%), equilibria were initialization-dependent, the variance decomposition
renormalized main effects, Cobb-Douglas grid points threw `DomainError`, and
the ROADMAP-vs-workplans contradiction existed.

**Still relevant -- genuinely, not just historically:**

- Its bug catalog is the reference description of *why* the old numbers were
  wrong. If the response letter includes a methodological-correction section
  (recommended by `definitive_guide.md` Result 1), this document supplies the
  technical content.
- The verification table mapping `vertdict.md` code-line claims to source
  lines is a completed audit trail.
- Its corrections section (overspend was 4.5% in the 71-sector model, not
  25.5%; the decisive fix is baseline-wage anchoring, not nominal-vs-real)
  prevents repeating two specific misreadings.

**Superseded:** nothing. It was the superseding document. The deleted
`tests/minimal_test/` diagnostics are noted; the organized test suite now
covers the same ground.

**Disposition:** keep as the audit record. Source for the response letter.

## 6. milestone_D_plan.md (2026-08-19)

A short task plan for the honest Sobol decomposition. Every step (first-order
and total-order indices, absolute shares, CSV output, robust missing-point
handling) was implemented and verified (WORKPLAN3 Milestone D plus the
2026-09-03 addendum).

**Still relevant:** only as the specification-of-record for what
`summary_table` and `SobolResult` are supposed to output. If the sectoral
Sobol extension (still pending, see doc 9 Part V) is built, this plan is the
template.

**Disposition:** archive; consult before extending the sensitivity code.

## 7. WORKPLAN3.md (2026-08-18, addendum 2026-09-03)

The repair workplan (Milestones A-F) with the critical 2026-09-03 re-run
addendum. This is where the *current canonical Part I numbers* live.

**Still relevant -- load-bearing:**

- **The 2026-09-03 wage-regime table** (flexible eta=0 turning negative with
  shock size; sticky `:fixed` positive and large, +19.3 pp at mult 10;
  eta=1 stalls from mult 0.5): this is Result 3 of the corrected story and
  the single most important empirical table for the revised paper.
- **The 2026-09-03 Sobol result** (theta = 0.395 dominant, sigma = 0.273,
  epsilon = 0.165, eta = 0.157 material-not-dominant): this is Result 2 and
  the numbers the manuscript will report -- *not* the 88.4%, *not* the "~2%
  negligible" reading, both of which are documented artifacts.
- The accounting canon: GDP P = I = 3,027,818; E = 2,864,724; residual
  5.387% documented valuation gap; `sum(lambda) = 2.1099`.
- Milestone A's list of equilibrium repairs (sector-1 zero-profit restored,
  normalized household demand, retcode checking, Tornqvist index) is the
  methodological-correction content again, from the fix side.

**Superseded:** the historical Milestone D (eta around 0.01%) and Milestone F
(74-112x ratios) numbers are explicitly superseded by the addendum; the
addendum's own supersession note is unambiguous.

**Disposition:** keep as the numerical record. The manuscript's empirical
section is built from this file plus `definitive_guide.md`.

## 8. accounting_consistency_plan.md (2026-09-02)

The section-4.1 accounting transformation (ROADMAP Phase 1): separate
imports, decompose value added, reconcile the three GDP sides, define shock
incidence, emit calibration artifacts. Executed, integrated into
`src/interface.jl`, and validated by the test suite.

**Still relevant:**

- It is the reference for the data pipeline: the verified source schema
  (rows 73-83 of the Destatis table), the proportional import-allocation
  assumption (stated explicitly), the 5.387% expenditure residual (a genuine
  raw-table valuation gap, not a bug -- after the off-by-one indexing fix),
  and the field mapping in the integration note.
- The manuscript's "Data and accounting" section (structure item 3 in both
  ROADMAP section 7 and definitive guide Part III) will be written almost
  directly from this document. It answers R2's "closed-economy miracle"
  demand with a documented, open-economy-consistent pipeline.
- The calibration-table demand of R1 is answered here (Step 6 artifacts,
  `output/AC_*.csv`, regenerable).

**Superseded:** nothing material. The ~0.8% residual figure was corrected in
place to 5.387%.

**Disposition:** keep as the data-pipeline reference.

## 9. definitive_guide.md (2026-09-03)

The current authoritative synthesis. It supersedes WORKPLAN/WORKPLAN2,
operationalizes the ROADMAP's Phase 6 go/no-go with evidence, and defines the
corrected three-part story:

1. **Methodological:** the original model was not an equilibrium; a
   corrected, verified baseline exists (machine-precision residuals).
2. **Negative but precise:** eta (intersectoral mobility) is a material but
   secondary channel (first-order Sobol share around 0.157, well below
   theta = 0.395); the mobility bridge is second-order (at most 0.07 pp);
   eta=1 solves only at multipliers 0.1/0.2.
3. **Positive:** the wage regime (sticky vs flexible) is the first-order
   labour-market margin (pp gaps +0.32 to +24.75); the question is "is the
   economy at full employment or in an unemployment regime?", not "how
   mobile is labour?".

**Still relevant -- this is the master document:**

- Part III's manuscript structure (9 sections), headline abstract results,
  the "remove/add" lists, and the ROADMAP compliance table.
- Part V's remaining-work list, which is the best available backlog:
  - Sobol on **sectoral** quantities (not yet run);
  - ROADMAP Phase 2 policy experiment (financing closure) -- **not started**;
  - Milestone C solver stability at the Cobb-Douglas limit -- open;
  - manuscript rewrite -- not started;
  - response letter -- not started;
  - `:fixed` with an explicit demand anchor (to show the sticky-wage result
    is not eta-driven) -- open.
- Part IV records that Closure D (IO endpoint) and Closure C
  (unemployment complementarity) are unimplemented.

**Superseded:** nothing; it internally supersedes the 2026-09-02 re-run
numbers with the 2026-09-03 ones.

**Disposition:** active. All manuscript writing starts here. Its remaining
work items should become the working task list.

## 10. labor_closures.md (2026-09-05)

The latest document, and the one that defines what the code *actually does*
now.

**Still relevant -- critical for correctness:**

- The three-dimension taxonomy: legacy exogenous `labor_slack` callback;
  geometric eta reallocation (`L_fixed^(1-eta) * L_costmin^eta`, eta in
  [0,1], beyond which is extrapolation); wage regime `:mobile`/`:fixed`.
  These must never be conflated in the manuscript -- and note that **the
  current eta is NOT the labor-supply elasticity the referees and the
  salvage plan talked about.** The paper's terminology must be chosen
  deliberately: either rename the model's parameter (e.g. "reallocation
  parameter") or implement the standard partial-mobility / elastic-supply
  formulation and map the referee's eta onto it.
- The proposed **standard partial-mobility formulation** (sectoral wages,
  wage-responsive sectoral supply, 3N system) is an *open design proposal*,
  not implemented. It is the main unresolved modeling decision: the current
  geometric-eta plus efficiency-penalty design is a project-specific reduced
  form whose allocative-wedge curvature is not established as the exact CES
  allocative-loss coefficient. This directly touches R2's demand for
  labor-supply elasticities ("the one obvious set of elasticities that
  really matters").
- The warning that `:fixed` is not a capped unemployment model (employment
  can exceed `labor_bar`; the gap is computed post-solve) -- essential for
  honest interpretation of the +19.3 pp sticky-wage result. The sticky-wage
  closure is currently a *fixed real wage with unconstrained employment*,
  i.e. closer to Closure D's spirit than to Closure C. The manuscript must
  not sell it as a calibrated unemployment closure without saying so.

**Disposition:** active. The eta-semantics decision and the
sticky-wage-interpretation caveat feed directly into Sections 4 and 6 of the
new manuscript.

# Synthesis: where the project stands

## Verified and usable today

- A corrected, equilibrium-consistent 71-sector model (Milestones A, B and
  the 2026-09-03 fixes), with organized tests in `tests/`.
- Section-4.1 accounting-consistent data pipeline integrated in
  `src/interface.jl`.
- Canonical Part I results: wage-regime table, Sobol decomposition
  (theta-led, eta material at about 15.7%), accounting reconciliation -- all
  from the 2026-09-03 re-run, regenerable via
  `julia --project=. rerun_results.jl`.

## Open items (candidates for the working plan)

1. **eta-semantics decision** (doc 10): keep geometric eta with honest
   renaming, or implement the standard partial-mobility / elastic-supply
   closure that R2 actually asked for. This blocks the framing of Sections
   4-6.
2. **Sticky-wage closure interpretation** (doc 10): `:fixed` = fixed real
   wage with endogenous employment, not a calibrated unemployment regime;
   decide whether Closure C (complementarity) is implemented or the
   limitation is stated.
3. **Financing closure** (ROADMAP Phase 2; definitive guide Part V.5): the
   demand shock is still unfinanced -- an explicit closure is required
   before CGE results are admissible by the ROADMAP's own rule. This is the
   largest *unstarted* modeling item.
4. **Sectoral Sobol** (definitive guide Part V.3): does eta matter for
   sectoral allocation even if not for aggregate GDP? R1 expects the
   aggregate/sectoral contrast to be quantified.
5. **Cobb-Douglas limit** (Milestone C): epsilon = 0.99 DomainError guard.
6. **eta=1 stall under autonomous demand**: needs the explicit
   demand/investment anchor or an honest limitation note.
7. **Manuscript rewrite** per definitive guide Part III / ROADMAP section 7
   -- `revised_manuscript/` chapters still contain only the old text
   ("% OLD TEXT FROM HERE ON").
8. **Response letter** (R1 + R2), to be checked against the R1.1-R1.9 /
   R2.1-R2.14 map in `REVISED_SALVAGE_PLAN.md` section 2.
9. Minor copy-edits from R1 (p.5 sentence, p.16 paragraph, Figure 3 axis).

## Suggested next decision

Item 1 (eta semantics) is the fork: it determines whether the paper says
"the reallocation parameter eta is material but secondary" (what the code
now supports) or "the labor-supply elasticity continuum behaves as follows"
(what R2's framework expects and what would require the standard
partial-mobility implementation from doc 10). Everything downstream --
figures, Sobol labeling, abstract -- inherits from this choice.

\textcolor{revisionV1}{(v2) These decisions are now laid out with full option tables in Foundation I and Foundation II at the top of this document; the workplan stages supersede the open-items list here.}

## Revision Log \textcolor{revisionV1}{\normalsize [added v2]}

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
  legacy callback, exogenous shift).
