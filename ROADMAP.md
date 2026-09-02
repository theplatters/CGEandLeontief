# Roadmap: Labour-Market Closure and Green-Investment Multipliers

## 1. Objective

Produce a substantially reframed paper that asks:

> How strongly do labour-market closures, relative to production and demand elasticities, determine the aggregate and sectoral effects of a sector-specific green-investment programme?

The paper will use the German housing transformation as a 71-sector application. It will treat the relationship between input-output (IO) multiplier models and computable general equilibrium (CGE) models as established knowledge, not as a newly discovered “bridge.” The empirical contribution will be to quantify where results move between the full-employment CGE and IO-type endpoints and to distinguish aggregate sensitivity from sectoral sensitivity.

### Working title

**Labour-Market Closure and the Estimated Output Effects of Green Investment: Evidence from a 71-Sector German Model**

### Candidate central result

The result to test—not assume—is:

- labour-market closure principally determines the feasible aggregate employment and real-GDP response; and
- substitution parameters principally affect sectoral allocation, relative prices, displacement, and bottlenecks conditional on that aggregate closure.

The manuscript must not report a variance share such as “88.4% explained” until the underlying model, parameter design, and decomposition have passed the validation gates below.

## 2. Scope and non-goals

### In scope

- A static, comparative-equilibrium analysis with its time horizon stated explicitly.
- A common German IO calibration used consistently across IO and CGE specifications.
- Intersectorally mobile labour and an economy-wide wage in the main CGE variants.
- A transparent continuum from fixed labour to fixed-real-wage/unconstrained labour, supplemented by an unemployment-oriented closure.
- Explicit fiscal, external, and factor-market closures.
- Aggregate and sectoral sensitivity analysis with interactions.
- Reproducible scripts, tests, tables, figures, and a response-to-reviewers document.

### Out of scope for the core paper

- Claiming that CES functions or alternative factor-market closures are novel.
- A general philosophical argument about Kuhnian incommensurability.
- A definitive forecast of the German transition path.
- A full dynamic CGE model unless the static model fails the policy-relevance gate in Phase 6.
- Baqaee–Farhi productivity-shock replications as a primary contribution. They may validate parts of the solver in an appendix, but do not validate the housing demand shock or its financing closure.

## 3. Current verified baseline

- [x] Correct the orientation and conditional normalization of the German intermediate-input matrix `Ω` to the Baqaee–Farhi user-by-supplier convention.
- [x] Correct the shock-input indexing and regenerate results based on the corrected IO calibration.
- [x] Replace the old fixed-base consumption calculation with a Törnqvist/Divisia real-GDP quantity index.
- [x] Add tests for the real-GDP index and run the Julia test suite.
- [x] Regenerate the existing GDP-dependent figures with the corrected measure.
- [ ] Obtain and review the co-author's claimed `mobile_labor.jl`, `variance_decomposition.jl`, and related replication code. These files are not present in the current repository and their reported results are therefore provisional.
- [x] Audit the supplied `mobile_labor.jl` and `variance_decomposition.jl` once obtained. **Audit outcome: rejected.** The two-sector diagnostic shows the reported mobile-labor system is not an equilibrium (shock-weighted household expenditure ≈ 1.25485 > 1; Walras' law cannot recover the omitted goods-market condition; a zero-profit equation was dropped instead of a redundant market equation; the solver accepted convergence without checking residuals or its own return code; labor supply used nominal rather than real wages; the high-η limit was misinterpreted as unlimited labor instead of L → 0; and the mobile model used a fixed-base sum instead of the shared Törnqvist index). All current mobile-labor sweep results—including the apparent price invariance and the 88.4% claim—are invalid until regenerated. See `roadmaps/vertdict.md`.

The current implementation in `src/` is the reproducible baseline. Work reported in external summaries is not part of the evidentiary record until it has been added, reviewed, and tested here.

## 4. Non-negotiable model requirements

These requirements are acceptance conditions for the revised quantitative analysis.

### 4.1 Accounting consistency

- Reproduce the benchmark IO accounts without a shock.
- Document the treatment of imports, exports, product taxes, production taxes, subsidies, labour compensation, capital income, and mixed income.
- Do not relabel total value added as labour without an explicit and prominently stated aggregation assumption.
- Reconcile production-, expenditure-, and income-side GDP in the benchmark.
- State whether the shock applies to domestically produced goods, imports, or both.
- Publish a calibration table containing at least intermediate-input shares, value-added shares, final-demand shares, import shares, Domar weights, and the elasticity values/ranges.

### 4.2 Real variables and the numeraire

- Specify labour supply as a function of the real wage, for example

  \[
  L^s = \bar L\left(\frac{w}{P}\right)^\eta,
  \]

  rather than as a function of the nominal wage.

- Preserve every sector's zero-profit/unit-cost condition.
- If a price index is used as numeraire, remove a redundant market-clearing equation justified by Walras' law; do not replace a zero-profit equation merely to balance the equation count.
- Verify homogeneity: changing the nominal unit must not alter real allocations.
- Use the shared Törnqvist/Divisia helper for real GDP in every new model variant.
- Report real GDP separately from welfare. Do not call a welfare index “GDP” without proving their equivalence under the maintained assumptions.

### 4.3 Financing the investment shock

The additional final demand must have an explicit counterparty. Implement and document at least one principal financing closure and one robustness closure, for example:

- lump-sum or specified tax financing;
- substitution away from other government expenditure;
- debt financing with a defined saving–investment balance; or
- foreign financing with a current-account rule.

An unfinanced increase in final expenditure is not an admissible CGE experiment. A budget-neutral preference reallocation must not be described as an autonomous investment multiplier.

### 4.4 Labour-market interpretation

Keep the following concepts distinct:

- voluntary labour supply along a labour–leisure margin;
- involuntary unemployment or surplus labour;
- intersectoral labour mobility;
- perfectly elastic labour supply at a fixed real wage; and
- exogenous expansion of the labour endowment.

Each closure must have an economic interpretation, equations, calibration, and an appropriate time horizon.

## 5. Model variants

Implement a common production and final-demand core with interchangeable labour-market closures.

### Closure A: Full employment

- Fixed aggregate labour endowment.
- Labour mobile across sectors.
- One economy-wide market-clearing wage.
- Purpose: standard full-employment CGE endpoint.

### Closure B: Finitely elastic labour supply

- `L = Lbar * (w / P)^η` or a fully specified labour–leisure problem.
- Use an empirically defensible range for `η`.
- Purpose: trace adjustment between fixed labour and a high-elasticity labour response.
- Interpretation: voluntary labour-supply response, not involuntary unemployment.

### Closure C: Unemployment or wage-curve closure

- Specify either a fixed-real-wage/surplus-labour regime up to a capacity bound or a documented wage curve relating real wages to unemployment.
- Include complementarity or regime-switching conditions if full employment can become binding.
- Purpose: represent demand-induced employment when initially unused labour exists.

### Closure D: IO-type endpoint

- Fixed real primary-factor price.
- Unconstrained primary-factor quantity within the relevant experiment.
- Fixed input coefficients, or a demonstrated special case in which unchanged relative prices make the CES coefficients reproduce the benchmark IO coefficients.
- Purpose: establish the exact conditions under which the CGE implementation reproduces the IO multiplier.

### Optional extension: heterogeneous labour and capital

After the one-factor implementation passes validation, assess a model with:

- at least skilled and unskilled labour;
- different closure rules by skill class; and
- fixed short-run sectoral capital or another explicit capacity constraint.

This extension becomes mandatory if the one-factor model is presented as a policy evaluation rather than as a methodological illustration.

## 6. Implementation plan

### Phase 0 — Recover and audit proposed work

- [ ] Obtain the co-author's implementation and its exact commit history or patch.
- [ ] Compare it with the current `main` branch rather than copying files blindly.
- [ ] Reproduce every reported headline statistic from a clean environment.
- [ ] Audit equation counts, omitted equations, normalization, convergence criteria, and residuals.
- [ ] Retain only modules that satisfy repository conventions and the model requirements above.

**Exit criterion:** every proposed result is linked to version-controlled code and a deterministic command.

### Phase 1 — Reconcile the German accounts

- [ ] Construct a benchmark-accounting report from the 2019 German IO table.
- [ ] Separate domestic and imported intermediate and final uses where the source table permits it.
- [ ] Represent taxes and subsidies consistently with producer and purchaser prices.
- [ ] Define the treatment of exports and the external balance.
- [ ] Verify benchmark GDP from all available accounting approaches.
- [ ] Generate the manuscript calibration table automatically.
- [ ] Document any aggregation or closed-economy mapping that remains necessary.

**Exit criterion:** benchmark replication errors are below stated numerical tolerances, and no empirical table item changes economic meaning silently.

### Phase 2 — Define the policy experiment

- [ ] Fix the sectoral investment vector, units, base year, and source.
- [ ] Separate domestic purchases from imported purchases.
- [ ] Choose the principal financing closure.
- [ ] Define at least one alternative financing closure for robustness.
- [ ] State whether the experiment represents construction-phase demand, operating-phase effects, or both.
- [ ] Confirm that the IO and CGE models receive economically comparable shocks.
- [ ] **Distinguish the four admissible experiment types explicitly:**
  - **preference reallocation** — a budget-neutral household compositional shift; renormalize CES weights as `β̃ᵢ = βᵢ dᵢ / Σⱼ βⱼ dⱼ` so `Σᵢ pᵢ cᵢʰ = Eʰ`, and do **not** describe it as an autonomous investment multiplier;
  - **tax-financed public investment** — an autonomous `gᵢ` financed by an explicit tax rule `Σᵢ pᵢ gᵢ = T`;
  - **expenditure-switching investment** — a reallocation of existing public expenditure across sectors; and
  - **debt/foreign-financed net expenditure** — financed by `B` (borrowing) or `F` (foreign) with a defined saving–investment or current-account rule.
  A single unnormalized `demand_shock` vector must not serve two economically different experiments. The financing closure (principal + robustness) must be decided here, before any labor closures are implemented.

**Exit criterion:** summing all institutional budgets and commodity accounts leaves no unexplained source of purchasing power.

### Phase 3 — Implement mobile-labour closures

- [ ] Add model types and options without duplicating the common CES core unnecessarily.
- [ ] **First build the shared final-demand and institutional-budget core** (household demand `cᵢʰ = Eʰ β̃ᵢ pᵢ⁻ᵃ / Z(p)`, government/investment budget `Σᵢ pᵢ gᵢ = T + B + F`, and all institutional balances) before implementing any labor closure. This core is what makes the four Phase-2 experiment types separately executable.
- [ ] Implement Closure A and show that it reproduces the fixed-total-labour benchmark.
- [ ] Implement Closure B using the real wage: `L = L̄[(w/P)/(w₀/P₀)]ᵝ`.
- [ ] Implement Closure C with an explicit wage/unemployment rule, using a complementarity formulation `0 ≤ L̄ − L ;⊥; w/P − ω̄ ≥ 0` rather than approximating the fixed-real-wage endpoint with `η = 10⁶`. The 2022 Baqaee–Farhi complementarity machinery may be adapted, but its N-sector sticky-labor structure must not be copied unchanged.
- [ ] Implement Closure D and document its exact IO-equivalence conditions.
- [ ] Use warm starts only as a numerical optimization; results must not depend on path or starting value.
- [ ] Expose residual diagnostics and solver termination status in analysis scripts.

Suggested source layout:

```text
src/
  labor_closures.jl
  mobile_labor.jl
  sensitivity.jl
scripts/
  build_calibration_table.jl
  run_closure_sweep.jl
  run_global_sensitivity.jl
  generate_paper_outputs.jl
tests/
  test_accounting.jl
  test_labor_closures.jl
  test_endpoint_equivalence.jl
  test_sensitivity.jl
```

**Exit criterion:** all four closures solve the no-shock benchmark and the policy shock with residuals below tolerance.

### Phase 4 — Validate the economics and numerics

Add automated tests for:

- [ ] no-shock benchmark replication for every closure and elasticity configuration;
- [ ] row normalization and orientation of `Ω`;
- [ ] zero profit in every sector, including any sector associated with a dropped market equation;
- [ ] goods-market clearing and aggregate Walras residual;
- [ ] labour-market clearing;
- [ ] government, household, saving–investment, and external balances;
- [ ] homogeneity with respect to the nominal numeraire;
- [ ] the numerical elasticity `d log L / d log(w/P) ≈ η`;
- [ ] convergence from multiple initial conditions;
- [ ] continuity as `η` and CES elasticities vary;
- [ ] `η = 0` equivalence to fixed aggregate labour;
- [ ] convergence of the high-`η` model to the fixed-real-wage closure;
- [ ] exact or qualified equivalence of Closure D and the IO model;
- [ ] consistency of real-GDP, nominal-GDP, employment, and price indices.
- [ ] **omitted-equation invariance** — repeat the diagnostic while omitting each possible commodity market in turn; real allocations must be invariant to that arbitrary choice, and every omitted goods-market residual must be small.
- [ ] **household expenditure exhaustion** — verify `Σᵢ pᵢ cᵢʰ = Eʰ` (or the appropriate institutional total) exactly, for every experiment type.

**Red-flag diagnostic:** if normalized prices and the real wage do not change across `η` but employment does, stop and identify the violated equation before using the results. (Note: the earlier reported "price invariance with changing employment" was an artifact of the rejected mobile-labor system, not a valid CGE finding.)

**Exit criterion:** `julia --project=. -e 'using BeyondHulten'` and `julia --project=. -e 'using Pkg; Pkg.test()'` pass from a clean environment, and an independent accounting/residual report passes.

### Phase 5 — Run sensitivity analysis

- [ ] Define literature-based or empirically estimated ranges/distributions for `η`, `ϵ`, `θ`, and `σ`.
- [ ] **Document the assumed product probability measure** for the factorial design; require complete designs for factorial (first-order and total-effect) indices. For a complete factorial product distribution, compute first-order indices exactly as `S_f = Var(E[Y|f]) / Var(Y)` under the specified weights; with equal weights, describe it as a "first-order sensitivity index under a discrete uniform distribution over the selected levels."
- [ ] Use a bounded transformation such as `η / (1 + η)` for plotting the continuum, but do not treat a uniform distribution over that transformation as an empirical prior without justification.
- [ ] Pre-register the main parameter grid or sampling distributions in the repository before inspecting final results.
- [ ] Run a transparent factorial design for diagnostic response surfaces.
- [ ] Run a variance-based global sensitivity analysis with first-order and total-effect indices; **report interactions rather than forcing all variation into four main-effect shares**, and **do not renormalize the percentage shares over main effects** (the earlier "88.4%" was η's share of the *sum of reported main effects*, not 85.9% of total variance under the chosen discrete distribution).
- [ ] Extend the current `SobolResult` sensitivity result (with deprecated `VarianceDecompositionResult` alias) with explicit factor levels and probability weights; selected interaction indices; solver status and residuals for every evaluation; and the assumed parameter probability measure.
- [ ] Do not silently discard solver failures: retry with independent initial values and continuation; if any required factorial cell remains invalid, abort the decomposition; report the failed region as part of the feasible parameter domain. Use an unbalanced regression/ANOVA only as a separately named analysis with an explicitly chosen sum-of-squares convention.
- [ ] Add essential sensitivity tests: a purely additive analytical function; a pure interaction function; nonuniform factor weights; a deliberately missing factorial cell (must raise an error); a solver result with a failed return code but finite `.u`; and invariance to grid traversal and warm-start order.
- [ ] Repeat the analysis for alternative financing closures and reasonable shock magnitudes.
- [ ] Assess solver failures as part of the feasible parameter domain, not as observations to discard silently.

Primary outcomes:

- Törnqvist real GDP;
- aggregate employment and real wage;
- CPI and sectoral prices;
- sectoral gross output, value added, and labour demand;
- displacement/crowding out outside the directly stimulated sectors;
- distance from the IO outcome; and
- sectoral rank correlation and weighted dispersion across model variants.

**Exit criterion:** headline rankings remain qualitatively stable across reasonable parameter distributions, financing closures, and shock sizes. Otherwise, the paper's result becomes that ranking instability itself is economically important.

### Phase 6 — Decide the final model scope

Apply the following go/no-go decisions after the validated sensitivity results exist:

1. **Does the labour closure affect aggregate outcomes materially?**
   - If no, abandon “closure dominance” and report the actual result.
2. **Do CES elasticities materially affect sectoral incidence?**
   - If yes, make the aggregate/sectoral contrast the central contribution.
   - If no, narrow the paper to the identification and calibration of labour-market closures.
3. **Does the one-factor model yield credible policy interpretation?**
   - If no, add labour skill classes and fixed short-run capital or present the application solely as a methodological illustration.
4. **Does the IO endpoint emerge under transparent limiting assumptions?**
   - If no, state precisely which other assumptions prevent equivalence and make that decomposition part of the result.

## 7. Manuscript plan

Rewrite rather than incrementally edit the current narrative.

### Proposed structure

1. **Introduction**
   - Establish that IO is a limiting/special case under known closures.
   - State the quantitative question and the aggregate-versus-sectoral hypothesis.
   - Summarize results only after validation.
2. **Related literature**
   - IO–CGE relationship and multiplier interpretation.
   - Factor-market closures, labour mobility, unemployment, and wage curves.
   - Systematic/global sensitivity analysis in CGE models.
   - German green-investment and housing-transition context.
3. **Common accounting framework and data**
   - German 2019 table, imports, taxes, calibration, and shock construction.
4. **Model and closure variants**
   - Shared CES structure, labour closures, financing closure, equilibrium, numeraire, and real-GDP measurement.
5. **Endpoint comparison**
   - Conditions for fixed-labour CGE and IO-type outcomes.
6. **Aggregate results**
   - GDP, employment, real wages, prices, and sensitivity decomposition.
7. **Sectoral results**
   - Reallocation, bottlenecks, relative prices, rankings, and sensitivity decomposition.
8. **Robustness and limitations**
   - Financing, parameter ranges, shock magnitude, factor aggregation, static horizon, and open-economy treatment.
9. **Conclusion**
   - Practical guidance on selecting and communicating closures; no paradigm-bridge claims.

### Required tables and figures

- [ ] Model-closure comparison table with endogenous/exogenous variables.
- [ ] German calibration and national-accounting summary table.
- [ ] Elasticity values, ranges, sources, and distributions.
- [ ] GDP/employment response across the labour-closure continuum.
- [ ] IO-distance or multiplier-attainment curve.
- [ ] Aggregate sensitivity indices with interactions.
- [ ] Sectoral sensitivity/distribution figure.
- [ ] Robustness table across financing closures.

### Reviewer-specific revisions

- [ ] Give a verbal model overview before equations.
- [ ] State clearly that the model is static and define its adjustment horizon.
- [ ] Explain benchmark calibration and exact replication.
- [ ] Remove the Type-I/Type-II multiplier discussion unless it directly supports the experiment.
- [ ] Remove the false supply-side-versus-demand-side characterization of general equilibrium.
- [ ] Correct the sign in the cost equation.
- [ ] Remove claims that Cobb–Douglas is the canonical CGE form or that CES use is recent/novel.
- [ ] Use economy-wide mobile labour in the main specification.
- [ ] Incorporate Robinson, Rose, Dervis–de Melo–Robinson, McGregor–Swales–Yin, Shoven–Whalley, Mansur–Whalley, and later labour-closure comparisons.
- [ ] Clarify every figure axis, unit, base value, and GDP measure.

## 8. Reproducibility and repository discipline

- Keep source data in `data/` and document download/source metadata without committing restricted files.
- Move final analysis out of notebooks into deterministic scripts.
- Generate manuscript tables and figures from a single top-level script.
- Record package versions in `Manifest.toml`.
- Seed randomized sensitivity designs.
- Save parameter designs and numerical summaries in machine-readable form when licensing permits.
- Keep economic model changes, tests, generated outputs, and manuscript rewrites in focused commits.
- Do not commit notebook checkpoints or machine-local paths.
- Tag the exact commit used for submission.

## 9. Suggested commit sequence

1. `Document revised research strategy and validation gates`
2. `Add German IO accounting reconciliation`
3. `Define financed housing investment shock`
4. `Implement mobile full-employment labour closure`
5. `Add elastic and unemployment labour closures`
6. `Validate IO endpoint and equilibrium identities`
7. `Add global sensitivity analysis`
8. `Generate revised calibration tables and figures`
9. `Rewrite manuscript around labour-market closures`
10. `Add response to reviewers and replication instructions`

Generated figures should normally be committed separately from the model change that produced them, in accordance with the repository guidelines.

## 10. Indicative schedule

| Stage                                | Indicative duration | Main output                                      |
| ------------------------------------ | ------------------: | ------------------------------------------------ |
| Recover proposed code and audit      |            2–3 days | Reproducible baseline or rejection report        |
| Accounting and shock closure         |           1–2 weeks | Reconciled calibration and financed experiment   |
| Labour closures and validation       |           1–2 weeks | Tested model variants and endpoint results       |
| Sensitivity analysis                 |              1 week | Response surfaces and global sensitivity indices |
| Scope decision and extensions        |            3–7 days | Final one-factor or heterogeneous-factor design  |
| Manuscript rewrite                   |             2 weeks | Complete revised manuscript                      |
| Response letter and replication pass |              1 week | Submission package                               |

The critical path is accounting consistency → definition and financing of the policy shock → shared equilibrium core → labor closures → full residual validation (including omitted-equation invariance and expenditure exhaustion) → sensitivity analysis → writing. A replication of unrelated productivity shocks is not on the critical path.

## 11. Definition of done

The paper is ready for submission only when:

- [ ] every quantitative claim can be reproduced from a clean checkout;
- [ ] the benchmark German accounts reconcile;
- [ ] the investment shock has an explicit financing and external closure;
- [ ] every labour closure is expressed in real terms and has a clear economic interpretation;
- [ ] all unit-cost, market-clearing, institutional-budget, and numeraire tests pass;
- [ ] the IO endpoint is reproduced or its failure is explained analytically;
- [ ] real GDP uses the common documented quantity index;
- [ ] sensitivity shares are conditional on stated parameter distributions and include interactions, computed from a documented complete product probability measure with no renormalization over main effects;
- [ ] residual validation includes omitted-equation invariance and household-expenditure exhaustion for every experiment type;
- [ ] aggregate and sectoral conclusions are reported separately;
- [ ] the manuscript no longer claims to discover the IO–CGE bridge or the use of CES functions;
- [ ] the limitations of a static, aggregated-factor model are explicit;
- [ ] figures, tables, manuscript, response letter, and replication instructions agree with the same code version; and
- [ ] package loading and the complete test suite pass.
