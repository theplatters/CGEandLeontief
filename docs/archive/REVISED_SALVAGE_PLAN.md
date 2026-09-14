---
title: "Revised Salvage Plan — Based on Complete Code Assessment"
project: Kapeller & Scharnreitner, submitted to Metroeconomica
date: July 2026
tags: [salvage-plan, io-modeling, cge, labor-supply, beyond-hulten, revised]
---

# Revised Salvage Plan — Based on Complete Code Assessment

## Executive Summary

This plan revises and supersedes the four earlier plan documents (`SALVAGE_PLAN.md`, `SALVAGE_PLAN Kopie.md`, `salvage_plan2.md`, `salvage1_comment.md`). Those plans were written based on an **incomplete assessment of the code base** — they examined only `(dep)BFReplicate` (a simple notebook with a non-converging solver, US data, and disconnected German data) and concluded the code was "~40% functional." In reality, the main code base `(3)BeyondHulten` is a **full Julia package that is ~70–80% functional**, with working solvers, connected German data, and labor slack already implemented as a function parameter.

This changes the salvage timeline from the 12–16 weeks estimated in earlier plans to **6–10 weeks**, because the core infrastructure — the model, the data pipeline, the solver, the labor slack framework, the elasticity gradient, and the plotting — already exists and works.

---

## 1. What the Earlier Plans Got Wrong About the Code

### 1.1 The wrong code base was assessed

The earlier plans all reference `BFReplicate` — a Julia notebook project (`MasterReplication.ipynb` + `src/Analysis.ipynb`) that:
- Uses `NLsolve` (which fails to converge)
- Contains US data (`us80dbasedata.mat`) with German data disconnected
- Has stale MATLAB references (`mvnrnd.jl`, `fmincon`, `knitromatlab`)
- The `problem()` and `Jacobian()` functions exist but the solver falls back to finite differences

**The actual code base is `(3)BeyondHulten`** — a proper Julia package (`BeyondHulten.jl`) that has:

| Component | Status | Details |
|-----------|--------|---------|
| Module structure | ✅ Complete | `BeyondHulten.jl` with 8 source files, exports, tests |
| CES model | ✅ Working | `solve(model::Model{CES})` using `NonlinearSolve.jl`, converges at `reltol=1e-8` |
| Leontief model | ✅ Working | `solve(model::Model{Leontief})` — standard IO multiplier |
| Cobb-Douglas model | ✅ Working | `solve(model::Model{CobbDouglas})` |
| German IO data | ✅ Connected | 71 sectors, 2019, `I-O_DE2019_formatiert.csv` |
| Hornykewycz impulses | ✅ Connected | `impulses.csv` loaded, `impulse_shock()` implemented |
| Labor slack | ✅ Implemented | `labor_slack::Union{Function, Symbol}` parameter in CES struct; three implementations: `full_labor_slack`, `full_labor_slack_alt`, `empirical_labor_slack` |
| Elasticity gradient | ✅ Implemented | `elasticity_gradient()` sweeps θ, ε, σ across range(0.99, 0.015, 1000) with threading |
| Labor slack gradient | ✅ Implemented | Linear interpolation between full/no slack in `CompareModels.ipynb` |
| Progressive sector shocks | ✅ Implemented | Each sector shocked individually, GDP recorded |
| Plotting | ✅ Working | Makie-based bar charts, panel plots, error bars, impulse responses |
| Plots generated | ✅ Done | 10 PNG files in `plots/` covering all scenarios |
| B&F US replication | ✅ Exists | `Notebooks/Translation.ipynb` (Julia), `Translation_Python.ipynb` |
| `main.jl` pipeline | ✅ Working | Full pipeline: load data → set shocks → solve → plot for all configurations |
| Demand shock mechanism | ✅ Working | `calculate_investment!()` and `impulse_shock()` |

### 1.2 What is genuinely still missing

1. **Continuous labor supply elasticity (η continuum)** — The code has *discrete* labor slack modes (full slack ≈ IO case, no slack = CGE, empirical = 3.1% unemployment). What's needed is a continuous parameter η ∈ [0, ∞) governing labor supply as a function of the real wage. This is the core of Reviewer 2's criticism.

2. **Open economy** — The model is closed economy. German IO table has imports, exports, and taxes, but these are shoehorned in (Rev2 calls this "brutal force").

3. **Mobile labor with economy-wide wage** — Labor is still sector-specific (each sector has its own wage `w[i]`). Reviewer 2 specifically recommends intersectorally mobile labor with a uniform wage.

4. **Variance decomposition** — The factorial design η × (ε, θ, σ) with partial R² or Sobol indices, as proposed in `salvage1_comment.md`.

5. **Some notebook API mismatches** — Notebooks like `DemandShocks.ipynb` use an older API (`real_gdp(sol, relative=true)`, `Shocks(supply, demand)` with 2 args) that doesn't match the current source code.

### 1.3 Impact on timeline

| Earlier plan estimate | Revised estimate | Reason |
|----------------------|-------------------|--------|
| 12–16 weeks (Strategy C) | **6–10 weeks** | No Python rewrite, no model rebuild, no data pipeline build |
| 2–3 weeks model modification | **1–2 weeks** | Labor slack infrastructure exists; just add η parameter |
| 2–3 weeks simulation | **1–2 weeks** | Elasticity gradient infrastructure exists; extend for η sweep |
| 3–4 weeks B&F replication | **1 week** (optional) | Translation notebook already exists |
| 2–3 weeks writing | **2–3 weeks** | Still needed — substantial rewrite of framing |

---

## 2. The Reviews: What Must Be Addressed

### Reviewer 1 (moderate, constructive — major revision)

Key demands:
- R1.1: Structure unclear, central argument hard to follow → **Major restructuring**
- R1.2: CGE model structure not described (static/dynamic, closure rules) → **Add clear model description**
- R1.3: Calibration process undocumented → **Add calibration documentation + summary table**
- R1.4: Purpose of "Keynesian" Section 2.1 unclear (Type I vs Type II multipliers?) → **Clarify or streamline**
- R1.5: Wage-setting rule missing for labor slack → **Resolved by endogenous labor supply**
- R1.6: Labor mobility assumption unclear → **Switch to mobile labor, state clearly**
- R1.7: Aggregate vs sectoral results not differentiated enough → **Make this the core contribution**
- R1.8: Engage with McGregor, Swales & Yin (1996) → **Add to literature**
- R1.9: Minor: p.5 convoluted sentence, p.16 incomplete, Figure 3 axis → **Copy-edit**

### Reviewer 2 (devastating — reject)

Core criticism: **The bridge between IO and CGE already exists.** The IO model is a special case of a CGE with fixed prices and unconstrained factor supply (Robinson 2006). The key parameter is the labor supply elasticity (η = 0 → CGE, η → ∞ → IO). The paper completely misses this and instead uses exogenous labor slack injection ("deus ex machina").

Specific points:
- R2.1: Labor supply elasticity not identified as key bridging parameter → **Make η the central variable**
- R2.2: Existing CGE literature on factor market closures not engaged (Robinson 2006, Rose 1995, Dervis et al. 1982, de Melo & Tarr 1992) → **Literature review + integration**
- R2.3: "Paradigmatic cleavage" framing contrived; IO is a special case of CGE → **Reframe as "domain of applicability"**
- R2.4: Closed-economy model applied to open economy Germany → **Open the economy or be transparent**
- R2.5: Labor treated as sector-specific immobile → **Switch to mobile labor, economy-wide wage**
- R2.6: Sign error in eq. 19 (− should be +) → **Trivial fix**
- R2.7: "Supply-side view" vs "demand-side view" nonsensical in GE → **Rewrite**
- R2.8: CES functions claimed as novel (standard since 1970s) → **Drop novelty claim**
- R2.9: Cobb-Douglas as "archetypical neoclassical" is arbitrary → **Drop or nuance**
- R2.10: Arbitrary uniform allocation of labor slack (Sec 5.2) → **Eliminated by endogenous labor supply**
- R2.11: "Allocation of labor remains fundamentally contested" unsubstantiated → **Add refs or drop**
- R2.12: Internal contradiction re factor mobility (p.14) → **Resolved by mobile factor**
- R2.13: Recommends McGregor, Swales & Yin (1996) → **Engage with this**
- R2.14: Suggests splitting labor by skill class → **Stretch goal**

---

## 3. The Salvage Strategy

### 3.1 Core idea: Closure dominance claim

**New framing:** "Given that the IO model is a special case of a CGE model with perfectly elastic labor supply (Robinson 2006, Rose 1995), we systematically *quantify* the relative importance of factor market closures versus production technology assumptions for both aggregate and sectoral outcomes. Using the German housing transformation as a policy-relevant application, we show that the labor supply elasticity explains an order of magnitude more variation in sectoral outcomes than CES substitution elasticities."

This does four things:
1. **Acknowledges** Rev2's core criticism (the bridge exists)
2. **Quantifies** how much the closure choice matters — this is genuinely new
3. **Preserves** the sectoral analysis as the novel empirical contribution
4. **Keeps** the pluralist angle as "domain of applicability" rather than "paradigmatic cleavage"

### 3.2 What stays, what goes, what changes

| Element | Action | Rationale |
|---------|--------|-----------|
| German housing transformation application | **Keep** | Well-motivated, policy-relevant |
| Sectoral-level analysis (prices, quantities, welfare) | **Keep & expand** | Core contribution |
| IO/Leontief model exposition | **Keep, streamline** | Benchmark |
| CES production function | **Keep** | Drop novelty claim |
| Exogenous labor slack injection (Secs 5.1, 5.2) | **Replace** | Endogenous labor supply with η |
| "Paradigmatic cleavage" framing | **Soften** → "domain of applicability" | Rev2 |
| Kuhnian introduction | **Shorten** to one paragraph | Both reviewers |
| Discovery narrative | **Remove** | State known results upfront |
| Cobb-Douglas as "archetypical neoclassical" | **Drop** | Rev2 |
| Closed-economy model | **Fix or justify** | Rev2 |
| Sector-specific immobile labor | **Replace** with mobile labor | Rev2 |
| Labor supply elasticity η (new) | **Add as central variable** | Rev2 core |
| Variance decomposition (new) | **Add** | Novel quantitative contribution |

---

## 4. Work Plan

### Phase 1: Literature & Reframing (Weeks 1–2)

**Goal:** Absorb criticism, engage missing literature, reframe contribution.

**Tasks:**

1. **Read & integrate key references:**
   - Robinson (2006) — *the* core reference; must be central
   - Rose (1995) — IO-CGE relationship, 30 years ago
   - Dervis, de Melo & Robinson (1982) — standard CGE labor market closures
   - de Melo & Tarr (1992) — labor-leisure choice microfoundation for elastic labor supply
   - McGregor, Swales & Yin (1996) — long-run IO interpretation with unemployment (both reviewers)
   - Shoven & Whalley (1984), Mansur & Whalley (1984) — CES usage in early CGE
   - Willenbockel (1994) — CGE as response to IO limitations

2. **Rewrite introduction:**
   - Open with the policy question (how large are the economic effects of green investment?)
   - State upfront: IO is a special case of CGE with fixed prices and unconstrained supply (Robinson 2006)
   - Frame contribution: "While the aggregate relationship is well understood, the *sectoral* implications and the *quantitative dominance* of the closure choice over production technology parameters have not been systematically measured"
   - Shorten Kuhnian preamble to a single paragraph
   - New subtitle: something like "Labor supply elasticities and sectoral dynamics in the IO–CGE continuum"

3. **Rewrite Section 2 (model exposition):**
   - Streamline IO exposition
   - Clearly state: (a) static model, (b) closure rules, (c) calibration procedure
   - Fix sign error in eq. 19 (− → +)
   - Remove "supply-side vs demand-side" distinction (nonsensical in GE)
   - Be transparent about open-economy treatment

**Deliverable:** New introduction draft + revised Section 2

### Phase 2: Code Extension (Weeks 2–4)

**Goal:** Add the continuous η parameter and fix technical issues. **Extend existing code — do not rewrite.**

**Tasks:**

1. **Implement continuous labor supply elasticity η:**
   The existing `labor_slack` function parameter is the hook. Currently:
   ```julia
   CES(elasticities::CESElasticities, labor_slack::Union{Function, Symbol}, labor_reallocation::Bool)
   ```
   The `labor_slack` function receives the model and returns a labor vector. Currently three discrete modes exist. Add a new mode:
   ```julia
   function elastic_labor_slack(model::Model, η::Float64)
       # η = 0: perfectly inelastic (current CGE, fixed supply)
       # η → ∞: perfectly elastic (IO multiplier model)
       # 0 < η < ∞: intermediate
       # L^s = L_bar * (w/w_bar)^η where w_bar is baseline wage
   end
   ```
   This requires:
   - Computing the equilibrium wage endogenously (currently wages are sector-specific `w[i]`)
   - Making labor mobile with an economy-wide wage
   - The labor supply function responds to the aggregate real wage

   **Key design decision:** To implement mobile labor + economy-wide wage, the `problem()` function in `ces.jl` needs modification:
   - Replace sector-specific `w[i]` with a single `w` determined by aggregate labor market clearing
   - Add a labor supply equation: `L^s = L_bar * (w/P)^η` where P is the CPI
   - The labor allocation across sectors is then determined by the model, not exogenously

2. **Switch to mobile labor with economy-wide wage:**
   - Replace sector-specific `l_max_k` constraints with aggregate `L_max` (or no constraint when η > 0)
   - Impose uniform economy-wide wage `w`
   - This eliminates the entire "sectoral allocation" problem

3. **Open the economy (or be transparent):**
   - Option A (preferred): Add imports, exports, distinguish domestic vs imported intermediates in the CGE model
   - Option B (fallback): Be fully transparent about how the German IO table is mapped into the closed-economy framework
   - The IO multiplier analysis (Leontief model) already handles the open economy correctly

4. **Implement variance decomposition:**
   - Factorial design: η ∈ {0, 0.1, 0.5, 1, 2, 5, 10, 100, ∞} × (ε, θ, σ) ∈ {low, mid, high}
   - For each combination: solve model, record aggregate GDP, sectoral prices, sectoral quantities
   - Decompose variance: how much is explained by η vs. (ε, θ, σ)?
   - Report partial R² or Sobol indices
   - The existing `elasticity_gradient()` function provides the template — extend it to a 2D sweep

5. **Fix calibration documentation:**
   - Verify that the CGE model replicates the IO table in baseline (it should — the code calibrates from the data)
   - Report α, β, λ, Domar weights in a summary table
   - Include elasticity values

**Deliverable:** Working code with η parameter, mobile labor, variance decomposition. **Go/no-go decision** at end of Phase 2: do the sectoral results show interesting variation along the η continuum?

### Phase 3: Simulation & Analysis (Weeks 4–6)

**Goal:** Generate the new results.

**Tasks:**

1. **Elasticity sweep (η continuum):**
   - For η ∈ {0, 0.1, 0.5, 1, 2, 5, 10, 100, ∞}:
     - Impose housing transformation demand shock
     - Solve for new equilibrium
     - Record: aggregate GDP, sectoral prices, sectoral quantities, sectoral welfare, Domar weight changes
   - Compare with IO multiplier benchmark and Cobb-Douglas benchmark

2. **Aggregate results (expected — state upfront):**
   - GDP increases monotonically with η
   - At η = 0: GDP ≈ 0 (or slightly negative)
   - At η → ∞: GDP ≈ IO multiplier result
   - State this as a known result (Robinson 2006), don't "discover" it

3. **Sectoral results (the novel contribution):**
   - Which sectors are most/least sensitive to the closure choice?
   - Price dynamics: at low η, strong price effects; at high η, prices stabilize
   - Are there threshold effects or non-monotonicities at the sectoral level?
   - The existing plots (panel, impulse, diff_lambda) can be re-cast as points along the continuum

4. **Variance decomposition:**
   - η × (ε, θ, σ) factorial
   - How much variation in GDP, sectoral prices, sectoral quantities is explained by η vs. CES elasticities?
   - This is the "size claim": closure choice dominates technology choice by factor X

5. **Skill-class disaggregation (stretch goal):**
   - Split labor into unskilled (unconstrained) and skilled (constrained)
   - Following Rev2's suggestion
   - More realistic for Germany

**Deliverable:** Complete simulation results, figures, and the variance decomposition table.

### Phase 4: Writing & Restructuring (Weeks 6–8)

**Goal:** Rewrite the paper around the new results.

**New structure:**

```
1. Introduction
   - Policy question (green investment impacts)
   - State upfront: IO is a special case of CGE (Robinson 2006)
   - Contribution: quantify how much closure choice dominates technology choice
   - Short Kuhnian paragraph (not pages)

2. Model frameworks
   2.1 IO multiplier model (streamlined)
   2.2 CGE model with CES functions (clearly specified, closure rules, calibration)
   2.3 The labor supply elasticity as bridging parameter

3. Data and empirical setup
   - German housing transformation
   - IO table description (summary table of key ratios)
   - Calibration parameters (α, β, λ, Domar weights, elasticities)

4. Results
   4.1 Aggregate: the known result (GDP increases with η — state upfront)
   4.2 Sectoral: the novel result (non-trivial dynamics, differential sensitivity)
   4.3 Variance decomposition: closure dominates technology by factor X
   4.4 Domain of applicability: which closure for which economic situation?

5. Discussion
   5.1 Comparison with McGregor et al. (1996) long-run interpretation
   5.2 Policy implications: sensitivity of green investment estimates to closure choice
   5.3 Relation to Baqaee-Farhi framework (as computational tool, not theoretical referent)

6. Conclusion
```

**Tasks:**
1. Write new introduction with "closure dominance" framing
2. Literature review engaging with all missing references
3. Cut theoretical exposition of sections 1–2 (pages 1–11) to ~4 pages
4. Precise model description (open economy, labor mobility, calibration)
5. Results: systematic decomposition with clear tables and figures
6. Discussion: policy implications, not philosophical musings
7. Fix all minor issues (p.5 sentence, p.16 paragraph, Figure 3 axis)

**Deliverable:** Fully revised manuscript.

### Phase 5: Response & Resubmission (Week 9–10)

**Goal:** Write a compelling response to both reviewers.

**Response to Reviewer 2 (the harsh one):**
- Acknowledge the core criticism fully and graciously
- "We thank the reviewer for pointing us to the literature on alternative factor market closures, particularly Robinson (2006). We have completely revised the paper to center the labor supply elasticity as the key bridging parameter."
- Show that endogenous labor supply replaces the exogenous injection
- Show the new sectoral results and the variance decomposition
- Address every specific point (sign error, closed economy, immobile labor, CES novelty, etc.)
- Be maximally gracious — don't defend the original framing

**Response to Reviewer 1 (the constructive one):**
- Address each point systematically
- Show restructured paper with clear model description
- Include calibration summary table
- Engage with McGregor et al. (1996)
- Clarify aggregate vs sectoral distinction (now quantified)

**Deliverable:** Response letter + final submission package.

---

## 5. Timeline Summary

| Phase | Duration | Key deliverable | Status change from earlier plans |
|-------|----------|-----------------|----------------------------------|
| 1. Literature & reframing | 1–2 weeks | New intro + revised Section 2 | Same |
| 2. Code extension | 1–2 weeks | η parameter + mobile labor + variance decomp | **Shorter**: extend, don't rebuild |
| Go/no-go | End Phase 2 | Confirm sectoral results are interesting | **Earlier**: was Phase 3 |
| 3. Simulation & analysis | 1–2 weeks | Elasticity sweep + variance decomposition | **Shorter**: infrastructure exists |
| 4. Writing & restructuring | 2–3 weeks | Fully revised manuscript | Same |
| 5. Response & resubmission | 1 week | Response letter | Same |
| **Total** | **6–10 weeks** | | **Down from 12–16** |

---

## 6. Risk Assessment

### What could still go wrong

1. **Rev2 irreconcilable:** The tone ("opposite of professional scholarship," "ridicule") suggests possible hostility regardless of revision quality. **Mitigation:** Be maximally gracious, acknowledge every criticism, don't defend the original framing. The variance decomposition is a genuinely new quantitative contribution that goes beyond "we did what you told us."

2. **Editor desk-reject:** Given Rev2's recommendation, the editor may not invite a revision. **Mitigation:** Strong response letter demonstrating full engagement and clear new contribution. Consider informal contact with the editor.

3. **Sectoral results are smooth/monotonic:** If sectoral dynamics along the η continuum are uninteresting, the "sectoral dynamics" contribution is weak. **Mitigation:** The "size claim" (variance decomposition) is the insurance policy — it works regardless of whether sectoral patterns are interesting. Even if sectoral results are smooth, showing that closure choice explains 90%+ of variance while CES elasticities explain <10% is a defensible finding.

4. **Mobile labor implementation is non-trivial:** Switching from sector-specific to mobile labor requires modifying the `problem()` function's structure (the wage equation changes from `w[i]` to a single `w`). **Mitigation:** This is a well-understood CGE modification. The existing code is clean and modular.

5. **Open economy is hard:** Properly opening the CGE model is the most technically demanding task. **Mitigation:** The fallback (be transparent about the closed-economy mapping) is acceptable and Rev2 acknowledges it. The IO multiplier analysis already handles the open economy correctly.

### Strategy-specific risks

- **If the editor won't invite revision:** The revised paper is publishable elsewhere — the size claim + German housing application is a standalone empirical contribution. Metroeconomica is not the only option.
- **If η implementation proves difficult in Julia:** The model is ~50 lines of mathematics. A Python implementation using `scipy.optimize` could be done in 1 week as a fallback. But this should not be the first choice — the Julia code works and extending it is faster.

---

## 7. The Code Asset: What Already Exists

For reference, here is the complete inventory of the `(3)BeyondHulten` code base:

### Source files (`src/`)

| File | Contents | Status |
|------|----------|--------|
| `BeyondHulten.jl` | Module declaration, exports | ✅ Complete |
| `interface.jl` | Data types (CESElasticities, Shocks, Model, Data), `read_data()`, `generate_data()` | ✅ Complete |
| `solution.jl` | Solution struct, accessor functions, `multiplier()` | ✅ Complete |
| `ces.jl` | `problem()` (CES equilibrium system), `solve()` for CES, `full_labor_slack()`, `full_labor_slack_alt()`, `empirical_labor_slack()` | ✅ Working |
| `leontief.jl` | `solve()` for Leontief model, `gdp()` | ✅ Working |
| `cobbdouglas.jl` | `solve()` for Cobb-Douglas, cost functions | ✅ Working |
| `impulses.jl` | `load_impulses()`, `standard_shock()`, `impulse_shock()`, `elasticity_gradient()`, `gradient()` | ✅ Working |
| `util.jl` | `cpi()`, plotting functions | ✅ Working |
| `plots.jl` | Makie-based visualization (bar charts, error bars, panels) | ✅ Working |
| `main.jl` | Full pipeline: data → shocks → solve → plot for all configurations | ✅ Working |
| `goodwin.jl` | Dynamic model sketch (ModelingToolkit.jl) | ⚠️ Experimental |
| `ces_temporal.jl` | Temporal CES model sketch | ⚠️ Incomplete |
| `mvnrnd.jl` | MATLAB compatibility (from old B&F code) | ✅ Legacy |
| `loaders.jl` | (Empty/stub) | ⚠️ Stub |

### Notebooks

| Notebook | Purpose | Status |
|-----------|---------|--------|
| `DemandShocks.ipynb` | German data + demand shocks (main application) | ⚠️ API mismatch with current source |
| `CompareModels.ipynb` | Compare CES/Leontief/Cobb-Douglas + labor slack gradient | ⚠️ API mismatch |
| `CobbDouglas.ipynb` | Cobb-Douglas specific tests | ⚠️ API mismatch |
| `Notebooks/Analysis.ipynb` | B&F US replication (original MATLAB translation) | ✅ Reference |
| `Notebooks/Translation.ipynb` | Clean B&F US replication (Julia) | ✅ Reference |
| `Notebooks/Translation_Python.ipynb` | B&F US replication (Python) | ✅ Reference |
| `Notebooks/Summary.ipynb` | Documentation of B&F code | ✅ Reference |
| `Notebooks/Covid_Code.ipynb` | B&F 2022 COVID model analysis | ✅ Reference |

### Data

| File | Description |
|------|-------------|
| `data/I-O_DE2019_formatiert.csv` | German IO table, 71 sectors, 2019 |
| `data/impulses.csv` | Hornykewysz housing transformation shock vector |
| `data/demand_shock.csv` | Progressive sector shock results (with labor slack) |
| `data/demand_shock_no_realloc.csv` | Progressive sector shock results (no labor slack) |
| `data/demand_shock_nominal.csv` | Nominal GDP progressive shocks |
| `data/demand_shock_nominal_no_realloc.csv` | Nominal GDP, no labor slack |
| `data/simulationData.mat` | B&F US data |
| `data/us80dbasedata.mat` | B&F US data (alternative) |
| `data/stfp.mat` | B&F TFP data |

### Generated plots

| Plot | Description |
|------|-------------|
| `panel.png` | Sectoral prices & quantities, CES without labor slack |
| `panel_ls.png` | Sectoral prices & quantities, CES with full labor slack |
| `panel_ls_empirical.png` | Sectoral prices & quantities, CES with empirical labor slack |
| `impulse.png` | GDP vs elasticity gradient, no labor slack |
| `impulse_ls.png` | GDP vs elasticity gradient, full labor slack |
| `impulse_ls_empirical.png` | GDP vs elasticity gradient, empirical labor slack |
| `comparison_between_labor_slacks.png` | All scenarios compared (Fig. 7 in paper) |
| `diff_lambda_imp.png` | Domar weight changes, no labor slack |
| `diff_lambda_imp_ls.png` | Domar weight changes, full labor slack |
| `diff_lambda_imp_ls_empirical.png` | Domar weight changes, empirical labor slack |

---

## 8. Decision Point

**Recommendation:** Proceed with the salvage. The code is much further along than the earlier plans assumed, the core infrastructure exists, and the revision can be done in 6–10 weeks.

**The key go/no-go decision is at the end of Phase 2** (code extension): run a pilot with the η parameter and check whether:
1. The aggregate result confirms the expected monotonic relationship (it should)
2. The sectoral results show interesting variation along the η continuum
3. The variance decomposition shows that η dominates (ε, θ, σ)

If all three check out, proceed to full simulation and writing. If the sectoral results are boring but the size claim holds, the paper still stands on the variance decomposition contribution.

**Bottom line:** This is salvageable as a genuinely useful methodological contribution. The paper's original framing was wrong, but the underlying computational infrastructure is sound and the policy application is strong. The revision replaces the "paradigmatic cleavage" framing with a quantitative measurement of what actually matters — and the code to do this already exists.
