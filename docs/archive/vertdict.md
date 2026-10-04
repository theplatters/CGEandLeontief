• ## Verdict

Do not replace the model wholesale with Baqaee–Farhi (2022). Use a corrected 2019-style competitive CES production network as the common core, then import a simplified version of the 2022
wage-rigidity/complementarity mechanism for the unemployment closure.

The 2019 model is the cleaner base for comparing full employment, elastic labor supply, and the IO endpoint. The 2022 model is better suited specifically to involuntary unemployment and
demand-constrained labor markets, but its sector-specific sticky labor, capital, heterogeneous households, nominal rigidities, and monetary closure would substantially change the paper’s
question.

Neither model, by itself, supplies the fiscal financing closure needed for a German green-investment experiment.

## Assessment of the current implementation

The two-sector diagnostic confirms the model is not an equilibrium. At the reported prices, shock-weighted household expenditure is approximately

[
\sum_i\beta_i d_i p_i^{1-\sigma}=1.25485,
]

so the household spends about 25.5% more than its income. Consequently, Walras’ law cannot recover the omitted goods-market condition, and the missing zero-profit equation permits the
sector-1 profit residual of about 0.39.

The problem is therefore a chain:

1. src/mobile_labor.jl:148 introduces the demand shock as an unnormalized household expenditure shifter.
2. src/mobile_labor.jl:158 removes a zero-profit equation rather than a redundant goods-market equation.
3. The solver accepts convergence without checking the omitted equilibrium identities or even the solver return code at src/mobile_labor.jl:197.
4. The model uses nominal w rather than w/P in labor supply.
5. The high-(\eta) limit is misinterpreted. With (w/P<1), (L=\bar L(w/P)^\eta) converges to zero, as the diagnostic demonstrates; it does not converge smoothly to “unlimited labor.”
6. The mobile model calculates a fixed-base sum instead of using the shared Törnqvist index at src/mobile_labor.jl:218.

Thus all current mobile-labor sweep results—including the apparent price invariance and the 88.4% claim—should be treated as invalid until regenerated.

## Recommended equilibrium redesign

### 1. Separate household demand from policy investment

Do not use one demand_shock vector for two economically different experiments.

For a budget-neutral household preference reallocation, construct adjusted CES weights, for example

[
\tilde\beta_i=\frac{\beta_i d_i}{\sum_j\beta_jd_j},
\qquad
Z(p)=\sum_j\tilde\beta_jp_j^{1-\sigma},
]

and household demand

[
c_i^h
=E_h\frac{\tilde\beta_i p_i^{-\sigma}}{Z(p)}.
]

This guarantees (\sum_i p_i c_i^h=E_h).

For green investment, introduce a separate vector (g_i). Its expenditure must appear in a government/investment budget:

[
\sum_i p_i g_i=T+B+F,
]

where (T), (B), and (F) represent the chosen tax, borrowing, or foreign-financing rules. The principal and robustness financing closures should be decided in ROADMAP Phase 2 before
implementing the labor closures.

The local 2022 replication already illustrates that compositional demand shifters must be renormalized at bf_replication2/src/model.jl:442. That treatment should not, however, be mistaken
for financing an autonomous public investment increase.

### 2. Use the correct square equilibrium system

For (N) prices, (N) outputs, and one wage, retain:

- (N) zero-profit equations;
- (N-1) goods-market equations;
- one labor-market/closure equation;
- one numeraire equation.

That gives (2N+1) equations for (2N+1) unknowns.

After solving, always calculate all (N) goods-market residuals, including the omitted one. The solution is valid only if:

- all zero-profit residuals are small;
- all goods-market residuals are small;
- the labor-market residual is small;
- every institutional budget balances;
- the solver reports successful termination;
- prices, outputs, wages, and consumption are strictly positive.

Repeat the diagnostic while omitting each possible commodity market in turn. Real allocations should be invariant to that arbitrary choice.

### 3. Implement labor closures as separate types

A shared production and final-demand core should support:

- FixedLabor: (L=\bar L).
- ElasticLaborSupply: (L=\bar L[(w/P)/(w_0/P_0)]^\eta).
- FixedRealWage: (w/P=\bar\omega), with labor demand endogenous.
- UnemploymentComplementarity, for example

[
0\leq \bar L-L
;\perp;
w/P-\bar\omega\geq0.
]

The final complementarity closure is where the 2022 framework is most useful. It captures either a binding labor capacity or a binding real-wage floor without approximating the endpoint
using (\eta=10^6).

The 2022 model itself has (N) sticky labor factors rather than one mobile aggregate labor market, as shown in bf_replication2/src/network.jl:45. Its complementarity machinery can be
adapted, but its labor structure should not be copied unchanged.

## Sensitivity-analysis redesign

The existing function computes a balanced-grid main-effect share, not statistical partial (R^2). In addition, src/variance_decomposition.jl:261 silently converts an incomplete design into
a potentially nonorthogonal one.

I recommend replacing VarianceDecompositionResult with a sensitivity result containing:

- explicit factor levels and probability weights;
- first-order indices;
- total-effect indices;
- selected interaction indices;
- solver status and residuals for every evaluation;
- the assumed parameter probability measure.

For a complete factorial product distribution, compute exactly

[
S_f=\frac{\operatorname{Var}(\mathbb E[Y\mid f])}
{\operatorname{Var}(Y)}
]

using specified weights. With equal weights, describe it as a “first-order sensitivity index under a discrete uniform distribution over the selected levels.”

Do not silently discard failures:

- Retry with independent initial values and continuation.
- If any required factorial cell remains invalid, abort the decomposition.
- Report the failed region as part of the feasible parameter domain.
- Use an unbalanced regression/ANOVA only as a separately named analysis with an explicitly chosen sum-of-squares convention.

The percentage displayed at src/variance_decomposition.jl:308 is additionally renormalized over main effects. Thus 0.859 would mean 85.9% of total variance under the chosen discrete
distribution, whereas “88.4%” is merely η’s share of the sum of reported main effects. That second normalization should be removed.

Essential sensitivity tests should include:

- a purely additive analytical function;
- a pure interaction function;
- nonuniform factor weights;
- a deliberately missing factorial cell, which must raise an error;
- a solver result with a failed return code but finite .u;
- invariance to grid traversal and warm-start order.

## Recommended ROADMAP changes

ROADMAP.md:44 is fundamentally pointing in the right direction. I would update it as follows:

1. Mark recovery and audit of the two files complete, but record the audit outcome as rejection of their current results.
2. Add an explicit statement that no mobile-labor figures, GO decision, or sensitivity shares are evidentially admissible.
3. In Phase 2, distinguish:
   - preference reallocation;
   - tax-financed public investment;
   - expenditure-switching investment;
   - debt/foreign-financed net expenditure.

4. In Phase 3, add shared final_demand and institutional-budget components before implementing labor closures.
5. Replace the finite approximation to (\eta=\infty) with the explicit fixed-real-wage closure.
6. In Phase 4, make omitted-equation invariance and household expenditure exhaustion separate tests.
7. In Phase 5, require a documented product probability measure and complete designs for factorial indices.
8. Mark WORKPLAN2.md:16 as superseded: its “complete and verified,” GO decision, price-invariance finding, and 88.4% result are not supportable.

The resulting critical path should remain:

> accounting reconciliation → definition and financing of the policy shock → shared equilibrium core → labor closures → full residual validation → sensitivity analysis.

## 2019 versus 2022

Criterion 2019 model 2022 model
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ ━━━━━━━━━━━━━━━━━━━━━ ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Competitive CES/IO core Strong fit Present but much larger
─────────────────────────────────────────── ───────────────────── ──────────────────────────────────────────────
Mobile-labor closure continuum Easier to implement Not the native labor structure
─────────────────────────────────────────── ───────────────────── ──────────────────────────────────────────────
Autonomous demand and unemployment Weak; must be added Much stronger conceptual treatment
─────────────────────────────────────────── ───────────────────── ──────────────────────────────────────────────
Fiscal green-investment closure Must be added Still must be added
─────────────────────────────────────────── ───────────────────── ──────────────────────────────────────────────
Transparent IO endpoint comparison Easier Obscured by nominal and household mechanisms
─────────────────────────────────────────── ───────────────────── ──────────────────────────────────────────────
Short-run sectoral unemployment extension Limited Strong fit

The 2019 paper studies nonlinear productivity-shock propagation in efficient production-network economies, making it a suitable production-side benchmark but not a demand-multiplier
model. Econometrica article (<https://onlinelibrary.wiley.com/doi/full/10.3982/ECTA15202>)

The 2022 paper explicitly adds multiple factors, downward nominal wage rigidities, credit constraints, a zero lower bound, and Keynesian unemployment, so it is the better reference for
Closure C. American Economic Review article (<https://pubs.aeaweb.org/doi/10.1257/aer.20201229>)

My final recommendation is therefore:

> Build a corrected, budget-complete 2019-style competitive core; make labor and financing closures modular; and use a stripped-down 2022-style complementarity closure as the unemployment
> extension. Switch entirely to the 2022 base only if the paper’s main question changes from “how closure choice moves results between CGE and IO endpoints” to “how sectoral wage rigidity
> and monetary conditions govern green-stimulus effects.”

I reproduced the two-sector diagnostic, including (p_1=1.2119), (p_2=0.8220), and the collapse of employment as (\eta) grows. No repository files were changed.
