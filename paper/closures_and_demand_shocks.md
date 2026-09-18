# Macroeconomic Closures, Neo-Keynesian Rigidities, and Demand Shocks

**Context:** Research note on macroeconomic closures, the treatment of demand shocks in Computable General Equilibrium (CGE) versus macroeconometric input-output models (specifically the EU JRC FIDELIO model), and their precise mapping into the `BeyondHulten` 5×3 closure matrix.

---

## 1. The Core Theoretical Cleavage: Demand Shocks in CGE vs. IO Models

A central question in applied input-output modeling—especially for evaluating green public investment programs such as the EU Recovery and Resilience Facility (RRF) or large-scale building renovation—is whether an exogenous demand stimulus:
1. **Expands real aggregate output** by mobilizing unutilized resources (the Keynesian / Leontief perspective), or
2. **Merely reallocates scarce resources**, causing factor prices to rise and crowding out private activity with virtually zero net change in aggregate real GDP (the canonical neoclassical CGE perspective).

### Why Standard CGE Models Produce Minimal Aggregate Effects from Demand Shocks
In textbook Walrasian CGE models:
* **Full factor utilization:** Aggregate primary factor endowments (labor $\bar{L}$, capital $\bar{K}$) are exogenously given and fully employed.
* **Instantaneous market clearing:** Prices and wages are perfectly flexible and clear all commodity and factor markets instantaneously.
* **Say's Law operates:** Output is supply-determined by the production possibility frontier $Y = F(K, L, \mathbf{X})$.
* **Crowding out:** An exogenous increase in public investment or government spending increases demand in targeted sectors. However, because total labor and capital cannot expand, the extra factor demand bids up wages and capital rentals. Relative price increases squeeze the margins of non-shocked sectors and crowd out private consumption, private investment, or net exports. Aggregate real GDP response is negligible ($\Delta Y \approx 0$).

---

## 2. The Macroeconometric Alternative: The EU JRC FIDELIO Model

The **FIDELIO** model (*Fully Interregional Dynamic Econometric Long-term Input-Output* model, JRC / DG GROW; Rocchi et al., 2019, 2025) was built specifically to overcome this limitation. Rather than adopting a Walrasian general equilibrium core, FIDELIO is a **dynamic macroeconometric, multi-sectoral Neo-Keynesian model** grounded in the official Eurostat/JRC FIGARO inter-country input-output tables (64 NACE industries, 46 countries/regions).

FIDELIO generates positive, non-trivial real output and employment multipliers from demand shocks through five structural features:

1. **Involuntary Unemployment & the Wage-Setting Curve:**
   Wages do *not* clear the labor market. The baseline economy features involuntary unemployment ($u > 0$). Nominal wages adjust sluggishly along an econometrically estimated **wage curve** (Blanchflower & Oswald tradition):
   $$\Delta \ln(w_{i,r,t}) = f(\text{inflation}_{t-1}, \text{productivity}_{t}, u_{r,t-1})$$
   Because labor supply is not vertical at baseline, firms can hire idle workers to expand production without hitting an immediate supply wall.

2. **Mark-up Pricing and Sluggish Cost Pass-Through:**
   Product prices are set via markups over unit costs (derived from flexible Translog cost functions over capital, labor, energy, and intermediate inputs) rather than jumping to Walrasian clearing levels. Output in the short-to-medium run is **demand-determined**.

3. **Keynesian Multiplier Loops (The ECM Consumption Block):**
   Household consumption is not governed by Ricardian equivalence or an intertemporal Euler equation. Instead, consumption follows a dynamic **Error Correction Model (ECM)** tied to current real disposable income:
   $$\text{Demand Shock } \Delta G \implies \uparrow Y \implies \uparrow L \implies \uparrow (w \cdot L) \implies \uparrow C^{\text{household}}$$
   This induces a secondary Keynesian expenditure loop (*Type II multiplier*), yielding aggregate output multipliers between **1.2 and 1.8** in JRC policy simulations.

4. **Dynamic Capital Accumulation (Endogenous Capacity):**
   Capital is not a fixed static endowment. Investment responds dynamically through an accelerator / user-cost formulation: sustained demand increases capacity utilization and expected returns, prompting firms to invest and expand capital stocks over time.

5. **Medium-to-Long-Run Re-equilibration:**
   In the longer run, as unemployment drops and capacity tightens, the wage curve drives up wages and unit costs. The resulting rise in output prices erodes export competitiveness, generating gradual crowding out and steering the economy toward a medium-run NAIRU—unless the stimulus was in public capital or R&D, in which case it permanently shifts potential output outward.

---

## 3. Mapping into the `BeyondHulten` 5×3 Closure Matrix

The `BeyondHulten` platform organizes macroeconomic closures across two dimensions (5 Labour Closures $\times$ 3 Financing Closures):

### The 5 Labour Closures
1. **`ALPHA` (Full-employment mobile CGE benchmark):** $\sum_i L_i = \bar{L}$, single flexible wage $w$, cost-minimizing intersectoral allocation ($\eta = 1$).
2. **`BF` (Immobile labour endpoint selector):** $\eta = 0$ freezes sectoral labor quantities to baseline shares $L_i = \bar{L}_i$, single flexible wage $w$.
3. **`BETA` (Elastic labour supply on the real wage):** $\sum_i L_i = \bar{L} \left(\frac{w/P}{w_0/P_0}\right)^{\eta_s}$ deflated by the CPI numeraire.
4. **`GAMMA` (Fixed real wage, uncapped extensive margin):** $w/P = \bar{w} = 1$, employment $L = \sum_i L_i$ is completely demand-determined.
5. **`DELTA` (Leontief IO corner):** `GAMMA` ($w/P = 1$) combined with the zero-substitution limit of the CES production core ($\theta, \epsilon, \sigma \to 0^+$), reproducing the fixed-price Leontief multiplier.
*(Complementarity extension `ZETA`: $0 \le \bar{L} - L \perp w/P - \bar{\omega} \ge 0$, regime-switching between `GAMMA` slack and `ALPHA` full employment).*

### The 3 Financing Closures
1. **`F1` (Preference reallocation / demand tilt):** Budget-neutral composition shift within fixed household expenditure ($\sum_i p_i c_i = E_h$). No net fiscal injection.
2. **`F2` (Tax-financed public investment):** Balanced government budget ($\sum_i p_i g_i = T(p)$). Autonomous public investment financed by an explicit domestic tax counterparty.
3. **`F3` (External / debt-financed programme):** Autonomous public spending financed externally or via debt ($\sum_i p_i g_i = F$), with imports and external balances absorbing the deficit.

### Matrix Location of Models

```
                      F1 (Demand Tilt)         F2 (Tax-Financed)       F3 (Debt / External)
  ┌─────────────────┬────────────────────────┬───────────────────────┬────────────────────────┐
  │ ALPHA (Full L)  │ Sectoral reallocation  │ Full Crowding Out     │ Crowding Out (w bids ↑)│  <-- Standard CGE
  │ BF    (Immobile)│ Bottlenecks / wedges   │ Contractionary wedge  │ Bottleneck limits      │
  │ BETA  (Elastic) │ Moderate shift         │ Mild expansion        │ Expansion via w/P      │
  │ GAMMA (Fixed w) │ Demand-driven shift    │ Balanced-budget mult. │ ★ KEYNESIAN MULTIPLIER │  <-- FIDELIO (Short run)
  │ DELTA (IO/Leont)│ IO composition shift   │ Leontief balanced mult│ ★ PURE IO MULTIPLIER   │  <-- FIDELIO (Modules 1-2)
  └─────────────────┴────────────────────────┴───────────────────────┴────────────────────────┘
```

* **Where Standard CGE Sits:** At **`ALPHA-F1`** or **`ALPHA-F2`**. Fixed aggregate labor $\bar{L}$ and flexible wages force crowding out; real output effects are close to zero.
* **Where FIDELIO Sits:**
  * **Static / Basic IO Mode (Modules 1–2):** Corresponds to **`DELTA-F3`** (a pure Leontief input-output multiplier driven by external funding).
  * **Short-Run Dynamic Mode:** Corresponds to **`GAMMA-F3`** (fixed real wage, unconstrained extensive margin, autonomous injection), activating the full Keynesian respending multiplier.
  * **Medium-to-Long Run:** As wages respond to falling unemployment via the wage curve, FIDELIO moves dynamically from **`GAMMA-F3` toward `BETA-F3`**, where real wage growth dampens further expansion.

---

## 4. How Neo-Keynesian Assumptions Map into Our Closures

In `BeyondHulten`, Neo-Keynesian assumptions are not auxiliary parameters—they define the mathematical structure of the closure:

### A. Wage Rigidity (Extensive Margin vs. Wage Inflation)
* **Real Wage Anchor:** In accordance with ADR-0014 and DE-0004, wage rigidity in our model is defined on the **real wage** ($w/P = \bar{w}$), deflated by the CPI numeraire ($P = 1$).
* **Mechanism:** In `GAMMA`, fixing $w/P = 1$ removes the wage as an equilibrating price. The labor market equation changes from clearing $L = \bar{L}$ to an unconstrained demand equation $L = \sum_i L_i(y_i, \mathbf{p})$. The extensive margin of employment absorbs the demand shock.
* **Asymmetric Rigidity (`ZETA`):** In the spirit of Baqaee & Farhi (2022), downward wage rigidity creates a regime switch: demand deficits cause unemployment (`GAMMA` regime), whereas demand surges that exhaust slack trigger wage inflation (`ALPHA` regime).

### B. Price Rigidity and Fixed Markups
* In standard New-Keynesian models, Calvo or Rotemberg pricing prevents prices from adjusting instantaneously to demand shocks.
* In `BeyondHulten`:
  * Under `GAMMA`, because wages are pegged ($w = 1$) and the shock is demand-driven, marginal costs remain roughly constant. Product prices naturally exhibit high inertia ($\mathbf{p} \approx 1$).
  * In **`DELTA`**, taking the zero-substitution limit ($\theta, \epsilon, \sigma \to 0^+$) prevents firms from adjusting technical factor proportions. Cost pass-through becomes purely linear and additive, reproducing the classic fixed-price Leontief/Robinson (2006) multiplier.

### C. Non-Ricardian Demand & Disposable Income Loops
* Standard CGE models enforce intertemporal optimization where households anticipate future tax liabilities, neutralizing debt-financed spending.
* In `BeyondHulten`, the consumption block routes household income $E_h = (1 - \tau) w \sum_i L_i + \text{capital income}$. Under `F3` (external/debt financing), the absence of an immediate tax hike $\tau$ allows the newly generated wage income $w \cdot \Delta L$ to directly fuel private consumption ($(1-s)E_h$), generating an endogenous Keynesian multiplier.

---

## 5. Summary and Methodological Takeaway for the Paper

1. **The Binary Wage Regime Dominates:** The literature often debates the elasticity of factor substitution or intersectoral mobility ($\eta$). As shown in `docs/definitive_guide.md` and our model runs, intersectoral mobility $\eta$ produces second-order reallocation effects, whereas **the wage regime (flexible `ALPHA` vs. sticky `GAMMA`) is the first-order determinant of the aggregate multiplier**.
2. **Reconciling Applied Discrepancies:** When policy studies using macroeconometric models (like FIDELIO or E3ME) report output multipliers of $1.5$ for green investments while standard CGE models report zero or negative multipliers, the divergence is not an empirical dispute over green technology. It is a direct mathematical consequence of evaluating the shock at **`GAMMA/DELTA-F3`** (slack resources, wage stickiness, debt financing) versus **`ALPHA-F2`** (full employment, market clearing, domestic tax financing).
