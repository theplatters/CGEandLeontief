# ADR-0022 -- Sectoral labour markets: demand-sensitive prices from the wage block

- **Status:** accepted (user ratification, 2026-09-19)
- **Date:** 2026-09-19
- **Supersedes:** —
- **Related:** ADR-0002 (naming), ADR-0004 (run records), ADR-0005 (closure
  plugins), ADR-0006 (preregistration), ADR-0010 (eta endpoints), ADR-0014
  (real-wage supply), ADR-0015 (polish), ADR-0017 (scale determinacy), ADR-0019
  (external account), ADR-0020 (the vector-wage apparatus at eta = 0), ADR-0021
  (the pinned wage vector); `docs/WORKPLAN_SENSITIVE_PRICES.md` (v4);
  `paper/equivalence.tex` (Lemma 1, the labour door, the nominal-wage no-go
  corollary); `experiments/probes/probe12_segmented_wages.jl`,
  `experiments/probes/probe13_sectoral_labour.jl`; `registry/closures.toml`
  (`labor.BETA`); `docs/ETAs.md` (candidate 1).

## Context

Under demand-only shocks and a **single** wage the price block is demand-free.
Zero profit is homogeneous of degree one in $(p, w)$ and contains no demand
term, so $p = w\,\pi(A)$, and with the CPI numeraire at $A = 1$ the calibration
identity pins $p = w = \Pi = 1$ (`paper/equivalence.tex`, Lemma 1). Every
single-wage closure — ALPHA, BETA, GAMMA, DELTA, and the retired $F = 0$ BF pin
— therefore reports the same prices in every financing cell. The only
price-moving row of the executed `matrix_5x3_v6` is the $\eta = 0$
sectoral-wage row of ADR-0020, and it moves because its per-sector wages are
tied to quantities by the sectoral first-order condition.

The labour door of `paper/equivalence.tex` (Section 6) and of `docs/ETAs.md`
(ranked candidate 1) is $N$ sectoral labour markets: a sectoral real-wage
supply curve per sector, so a wage that depends on that sector's quantity
injects demand into the price block. ADR-0020 built the apparatus (the kernel
accepts a wage vector; `problem_sectoral` is the $\eta = 0$ system). ADR-0021
proposed pinning a wage **vector** in the fixed-wage row — an *exogenous* pin,
hence price heterogeneity but not price sensitivity.

The restricted alternative — a sticky set $S$ with pinned wages plus one free
wage $w_f$ — was specified in Version 3 of the workplan and is **refuted** by
probe 12 (2026-09-19). With a single free wage, zero profit ($N$ equations)
and the CPI numeraire (one equation) already determine the $N$ prices and
$w_f$; the flexible labour-market condition never reaches the wage, it only
fixes the flexible segment's employment. Measured on the full-71 calibration:
$S = \emptyset$ nests ALPHA/BETA bit-exactly and $S$ = all nests GAMMA
bit-exactly, the external account closes to $\sim 10^{-16}$, and employment
responds to the sticky share — but $w_f \equiv 1$, `max |p - 1|` $\sim 10^{-14}$
in every financing cell, and a $+10$ % tilted pin moves prices by
$5.5\times10^{-2}$ **identically** across F1, F2 and F3.

Probe 13 (2026-09-19) prototypes the general closure on the same calibration,
with no `src/` change, and measures it (evidence below).

## Options

| Option | Mechanism | Assessment |
| --- | --- | --- |
| A. General sectoral labour markets | $N$ real-wage supply curves $L^{cm}_i = \bar L_i\,((w_i/\Pi)/(\bar w_i/\bar\Pi))^{\eta_{s,i}}$, unknowns $[p; y; w(1{:}N); F]$ | **Recommended.** Creates a demand channel to prices, nests bit-exactly onto the executed $\eta = 0$ row at $\eta_{s,i} = 0$, closes the account, and is well conditioned at moderate $\eta_s$ |
| B. Restricted sticky set (one free wage) | a pinned set $S$ plus one $w_f$ | **Rejected by measurement** (probe 12): the free wage is pinned by zero profit and the numeraire, so prices stay demand-free |
| C. Exogenous pinned wage vector | ADR-0021 option A | Kept as a design dimension (heterogeneity, allocation margin); not a sensitivity route |
| D. Supply-shock arm | ADR-0018-style $A \neq 1$ scenario | Kept as the identification complement for $\eta_s$ (`docs/ETAs.md`); a different experiment, not a repair of the demand-only matrix |

## Decision

Adopt option A. The labour axis gains a closure whose parameter is a vector of
sectoral real-wage supply elasticities, so that demand can reach prices through
the wage block. The closure is the sectoral generalisation of BETA: at
$\eta_{s,i} = 0$ for every sector it is exactly the ADR-0020 option C system
(the executed `BF` row); as $\eta_{s,i} \to \infty$ it becomes the fixed real
wage (the GAMMA corner) in prices.

Two semantics are fixed, because they are properties of the construction:

- Only the **ratios** $\bar w_i / \bar\Pi$ enter, so a uniform rescale of the
  anchors is a change of units and not an instrument. Scenarios must vary the
  elasticity **vector** (and, for the low-dimensional report, the rigid group).
- The elasticity vector is an exogenous scenario parameter. A welfare
  difference across $\eta_s$ is not a statement about wage policy, and the arm
  does not identify a labour-supply elasticity from the demand-only design
  alone: it is a sensitivity ladder over an assumed value.

## Specification

Unknowns $X = [\,p(1{:}N);\; y(1{:}N);\; w(1{:}N);\; F\,]$ — $3N + 1$ (the same
count as the $\eta = 0$ endpoint of ADR-0020).

1. $N$ zero-profit conditions, $p_i = cost_i(p, w_i)$;
2. $N$ goods-market clearings, all of them, with the cost-minimizing allocation
   $L_i = L^{cm}_i(p, y_i, w_i)$ and household wage income $\sum_i w_i L_i$;
3. $N$ sectoral labour-market conditions,
   $L^{cm}_i(p, y_i, w_i) = \bar L_i\,[\,(w_i/\Pi(p))/(\bar w_i/\bar\Pi)\,]^{\eta_{s,i}}$,
   anchors $\bar w_i/\bar\Pi = 1$ at the calibration;
4. one numeraire, $\Pi(p) = 1$.

$F$ is free, exactly as at $\eta = 0$: replacing the mobile system's single
aggregate labour equation by $N$ sectoral conditions adds $N - 1$ equations to
a block homogeneous of degree one in $(p, w, F)$, so one free scalar is needed
for consistency with the demand block (dropping a clearing equation instead is
the retired ADR-0010 shortcut).

### Endpoints

- $\eta_{s,i} = 0$ for every $i$: $L^{cm}_i = \bar L_i$, which is the
  `problem_sectoral` system (ADR-0020 option C) — the executed `matrix_5x3-v6`
  `BF` row. Probe 13 reproduces it **bit-exactly** (below).
- $\eta_{s,i} \to \infty$: the real wage tends to its anchor, so $\pi$ and $w$
  tend to the calibration and the price response vanishes — the row is
  GAMMA-like *in prices*. The limit is approached only asymptotically: the
  supply slope diverges, so the direct $3N+1$ formulation becomes singular for
  very large $\eta_s$ (measured: $\eta_s = 10^6$ fails to converge; the
  employment limit to the GAMMA closure is not attained at finite $\eta_s$).

## Probe evidence (2026-09-19, full-71 A-bill, `experiments/probes/probe13_sectoral_labour.jl`)

**Nesting at $\eta_{s,i} = 0$ against the executed `S1` (v6 `BF`) cells:**
`max |dp| = max |dy| = max |dw| = |dF| = 0.00e+00` at all three financing
columns (residuals 9.84e-12 / 4.97e-13 / 5.94e-13). The current-generation
canary holds exactly.

**Uniform-$\eta_s$ ladder** (`max |p - 1|` for F1 / F2; employment, wage range
and Jacobian condition for F2/F3):

| $\eta_s$ | `max abs(p-1)` F1 | `max abs(p-1)` F2/F3 | $\sum L$ | wage min/max | cond(J) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 0.2214 | 0.2790 | 1.00000000 | 0.970 / 1.672 | 9.35e7 |
| 0.25 | 0.1336 | 0.1529 | 1.00086363 | 0.982 / 1.357 | 6.22e7 |
| 0.5 | 0.09584 | 0.1057 | 1.00125507 | 0.987 / 1.244 | 4.65e7 |
| 1 | 0.06128 | 0.06542 | 1.00162735 | 0.992 / 1.149 | 3.10e7 |
| 2 | 0.03563 | 0.03717 | 1.00191380 | 0.995 / 1.084 | 1.85e7 |
| 5 | 0.01580 | 0.01620 | 1.00214171 | 0.998 / 1.036 | 8.41e6 |
| 10 | 0.008197 | 0.008351 | 1.00223069 | 0.999 / 1.019 | 4.40e6 |

Three readings: prices move with the financing cell at every rung (F1 distinct
from F2 ≡ F3) — the first closure in the kernel to do so; the ladder is
monotone, from the `S1` maximum toward zero; and F2/F3 financing neutrality
holds throughout (identity gap $\le 1.3\times10^{-11}$, at the solver-residual
level; $\le 4\times10^{-12}$ except the stiffest rung).

**Two-group variant** (rigid group at $\eta_{s,i} = 0$, flexible group at
$\eta_s = 0.5$):

| rigid group | share | `max abs(p-1)` F1 / F2 / F3 | $\sum L$ F1 / F2 |
| --- | ---: | --- | --- |
| programme sectors | 7/71 | 0.2151 / 0.2720 / 0.2720 | 0.996999 / 0.996621 |
| largest half by baseline employment | 36/71 | 0.1096 / 0.1234 / 0.1234 | 1.000188 / 1.000196 |

The rigidity share governs the price response and is preregisterable by rule;
identity gap $\le 1.2\times10^{-12}$.

## Consequences

- The closure is a `src/` change, so the current generation's provenance is
  superseded and the matrix is re-minted as a new generation (`matrix_5x3_v7`)
  before any sectoral cell is cited (ADR-0004). The fifteen existing cells are
  re-run in the same batch and must reproduce their manifests (the twelve
  non-BF cells to solver precision; `GAMMA`/`DELTA` to the documented
  warm-start sensitivity).
- The design schema must carry a per-**cell** elasticity vector; today the
  design pins a scalar `eta_s`. This is the same schema change that
  `docs/ETAs.md` records as the blocking item of the supply arm (ADR-0018 was
  the real-GDP decision, so the supply-shock schema ADR must take a later
  number). It should be done once, for both routes.
- Determinacy: the sectoral system is stiff at the rigid corner
  (cond $\sim 9.4\times10^7$, the documented $\eta = 0$ figure) and relaxes as
  $\eta_s$ rises ($4.4\times10^6$ at $\eta_s = 10$). The acceptance gate and the
  polish target must be set from the measured residual floor (ADR-0015),
  not assumed.
- The reported `wage` metric stays the wage-bill-weighted aggregate, with
  `wage_min` / `wage_max` carrying the dispersion (as at $\eta = 0$).
- `GAMMA`'s role: the sectoral family's $\eta_s \to \infty$ corner is the fixed
  real wage, so GAMMA can be presented as the family's elastic endpoint and the
  `BF` row as its rigid endpoint; whether that replaces the two rows or keeps
  them as matrix entries is a presentational decision for the operator.

## Enforcement

- The endpoint nesting: the general closure at $\eta_{s,i} = 0$ on the fixtures
  and on one full-71 cell must reproduce the committed `S1` (`matrix_5x3-v6`
  `BF`) manifests — the promotion's canary.
- The demand-sensitivity assertion: the price deviation must move with the
  financing cell, which no single-wage closure satisfies.
- The level invariance: a uniform rescale of the anchors leaves the real
  allocation invariant (the anchors enter only as ratios).
- The identity gate (`assert_external_account`) and the determinacy check at
  the admitted elasticity vector.
- A preregistered design carrying the `eta_s` ladder and the two-group variant
  across F1/F2/F3, plus the re-run of the fifteen existing cells.

## Open items

- The closure's registry id and its relation to `BETA` (`labor.BETA` currently
  documents the single-wage elastic-supply closure). A name must be fixed at
  promotion; the ADR-0002 taxonomy applies.
- The rigid-group rule for the low-dimensional report (a rule is preferred to a
  hand-picked list, to forestall the tuning objection).
- Whether the sectoral price response becomes a headline dimension of the
  revision or a follow-up paper.

## Implementation (promoted 2026-09-19)

The sectoral form is in the kernel (`src/`), promoted dispatch-only, so every
existing path is unchanged:

- `MobileLaborCESElasticities` gains `eta_s_vec::Union{Nothing,Vector{Float64}}`
  (default `nothing`). The 4-arg and scalar 5-arg constructors are preserved and
  set it to `nothing`, so every existing call site and test is bit-identical.
- `SectoralElasticLaborClosure(eta_s_vec)` is the new closure description
  (`src/closures/labor/types.jl`); `_closure_symbol(...) = :beta` keeps the
  closure taxonomy and the registry id unchanged (the sectoral form is BETA's
  N-market generalisation, not a new axis entry).
- `problem_sectoral` keeps its 3N+1 structure: equation 2 is now
  `log L^cm_i - log Lbar_i - eta_s,i * log(w_i / CPI) = 0`. The
  `eta_s_vec === nothing` branch is arithmetically separate, so the eta = 0
  endpoint (ADR-0020 option C) and its committed v6 cells are untouched.
- `solve`, `equilibrium_residuals` and `_equilibrium_residuals` dispatch on
  `eta == 0 || eta_s_vec !== nothing` (the same 3N+1 canonical vector, the same
  stiff tolerance ladder); a wrong-length vector throws `DimensionMismatch`
  rather than truncating. `mobile_labor_model` gains the `eta_s_vec` keyword.
- `external_balance_canary`, `market_clearing_residuals` and `gdp_components`
  already accepted the 3N+1 vector (ADR-0020) and are unchanged.

Tests (`tests/test_sectoral_labour.jl`, wired into `tests/runtests.jl`): the
scalar path still builds the 2N+2 system; the sectoral path builds the 3N+1
system and its closure type; a wrong-length vector throws; and on the full-71
calibration (guarded on the gitignored tables) the **canary** holds — the
`eta_s,i = 0` corner reproduces the eta = 0 endpoint to `< 1e-10` in
`(p, y, w, F)` — together with the demand-sensitivity assertion (a uniform
`eta_s = 0.5` moves `max|p-1|` above `1e-3`, F2/F3 remain neutral, and the
single-wage BETA closure at the same elasticity does not move prices at all)
and the external-account identity.

Registry: `labor.BETA` carries the sectoral formulation, the new symbol and
test file, the ADR link, and an open-gate entry recording that the scalar row's
demand-only degeneracy is a property of the *single-wage* form which the
sectoral form resolves.

Step 5 (done, 2026-09-19): `experiments/run.jl` carries a per-cell
`eta_s_rigid_group` (`"scalar"` / `"none"` / `"programme"` / `"largest_half"`,
with an explicit `eta_s_vec` override), the manifest records the resolved
vector (`scenario.eta_s_vec`, `scenario.eta_s_rigid_group`), and a sectoral
cell solves its own 3N+1 system. **Two promotion-plumbing defects surfaced in
the batch, each caught by an existing gate:**

- `matrix_5x3_v7`: `evaluate_gates` had no sectoral branch -- it collapsed the
  wage to a scalar and called the aggregate `labor_market_residual` -- so the
  eighteen sectoral cells failed at *gate evaluation* after a successful solve.
- `matrix_5x3_v8`: the acceptance gate `assert_external_account` collapsed the
  sectoral wage vector the same way (as did the reported `wage` metric), so the
  ADR-0019 goods-market clearing check was evaluated at a wrong common wage and
  refused the solutions (`max|market_clearing_residuals| ~ 1e-3`).

The kernel gained `sectoral_supply_gap` (with `_cpi_from_prices`); the three
scalar-wage sites in the cell path (`evaluate_gates`, `assert_external_account`,
the `wage` metric) became sectoral-aware. The **corrected generation is
`matrix_5x3_v9`** (33 cells), minted only after a full-path pre-flight that
exercises `build_cell_model -> solve_cell -> evaluate_gates` for every cell type
(twelve cells across all five labour ids and all six sectoral variants, all
`pass`). `matrix_5x3_v7` and `matrix_5x3_v8` are retained as the record of the
aborted attempts (15 executed + 18 failed each).

## Why the nesting is exact (analytical)

The promotion's canary -- the `eta_s,i = 0` corner of the sectoral family
reproduces the `η = 0` endpoint (ADR-0020 option C) -- is not a numerical
coincidence. The two closures are the same equations.

Write the sectoral system in the unknowns `X = [p; y; w; F]` as four blocks:

  P (N):  p_i - c_i(p, w_i) = 0                            [zero profit, unit cost]
  Y (N):  y_i - int_i(p, y) - fin_i(p, y, L; F) - g_i = 0  [goods-market clearing]
  W (N):  log L^cm_i(p, y, w_i) - log Lbar_i - eta_s,i * log(w_i / P) = 0
  E (1):  Pi(p) - 1 = 0                                    [CPI numeraire]

The `η = 0` endpoint is the same system with two changes: its W block drops the
supply term, and its Y block evaluates labour income at the *frozen*
allocation, `L_i := Lbar_i`, instead of at the cost-minimizing demand `L^cm_i`.

**Proposition.** With `eta_s,i = 0` for every sector the two systems have the
same solution set.

**Proof.** (i) Let `X` solve the sectoral system at `eta_s ≡ 0`. Its W block
gives `L^cm_i(X) = Lbar_i` for every i, the supply term vanishing identically.
Substituting `L^cm_i = Lbar_i` into its Y block makes that block equal the
endpoint's Y block (P and E were already identical), so `X` solves the endpoint
system. (ii) Let `X` solve the endpoint system. Its W block gives
`L^cm_i(X) = Lbar_i`, which satisfies the sectoral W block at `eta_s ≡ 0`, and
substituting into the Y block yields the sectoral Y block; so `X` solves the
sectoral system. QED

Both systems are square (`3N+1` unknowns, `3N+1` equations) with the same
residual map on the common solution set, so the solutions coincide -- the
canary measured `0.00e+00` in all four blocks when the sectoral solve is
warm-started from the endpoint (and `< 1e-10` cold).

**The intuition, stated once for the record.** In both closures labour is
immobile, but for different reasons. The `η = 0` endpoint freezes the
*allocation*: `L_i` is a datum, and the sectoral wage is whatever clears
cost-minimizing demand against it. The sectoral closure at `eta_s,i = 0` has a
*vertical supply curve* in every sector (`L^s_i = Lbar_i` for every `w_i`), so
the market clears at the same allocation with a demand-determined wage. A
constraint imposed on a quantity and an infinitely inelastic price response are
the same row of algebra. The common content is that *quantities are rigid and
prices absorb*: the N wages -- and hence `p = pi(A, wage structure)` -- are
free, which is exactly why this corner is demand-sensitive while the
single-wage BETA form is not.

**The family as an interpolation.** Sweeping `eta_s,i` moves the adjustment
from quantities to prices: at `eta_s,i = 0` the allocation is pinned and the N
wages absorb; as `eta_s,i -> inf` the supply curves become infinitely elastic,
the real wage pins to its anchor and employment absorbs. The two corners are
therefore the matrix's own labour-rigidity rows -- `BF` (allocation frozen,
wages free) and `GAMMA` (real wage pinned, allocation free) -- so the sectoral
family is a one-parameter path *between* them, with the uniform `eta_s` as its
canonical representative. The `GAMMA` end is a limit, not a cell: the
formulation stiffens and the direct solve eventually fails (`eta_s = 1e6`
singular, probe13), so it is approached asymptotically, never attained.

**Parallel with `BETA ≡ ALPHA`.** That equivalence collapses BETA onto ALPHA
because the *real-wage anchor binds* (`w/P = 1`), making the supply curve return
`Lbar` whatever its elasticity. Here the same curve returns `Lbar` because the
*elasticity vanishes* rather than the price being pinned. Two degenerations of
one supply family onto one fixed quantity -- both exact, both visible in the
kernel as a branch that never executes.

## Amendment 2026-09-20 — the C1 measurement defect and the matrix_5x3_v10 re-mint

**Root cause.** `gdp_components` (`src/core/diagnostics.jl`) selected the
sectoral wage vector only at `η == 0`
(`w = elasticities.η == 0.0 ? sol.wages_raw : sol.wages_raw[1]`), but every
sectoral cell of this ADR runs at `η = 1` with `eta_s_vec !== nothing` — so
all eighteen `matrix_5x3-v9-BETA-*` sectoral manifests measured GDP at the
scalar wage `w_1` with `L_i` recomputed at that wage. Wrong in those
manifests: `gdp`, `gdp_rel`, `gdp_expenditure(_rel)`, `gdp_deflator`,
`gdp_wedge` and the `gdp_c` / `gdp_m_final` diagnostics. The solve was
unaffected (prices, quantities, employment, consumption identical).

**Disproved claim.** The Implementation section above states that
"`external_balance_canary`, `market_clearing_residuals` and `gdp_components`
already accepted the 3N+1 vector (ADR-0020) and are unchanged" — true of the
first two, false of the third: accepting the vector is not using it.

**Fix (commit `8774595`).** The wage selection now mirrors `sect`
(`η == 0 || eta_s_vec !== nothing`); `experiments/run.jl` `evaluate_gates`
hard-asserts the companion identity `gdp_wedge = -canary_diff` at 1e-9
(measured max 4.4e-16 over the new generation); tests gain a data-free
open-fixture regression, the real-table identity on the solved sectoral cell,
and a sectoral smoke cell through the harness.

**Re-mint (`matrix_5x3_v10`).** Design committed at `979a58e` (preregistered,
sha256 `ab5974a3…`); all 33 cells executed at `c2c5b7d` (actor `calculato`),
every gate passes. Measured corrections (v10 vs v9 manifests): `gdp_rel`
flips sign in 15 of the 18 sectoral cells — `BETA-F1-etas05` −0.4414 % →
+0.1217 %, `BETA-F2-rigidprog` −1.4525 % → −0.3414 %,
`BETA-F1-rigidhalf` −0.7748 % → +0.0144 %; `gdp_wedge` drops from up to
3.8e-3 to ≤ 1.5e-12, i.e. the ADR-0019 companion identity now closes in the
sectoral cells; deflator corrections are small but real (F2-etas025 1.004096
→ 1.004088; F2-rigidhalf 1.006294 → 1.006281). The fifteen matrix cells
reproduce v9 to solver precision, not bit-for-bit: v6 and v9 are bit-identical
to each other (both executed in another working copy), while v10 ran in this
working copy, so cross-copy floating-point noise appears at ≤ 9.4e-9 (the
stiff BF-F1 `wage_max` diagnostic), ≤ 1.7e-11 in the ALPHA/BETA/GAMMA/DELTA
rows, and ≤ 4.6e-15 where the equilibrium is the exact baseline — against the
v10 design file's "expected bit-for-bit" comment, which was written from the
same-environment v6→v9 identity and is not met across environments.
Solutions are unchanged: v10 reproduces this working copy's own pre-batch
solves bit-for-bit, and all non-gdp metrics differ only at the noise level.

**Citation rule.** `matrix_5x3_v9` stays citable for everything except the
eighteen sectoral cells' gdp-family metrics; for those cite `matrix_5x3_v10`
(flow table `paper/tables/matrix_5x3_v10_flows.md`, created in parallel; the
v9 flow table carries a supersession banner).