# Labor closures

The package exposes three dimensions that should not be conflated:

* **Legacy `CES` and `CobbDouglas`** store an exogenous `labor_slack` callback
  (or their compatibility `Symbol`). `labor_closure` reports this as an
  `ExogenousLaborClosure`; it is not a wage-regime selector.
* **Mobile η** controls geometric intersectoral allocation:
  `L_fixed^(1-η) * L_costmin^η`.
* **The wage regime** is flexible or fixed.  `:mobile` (or
  `FlexibleWageClosure()`) clears the labor market with a common wage; `:fixed`
  (or `FixedWageClosure()`) fixes wages at one.  In the fixed regime the
  employment gap (`labor_bar - sum(L_i)`) can be computed after solving but is
  not stored in `Solution` and is not an equilibrium constraint. Employment
  can therefore exceed `labor_bar`; `:fixed` is not a capped unemployment model.

Use `labor_closure(model)` (or `labor_closure(model.options)`) to inspect the
canonical taxonomy and
`equilibrium_residuals(solution)` for closure-appropriate diagnostics.

The legacy callback, η, and the wage regime are independent concepts: η=0 is
immobility and η=1 is full cost-minimizing allocation. Values outside `[0,1]`
are extrapolations; η is not a labor-supply elasticity or `L̄*w^η`.

## Standard partial-mobility formulation

The conventional static treatment of limited intersectoral labor mobility is
to give each sector its own wage and make sectoral labor supply respond to
relative wages. For example,

```math
L_i^s = \bar L
\frac{s_i^0 (w_i/\bar w_i)^\nu}
     {\sum_j s_j^0 (w_j/\bar w_j)^\nu},
```

where `s_i^0` is the baseline employment share, `\bar w_i` is the baseline
sectoral wage, and `\nu >= 0` is an intersectoral mobility elasticity. This
allocation always satisfies `\sum_i L_i^s = \bar L`. At `\nu = 0`, employment
shares are fixed. As `\nu` becomes large, employment becomes highly responsive
to wage differences and adjusted sectoral wages approach equality in an
interior equilibrium.

Each sector then minimizes production cost at its own wage. Labor demand follows
from Shephard's lemma,

```math
L_i^d = y_i \frac{\partial c_i(\mathbf p,w_i,A_i)}{\partial w_i},
```

and the model imposes `L_i^d = L_i^s` for every sector. With prices, quantities,
and sectoral wages as unknowns, the equilibrium has `N` zero-profit equations,
`N-1` goods-market equations, `N` sectoral labor-market equations, and one
numeraire equation: `3N` equations for `3N` unknowns.

The standard endpoints are:

| Closure | Labor allocation | Wages |
| --- | --- | --- |
| Immobile | `L_i = \bar L_i` | Sector-specific shadow wages |
| Partially mobile | Wage-responsive sectoral supply | Sector-specific wages |
| Fully mobile | `\sum_i L_i = \bar L` | Common wage |

This differs from the current `MobileLaborCES` interpolation. The current model
uses a common wage while geometrically interpolating between fixed and
cost-minimizing labor, then applies an exponential-quadratic efficiency factor
to the flexible-wage cost equation. That factor is a project-specific
reduced-form assumption, not a labor wedge derived from the CES first-order
conditions. In particular, its curvature
`factor_share_i * (1-factor_share_i) * abs(1-ϵ)/ϵ` guarantees a nonnegative
penalty but is not established as the exact CES allocative-loss coefficient.

For a structural partial-mobility interpretation, prefer sectoral wages and the
labor-supply system above; the output loss then follows from constrained
reallocation without an additional productivity factor. If the mechanism is a
tax or another explicit distortion instead, put a wedge in the labor condition,
for example `p_i MPL_i = (1+τ_i)w`, and account for the resulting revenue or
rents. For adjustment over time, use dynamic worker mobility or convex
employment-adjustment costs rather than a static productivity penalty.

Related references include Hsieh and Klenow (2009),
[doi:10.1162/qjec.2009.124.4.1403](https://doi.org/10.1162/qjec.2009.124.4.1403),
on factor wedges and misallocation; Baqaee and Farhi (2020),
[doi:10.1093/qje/qjz030](https://doi.org/10.1093/qje/qjz030), on productivity
and misallocation in general equilibrium; and Artuç, Chaudhuri, and McLaren
(2010), [doi:10.1257/aer.100.3.1008](https://doi.org/10.1257/aer.100.3.1008),
on dynamic labor mobility costs.

```julia
legacy = CES(elasticities; labor_slack=full_labor_slack)
mobile = mobile_labor_model(data, shocks, .5, .5, .9, 1.0;
                            closure=:mobile)
fixed = MobileLaborCES(mobile.options.elasticities,
                       sum(data.labor_share), FixedWageClosure())
```

`VarianceDecompositionResult` remains as a deprecated alias of
`SobolResult`; new code should use `SobolResult` and
`eta_sweep_diagnostics`.
