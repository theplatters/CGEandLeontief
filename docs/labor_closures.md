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
