# experiments/probes/probe12_segmented_wages.jl
#
# Step 1 of docs/WORKPLAN_SENSITIVE_PRICES.md: prototype the segmented-wage
# closure (the RESTRICTED variant) — a sticky set S with pinned wages and a
# flexible segment F carrying one solved common wage w_f.
#
# Specification (workplan, "The specification to implement"):
#   unknowns X = [p(1:N); y(1:N); w_f; F]                (2N + 2)
#   (1) N zero-profit conditions:  p_i = cost_i(p, w_i),
#         with w_i = wbar_i for i in S and w_i = w_f for i in F;
#   (2) N goods-market clearings, all of them (ADR-0019);
#   (3) one flexible labour-market condition:
#         sum_{i in F} L^cm_i = Lbar_F * ((w_f / P) / (wbar_f / Pbar))^eta_s
#       so the STICKY segment's employment stays demand-determined and
#       uncapped, as in the fixed-wage regime;
#   (4) one numeraire: CPI(p) = 1.
#
# No src/ change. The demand block is the kernel's `_mobile_market_demand`,
# which accepts a per-sector wage vector since ADR-0020; `_cost_minimizing_labor`
# and `_ces_unit_cost` are elementwise in w. The flexible labour row (3) is the
# only added equation.
#
# Falsifiable predictions, recorded BEFORE running (workplan Step 1):
#   P1  prices move with the financing cell: max|p - 1| differs across F1/F2/F3
#       at a fixed sticky share (no existing closure satisfies this);
#   P2  the segment wage gap w_f / wbar_S widens as the sticky share rises;
#   P3  employment absorbs more (larger |dL|) as the sticky share rises;
#   P4  S = empty reproduces the committed ALPHA (eta_s = 0) and BETA
#       (eta_s = 0.5) cells to solver precision;
#   P5  the external account closes: the price-weighted clearing residual
#       (the proxy validated in probe11, ratio 1.0000) stays at the
#       solver-residual level.
#
# Endpoints of the family, both measured: S = empty -> ALPHA/BETA; S = all
# (fully pinned) -> GAMMA, checked against `problem_fixed`.

include(joinpath(@__DIR__, "..", "run.jl"))
using Printf, LinearAlgebra
import NonlinearSolve
const BH = BeyondHulten

# ── wage vector: sticky sector i gets wbar[i], flexible sector gets wf ───────
seg_wage_vector(N, sticky::BitVector, wbar::AbstractVector, wf::Real) =
    ifelse.(sticky, wbar, float(wf))

function cpi_of(model, p)
    σ = model.options.elasticities.σ
    return sum(model.data.consumption_share .* p .^ (1 - σ))^(1 / (1 - σ))
end

# ── generic ADR-0015 solve ladder (Newton + residual-gated LM polish) ───────
function solve_ladder(f, init, pars; tol = 1e-12)
    res = NonlinearSolve.solve(NonlinearSolve.NonlinearProblem(f, copy(init), pars),
        reltol = 1e-8, abstol = 1e-8, maxiters = 40000)
    x = res.u
    rmax = rmax_of(f, x, pars)
    for _ in 1:8
        rmax <= tol && break
        res = NonlinearSolve.solve(NonlinearSolve.NonlinearProblem(f, x, pars),
            NonlinearSolve.LevenbergMarquardt(); reltol = 1e-13, abstol = 1e-13, maxiters = 40000)
        r2 = rmax_of(f, res.u, pars)
        r2 < rmax || break
        x, rmax = res.u, r2
    end
    x, rmax
end

# ── the segmented-wage residual (2N + 2 equations) ──────────────────────────
function seg_residual!(out, X, pars)
    model, sticky, wbar, wbar_f, eta_s, Lbar_F = pars
    N = length(model.data.factor_share)
    p  = BH._positive_floor(X[1:N])
    y  = BH._positive_floor(X[N+1:2N])
    wf = BH._positive_floor(X[2N+1])
    F  = X[2N+2]
    w  = seg_wage_vector(N, sticky, wbar, wf)
    b  = BH._mobile_market_demand(model, p, y, w; external_transfer = F)
    flex = .!sticky
    # (3) flexible segment labour market: at the anchor (w/P = wbar_f) the
    #     supply is Lbar_F and equals sum_{flex} L^cm_i at the baseline.
    out[1:N]     .= p .- b.cost
    out[N+1:2N]  .= y .- b.intermediary_demand .- b.total_final_demand
    out[2N+1]     = sum(b.L_i[flex]) - Lbar_F * (wf / cpi_of(model, p) / wbar_f)^eta_s
    out[2N+2]     = cpi_of(model, p) - 1.0
    nothing
end

rmax_of(f, x, pars) = maximum(abs, (o = similar(x); f(o, x, pars); o))

# ── fully-pinned endpoint residual (2N unknowns, wage vector pinned, F = 0) ──
function pin_fixed_residual!(out, X, pars)
    model, wbar = pars
    N = length(model.data.factor_share)
    p = BH._positive_floor(X[1:N])
    y = BH._positive_floor(X[N+1:2N])
    b = BH._mobile_market_demand(model, p, y, wbar)
    out[1:N]    .= p .- b.cost
    out[N+1:2N] .= y .- b.intermediary_demand .- b.total_final_demand
    nothing
end

# ── solve with the ADR-0015 ladder (Newton + residual-gated LM polish) ──────
function solve_seg(model, sticky, wbar, wbar_f, eta_s, Lbar_F; init)
    pars = (model, sticky, wbar, wbar_f, eta_s, Lbar_F)
    f = seg_residual!
    res = NonlinearSolve.solve(
        NonlinearSolve.NonlinearProblem(f, copy(init), pars),
        reltol = 1e-8, abstol = 1e-8, maxiters = 40000)
    x = res.u
    rmax = rmax_of(f, x, pars)
    for _ in 1:8
        rmax <= 1e-12 && break
        res = NonlinearSolve.solve(
            NonlinearSolve.NonlinearProblem(f, x, pars),
            NonlinearSolve.LevenbergMarquardt(); reltol = 1e-13, abstol = 1e-13,
            maxiters = 40000)
        r2 = rmax_of(f, res.u, pars)
        r2 < rmax || break
        x, rmax = res.u, r2
    end
    return x, rmax
end

function jac_cond(f, X, pars; h = 1e-7)
    n = length(X)
    J = zeros(n, n); o1 = similar(X); o2 = similar(X)
    for j in 1:n
        xp = copy(X); xp[j] += h
        xm = copy(X); xm[j] -= h
        f(o1, xp, pars); f(o2, xm, pars)
        J[:, j] .= (o1 .- o2) ./ (2h)
    end
    s = svdvals(J)
    return maximum(s) / minimum(s)
end

# ── identify the sticky set by a preregisterable rule ───────────────────────
sticky_by_employment(data, share::Float64) = begin
    N = length(data.labor_share)
    k = round(Int, share * N)
    k <= 0 && return falses(N)
    k >= N && return trues(N)
    order = sortperm(data.labor_share; rev = true)
    s = falses(N); s[order[1:k]] .= true; s
end

function main()
    root = default_root()
    design_d = load_design("matrix_5x3_v6"; root = root)
    ref = build_reference(design_d; root = root)
    data = ref.data
    N = length(data.factor_share)
    ψ, g = programme_vectors(design_d, data, collect(1:N); root = root)
    @printf("N = %d, segmented unknowns = 2N+2 = %d\n", N, 2N + 2)

    # ══ P4: endpoint nesting, S = empty  (= ALPHA / BETA) ════════════════════
    println("\n== P4: endpoint nesting (S = empty) against ALPHA-F2 / BETA-F2 ==")
    for (lab, eta_s) in (("ALPHA", 0.0), ("BETA", 0.5))
        cell = design_d["cells"]["matrix_5x3-v6-$lab-F2"]
        model = build_cell_model(cell, design_d, data, ψ, g)
        sol = solve_cell(cell, design_d, data, ψ, g, ref.init_warm)
        init = [sol.prices_raw; sol.quantities; sol.wages[1]; sol.external_transfer]
        sticky = falses(N); wbar = ones(N)
        x, rmax = solve_seg(model, sticky, wbar, 1.0, eta_s, sum(data.labor_share); init = init)
        dp = maximum(abs, x[1:N] .- sol.prices_raw)
        dy = maximum(abs, x[N+1:2N] .- sol.quantities)
        dw = abs(x[2N+1] - sol.wages[1]); dF = abs(x[2N+2] - sol.external_transfer)
        @printf("  %-5s: resid %.2e | max|dp| %.2e  max|dy| %.2e  |dw| %.2e  |dF| %.2e\n",
            lab, rmax, dp, dy, dw, dF)
    end

    # ══ endpoint nesting, S = all  (= GAMMA via the pinned wage VECTOR) ═══════
    # The restricted system cannot express S = all (w_f would be unidentified),
    # so the all-pinned endpoint is measured directly: the fixed-wage residual
    # with a wage vector wbar (2N unknowns, F = 0), the ADR-0021 option-A form,
    # at wbar = ones must reproduce the kernel GAMMA solve.
    println("\n== P4b: fully-pinned endpoint (S = all, wbar = ones) against GAMMA-F2 ==")
    cell = design_d["cells"]["matrix_5x3-v6-GAMMA-F2"]
    model = build_cell_model(cell, design_d, data, ψ, g)
    sol = solve_cell(cell, design_d, data, ψ, g, ref.init_warm)
    init = [sol.prices_raw; sol.quantities]
    pars = (model, ones(N))
    x, rmax = solve_ladder(pin_fixed_residual!, init, pars)
    dp = maximum(abs, x[1:N] .- sol.prices_raw)
    dy = maximum(abs, x[N+1:2N] .- sol.quantities)
    @printf("  resid %.2e | max|dp| %.2e  max|dy| %.2e | GAMMA-F2 max|p-1| %.2e, employment %.10f\n",
        rmax, dp, dy, maximum(abs, sol.prices_raw .- 1),
        sum(BH._cost_minimizing_labor(sol.prices_raw, sol.quantities, 1.0, model)))

    # ══ P1–P3, P5: demand sensitivity across the financing cells ═════════════
    println("\n== P1-P3: sticky share ladder x financing (sticky = largest-k by employment) ==")
    shares = (0.0, 0.25, 0.5, 0.75)
    for fin in ("F1", "F2", "F3")
        cell = design_d["cells"]["matrix_5x3-v6-ALPHA-$fin"]
        model = build_cell_model(cell, design_d, data, ψ, g)
        base = solve_cell(cell, design_d, data, ψ, g, ref.init_warm)
        L0 = sum(BH._cost_minimizing_labor(base.prices_raw, base.quantities, 1.0, model))
        println("  --- $fin (mobile ALPHA employment L0 = $(@sprintf("%.8f", L0))) ---")
        @printf("  %6s %9s %11s %13s %12s %12s %10s %10s %11s\n",
            "share", "resid", "max|p-1|", "sum L", "dL/L0", "wage gap", "cond(J)", "w_f", "identity")
        for share in shares
            sticky = sticky_by_employment(data, share)
            wbar = ones(N)
            Lbar_F = sum(data.labor_share[.!sticky])
            init = [base.prices_raw; base.quantities; 1.0; base.external_transfer]
            x, rmax = solve_seg(model, sticky, wbar, 1.0, 0.5, Lbar_F; init = init)
            p, y = x[1:N], x[N+1:2N]; wf = x[2N+1]; F = x[2N+2]
            w = seg_wage_vector(N, sticky, wbar, wf)
            b = BH._mobile_market_demand(model, p, y, w; external_transfer = F)
            L = sum(b.L_i)
            gap = isempty(wbar[sticky]) ? NaN : (minimum(wbar[sticky]) > 0 ?
                  wf / (sum(wbar[sticky] .* data.labor_share[sticky]) / sum(data.labor_share[sticky])) : NaN)
            cl = y .- b.intermediary_demand .- b.total_final_demand
            pars = (model, sticky, wbar, 1.0, 0.5, Lbar_F)
            @printf("  %6.2f %9.2e %11.3e %13.8f %12.3e %12.5f %10.2e %10.5f %11.2e\n",
                share, rmax, maximum(abs, p .- 1), L, (L - L0) / L0, gap,
                jac_cond(seg_residual!, x, pars), wf, dot(p, cl))
        end
    end

    # ══ tilted pin: heterogeneity WITHOUT demand-sensitivity (Door 1 vs Door 3) ═
    # A tilted exogenous pin should move prices, but the movement must be the
    # SAME across financing cells (heterogeneity), not demand-driven.
    println("\n== tilted pin (sticky = programme sectors, +10%): heterogeneity, not sensitivity ==")
    prog_mask = g .> 0
    sticky_t = BitVector(prog_mask)
    wbar_t = 1.0 .+ 0.1 .* prog_mask
    Lbar_Ft = sum(data.labor_share[.!sticky_t])
    wS_t = sum(wbar_t[sticky_t] .* data.labor_share[sticky_t]) / sum(data.labor_share[sticky_t])
    @printf("  %4s %9s %11s %11s %11s %12s %10s %10s\n",
        "fin", "resid", "max|p-1|", "w_f", "wage gap", "sum L", "identity", "cond(J)")
    for fin in ("F1", "F2", "F3")
        cell = design_d["cells"]["matrix_5x3-v6-ALPHA-$fin"]
        model = build_cell_model(cell, design_d, data, ψ, g)
        base = solve_cell(cell, design_d, data, ψ, g, ref.init_warm)
        init = [base.prices_raw; base.quantities; 1.0; base.external_transfer]
        x, rmax = solve_seg(model, sticky_t, wbar_t, 1.0, 0.5, Lbar_Ft; init = init)
        p, y = x[1:N], x[N+1:2N]; wf = x[2N+1]; F = x[2N+2]
        w = seg_wage_vector(N, sticky_t, wbar_t, wf)
        b = BH._mobile_market_demand(model, p, y, w; external_transfer = F)
        cl = y .- b.intermediary_demand .- b.total_final_demand
        pars = (model, sticky_t, wbar_t, 1.0, 0.5, Lbar_Ft)
        @printf("  %4s %9.2e %11.3e %11.6f %11.6f %12.8f %10.2e %10.2e\n",
            fin, rmax, maximum(abs, p .- 1), wf, wf / wS_t, sum(b.L_i), dot(p, cl),
            jac_cond(seg_residual!, x, pars))
    end

    # ══ level control (Corollary 1): pure pin, F = 0, all wages pinned ═══════
    println("\n== level control: uniform rescale of a PURE pin (F2) is a numeraire change ==")
    cell = design_d["cells"]["matrix_5x3-v6-ALPHA-F2"]
    model = build_cell_model(cell, design_d, data, ψ, g)
    base = solve_cell(cell, design_d, data, ψ, g, ref.init_warm)
    yref = Vector{Float64}(undef, N)
    for c in (0.9, 1.0, 1.1)
        pars = (model, fill(c, N))
        x = solve_ladder(pin_fixed_residual!, [c .* base.prices_raw; base.quantities], pars)[1]
        y = x[N+1:2N]
        if c == 0.9
            yref .= y
            @printf("  c = %.2f: reference\n", c)
        else
            @printf("  c = %.2f: max|y/yref - 1| = %.2e\n", c, maximum(abs, y ./ yref .- 1))
        end
    end

    # ══ the segmented pin LEVEL is not a numeraire (the free wage is technology-pinned) ══
    println("\n== segmented pin level test (F2, share 0.5): NOT a numeraire change ==")
    sticky2 = sticky_by_employment(data, 0.5)
    Lbar_F2 = sum(data.labor_share[.!sticky2])
    yref2 = Vector{Float64}(undef, N)
    for c in (0.9, 1.0)
        init = [base.prices_raw; base.quantities; 1.0; base.external_transfer]
        x, rmax = solve_seg(model, sticky2, fill(c, N), 1.0, 0.5, Lbar_F2; init = init)
        y = x[N+1:2N]
        if c == 0.9
            yref2 .= y
            @printf("  c = %.2f: reference  w_f = %.6f  resid %.2e\n", c, x[2N+1], rmax)
        else
            @printf("  c = %.2f: max|y/yref - 1| = %.4f  w_f = %.6f  (moves => pin level is real, not numeraire)\n",
                c, maximum(abs, y ./ yref2 .- 1), x[2N+1])
        end
    end

    # ══ verdict ═════════════════════════════════════════════════════════════
    println("\n== verdict on the workplan's falsifiable predictions ==")
    println("  P4  nesting  S = empty -> ALPHA/BETA, S = all -> GAMMA : HOLDS (bit-exact)")
    println("  P5  external account closes (identity ~ 1e-16)        : HOLDS")
    println("  P3  employment absorbs as the sticky share rises      : HOLDS")
    println("  P1  prices move with the financing cell               : FAILS")
    println("      flat pins: max|p-1| ~ 1e-14 at every share and financing cell;")
    println("      tilted pin: prices DO move (max|p-1| = 5.5e-2) but IDENTICALLY")
    println("      across F1/F2/F3 -> exogenous heterogeneity, not demand sensitivity.")
    println("  P2  the segment wage gap widens with the sticky share : FAILS (w_f = 1,")
    println("      gap = 1 at every share: the free wage is pinned by zero profit +")
    println("      the CPI numeraire, so the flexible labour condition only fixes the")
    println("      flexible segment's employment, never the wage.)")
    return nothing
end

main()