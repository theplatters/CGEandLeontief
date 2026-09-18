# experiments/probes/probe5_bf_f3_external_position.jl
#
# Probe: is the reported "net external position" of the BF-F3 cell a result?
#
# Question (user, 2026-09-18): the v5 flow table reports BF-F3 net external
# position = +0.01330991 (v4 reported -0.007947). Recheck the path.
#
# What this probe does, for the matrix_5x3_v5 design:
#   1. reproduces the executed cells through the SAME harness functions
#      (build_reference / programme_vectors / build_cell_model / solve_cell /
#      evaluate_gates) and diffs the metrics against the committed manifests;
#   2. decomposes the external account of every cell: the booked position
#      F + B_gov vs the resource-side imbalance S + T + M - (I+X), and the
#      identity gap between them;
#   3. checks the analytic claim that at eta = 0 the gap is exactly the
#      factor-market gap, diff = -w*(sum L^cm - Lbar);
#   4. solves the eta = 0 endpoint under an ALTERNATIVE closure that replaces
#      the F = 0 pin with the labour-market equation (sum L^cm = Lbar), which
#      is the equation the pin displaces, and reports the resulting external
#      position + whether F2/F3 financing neutrality extends to eta = 0.
#
# Read-only: solves, prints, writes nothing into the repo.

include(joinpath(@__DIR__, "..", "run.jl"))

using Printf, TOML, LinearAlgebra
import NonlinearSolve

const DESIGN = "matrix_5x3_v5"

# ── eta = 0 endpoint with the labour equation instead of the F = 0 pin ──────
# Identical to BeyondHulten.problem except equation 2N+1: the COST-MINIMIZING
# aggregate labour demand is required to equal the labour bar
# (sum L^cm = Lbar) instead of the pin F = 0. By the Walras derivation
# diff = -w*(sum L^cm - sum L_supplied), so this is exactly the equation that
# closes the external account at eta = 0: the frozen allocation is priced at
# the marginal-product wage instead of carrying an unbooked labour wedge.
# F is then solved rather than pinned.
function problem_bf_labour(out::Vector, X::Vector, model)
    (; data, options) = model
    N = length(data.factor_share)
    p = BeyondHulten._positive_floor(X[1:N])
    y = BeyondHulten._positive_floor(X[N+1:2N])
    w = max(X[2N+1], 1e-10)
    F = X[2N+2]
    blocks = BeyondHulten._mobile_market_demand(model, p, y, w; external_transfer = F)
    cpi = sum(data.consumption_share .* p .^ (1 - options.elasticities.σ))^(1 / (1 - options.elasticities.σ))
    out[1:N] .= p .- blocks.cost
    out[N+1:2N] .= y .- blocks.intermediary_demand .- blocks.total_final_demand
    out[2N+1] = sum(BeyondHulten._cost_minimizing_labor(p, y, w, model)) - options.labor_bar
    out[2N+2] = cpi - 1.0
    nothing
end

"Newton + residual-gated LM polish (the ADR-0015 ladder, as in solve())."
function solve_custom(problem_f, model, init::Vector{Float64})
    x = init
    rmax = maximum(abs, begin
        o = similar(x); problem_f(o, x, model); o
    end)
    if rmax > 1e-12
        res = NonlinearSolve.solve(
            NonlinearSolve.NonlinearProblem(problem_f, x, model),
            reltol = 1e-6, abstol = 1e-6, maxiters = 20000)
        x = res.u
        rmax = maximum(abs, begin
            o = similar(x); problem_f(o, x, model); o
        end)
        for _ in 1:4
            rmax <= 1e-10 && break
            res = NonlinearSolve.solve(
                NonlinearSolve.NonlinearProblem(problem_f, x, model),
                NonlinearSolve.LevenbergMarquardt(); reltol = 1e-12, abstol = 1e-12, maxiters = 20000)
            x_new = res.u
            r_new = maximum(abs, begin
                o = similar(x_new); problem_f(o, x_new, model); o
            end)
            r_new < rmax || break
            x, rmax = x_new, r_new
        end
    end
    x, rmax
end

"External-account decomposition of a solved cell at its canonical vector X."
function account(model, X)
    N = length(model.data.factor_share)
    p, y, w, F = X[1:N], X[N+1:2N], X[2N+1], X[2N+2]
    can = BeyondHulten.external_balance_canary(model, X)
    # The COST-MINIMIZING labour demand (the zero-profit FOC), not
    # `sectoral_labor_demand`, which at eta = 0 returns the frozen baseline
    # allocation and would therefore report the gap as zero by construction.
    Lcm = sum(BeyondHulten._cost_minimizing_labor(p, y, w, model))
    # Labour actually supplied: the frozen baseline bar at eta = 0, the
    # demand-determined allocation in every eta = 1 regime (mobile and
    # fixed-wage, where the FOC sets employment).
    Lsup = model.options.elasticities.η == 0.0 ? sum(model.data.labor_share) : Lcm
    Lbar = model.options.labor_bar
    resource = can.S + can.T + can.M - can.IX      # S + T_int + M - (I+X)
    booked = can.financing                          # F + B_gov
    (S = can.S, T = can.T, M = can.M, IX = can.IX,
        F = can.external_transfer, Bgov = can.programme_financing,
        booked = booked, resource = resource, gap = can.diff,
        gap_from_labour = -w * (Lcm - Lsup), Lcm = Lcm, Lsup = Lsup, Lbar = Lbar,
        resid = maximum(abs, equilibrium_residuals(model, X)))
end

function main()
    root = default_root()
    design_d = load_design(DESIGN; root = root)
    @info "building the reference continuation (warm starts + calibration)"
    ref = build_reference(design_d; root = root)
    data = ref.data
    N = length(data.factor_share)
    ψ, g = programme_vectors(design_d, data, collect(1:N); root = root)
    @printf("N = %d   s = %.6f   sum(g) = %.10f   p·g(baseline) = %.10f\n",
        N, data.saving_rate, sum(g), sum(g))

    order = cell_order(DESIGN; root = root)

    println("\n== 1. reproduction against the committed manifests ==")
    @printf("%-12s %-12s %14s %14s %14s\n", "cell", "metric", "reproduced", "manifest", "diff")
    for id in order
        cell = design_d["cells"][id]
        model = build_cell_model(cell, design_d, data, ψ, g)
        sol = solve_cell(cell, design_d, data, ψ, g, ref.init_warm)
        ev = evaluate_gates(cell, design_d, model, sol, ref.sol)
        man = TOML.parsefile(joinpath(root, "runs", id, "manifest.toml"))
        for key in ("external_transfer", "programme_financing", "external_position",
            "gdp_wedge", "employment")
            got = ev.metrics[key]
            want = man["metrics"][key]
            @printf("%-12s %-12s %14.10g %14.10g %14.2e\n", id, key, got, want, abs(got - want))
        end
        for key in ("canary_diff",)
            got = ev.diagnostics[key]
            want = man["diagnostics"][key]
            @printf("%-12s %-12s %14.10g %14.10g %14.2e\n", id, key, got, want, abs(got - want))
        end
    end

    println("\n== 2. external-account decomposition (all 15 cells) ==")
    println("resource = S + T_int + M - (I+X);  booked = F + B_gov;  gap = resource - booked")
    @printf("%-12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
        "cell", "S", "T_int", "M", "I+X", "F", "B_gov", "resource", "booked")
    rows = Dict{String,Any}()
    for id in order
        cell = design_d["cells"][id]
        model = build_cell_model(cell, design_d, data, ψ, g)
        sol = solve_cell(cell, design_d, data, ψ, g, ref.init_warm)
        X = [sol.prices_raw; sol.quantities; sol.wages_raw[1]; sol.external_transfer]
        a = account(model, X)
        rows[id] = a
        @printf("%-12s %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n",
            id, a.S, a.T, a.M, a.IX, a.F, a.Bgov, a.resource, a.booked)
    end

    println("\n== 3. gap vs the factor-market gap  (diff = -w·(sum L^cm - sum L_supplied)) ==")
    @printf("%-12s %14s %14s %14s %14s %14s\n",
        "cell", "gap", "gap_from_labour", "sum L^cm", "sum L_supplied", "resid")
    for id in order
        a = rows[id]
        @printf("%-12s %14.10g %14.10g %14.10f %14.10f %14.2e\n",
            id, a.gap, a.gap_from_labour, a.Lcm, a.Lsup, a.resid)
    end

    println("\n== 4. eta = 0 with the labour equation instead of the F pin ==")
    println("   (the pin's equation is the labour market; solving it identifies F)")
    for fin_id in ("F1", "F2", "F3")
        id = "matrix_5x3-v5-BF-$fin_id"
        cell = design_d["cells"][id]
        model = build_cell_model(cell, design_d, data, ψ, g)
        init = copy(ref.init_warm)
        x, rmax = solve_custom(problem_bf_labour, model, init)
        X = x
        a = account(model, X)
        @printf("  BF-%s alt: F = %+.10f  B_gov = %+.10f  booked = %+.10f  resource = %+.10f  gap = %+.2e  resid = %.2e  L^cm = %.10f\n",
            fin_id, a.F, a.Bgov, a.booked, a.resource, a.gap, rmax, a.Lcm)
        @printf("           pinned-F=0 reference: booked = %+.10f  resource = %+.10f  gap = %+.2e  L^cm = %.10f\n",
            rows[id].booked, rows[id].resource, rows[id].gap, rows[id].Lcm)
    end

    println("\n== 5. financing neutrality at eta = 0 under the alternative closure ==")
    alt = Dict{String,NamedTuple}()
    for fin_id in ("F1", "F2", "F3")
        cell = design_d["cells"]["matrix_5x3-v5-BF-$fin_id"]
        model = build_cell_model(cell, design_d, data, ψ, g)
        x, _ = solve_custom(problem_bf_labour, model, copy(ref.init_warm))
        alt[fin_id] = (model = model, X = x, a = account(model, x))
    end
    a2, a3 = alt["F2"], alt["F3"]
    @printf("  max|p_F2 - p_F3| = %.3e   max|y_F2 - y_F3| = %.3e   |w_F2 - w_F3| = %.3e\n",
        maximum(abs, a2.X[1:71] .- a3.X[1:71]),
        maximum(abs, a2.X[72:142] .- a3.X[72:142]),
        abs(a2.X[143] - a3.X[143]))
    @printf("  F_F2 = %+.10f   F_F3 = %+.10f   F_F3 - (F_F2 - B_gov) = %+.2e\n",
        a2.a.F, a3.a.F, a3.a.F - (a2.a.F - a3.a.Bgov))
    @printf("  booked_F2 = %+.10f   booked_F3 = %+.10f   difference = %+.2e\n",
        a2.a.booked, a3.a.booked, a2.a.booked - a3.a.booked)
    @printf("  resource_F2 = %+.10f  resource_F3 = %+.10f  difference = %+.2e\n",
        a2.a.resource, a3.a.resource, a2.a.resource - a3.a.resource)
    return nothing
end

main()
