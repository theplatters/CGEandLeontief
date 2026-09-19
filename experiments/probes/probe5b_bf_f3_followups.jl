# experiments/probes/probe5b_bf_f3_followups.jl
#
# Follow-ups to probe5 (same session, same question: the BF-F3 external position).
#
# (A) Identification: the eta = 0 endpoint pins F = 0. Replace the pin by
#     F = c for a range of c and solve. If every c gives an exact root with a
#     different resource imbalance, then at eta = 0 the external position is
#     whatever the pin says: F is unidentified and the reported
#     "net external position" = c + B_gov is not an equilibrium object.
#
# (B) The closing closure: at eta = 0 replace the pin by the labour equation
#     sum L^cm = Lbar (the equation the pin displaces). Then compare the
#     solution with the corresponding ALPHA cell at full precision.

include(joinpath(@__DIR__, "..", "run.jl"))

using Printf, TOML
import NonlinearSolve

const DESIGN = "matrix_5x3_v5"

"eta = 0 system with the pin replaced by F = c."
function problem_bf_pin(out::Vector, X::Vector, model, c::Float64)
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
    out[2N+1] = F - c
    out[2N+2] = cpi - 1.0
    nothing
end

"eta = 0 system with the labour equation (sum L^cm = Lbar) instead of the pin."
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

function solve_f(problem_f, model, init::Vector{Float64}; c = nothing)
    f = c === nothing ? problem_f : (o, x, m) -> problem_f(o, x, m, c)
    x = copy(init)
    res = NonlinearSolve.solve(
        NonlinearSolve.NonlinearProblem(f, x, model),
        reltol = 1e-8, abstol = 1e-8, maxiters = 20000)
    x = res.u
    rmax = maximum(abs, begin
        o = similar(x); f(o, x, model); o
    end)
    for _ in 1:4
        rmax <= 1e-12 && break
        res = NonlinearSolve.solve(
            NonlinearSolve.NonlinearProblem(f, x, model),
            NonlinearSolve.LevenbergMarquardt(); reltol = 1e-13, abstol = 1e-13, maxiters = 20000)
        x_new = res.u
        r_new = maximum(abs, begin
            o = similar(x_new); f(o, x_new, model); o
        end)
        r_new < rmax || break
        x, rmax = x_new, r_new
    end
    x, rmax
end

function main()
    root = default_root()
    design_d = load_design(DESIGN; root = root)
    ref = build_reference(design_d; root = root)
    data = ref.data
    N = length(data.factor_share)
    ψ, g = programme_vectors(design_d, data, collect(1:N); root = root)
    init0 = [ones(N); data.λ; 1.0; 0.0]

    println("== A. eta = 0, pin F = c for a range of c (BF-F3) ==")
    println("   if every c is an exact root with a different account, F is unidentified")
    @printf("%12s %12s %14s %14s %14s %14s\n",
        "c", "resid", "F", "B_gov", "resource", "booked=c+B_gov")
    cell = design_d["cells"]["matrix_5x3-v5-BF-F3"]
    model = build_cell_model(cell, design_d, data, ψ, g)
    for c in (-0.02, -0.01, 0.0, 0.01, 0.02)
        x, rmax = solve_f(problem_bf_pin, model, init0; c = c)
        can = BeyondHulten.external_balance_canary(model, x)
        res = can.S + can.T + can.M - can.IX
        @printf("%12.4f %12.2e %14.10f %14.10f %14.10f %14.10f\n",
            c, rmax, can.external_transfer, can.programme_financing, res, can.financing)
    end
    # the same sweep under F2 (no B_gov): the account follows the pin 1:1 too
    println("   same sweep on BF-F2 (B_gov = 0):")
    cell2 = design_d["cells"]["matrix_5x3-v5-BF-F2"]
    model2 = build_cell_model(cell2, design_d, data, ψ, g)
    for c in (-0.01, 0.0, 0.01)
        x, rmax = solve_f(problem_bf_pin, model2, init0; c = c)
        can = BeyondHulten.external_balance_canary(model2, x)
        res = can.S + can.T + can.M - can.IX
        @printf("%12.4f %12.2e %14.10f %14.10f %14.10f %14.10f\n",
            c, rmax, can.external_transfer, can.programme_financing, res, can.financing)
    end

    println("\n== B. eta = 0 with the labour equation vs the ALPHA cell (full precision) ==")
    for fin_id in ("F1", "F2", "F3")
        bf = build_cell_model(design_d["cells"]["matrix_5x3-v5-BF-$fin_id"], design_d, data, ψ, g)
        al = build_cell_model(design_d["cells"]["matrix_5x3-v5-ALPHA-$fin_id"], design_d, data, ψ, g)
        xa, _ = solve_f(problem_bf_labour, bf, copy(ref.init_warm))
        sola = solve(al; init = ref.init_warm)
        Xa = [sola.prices_raw; sola.quantities; sola.wages_raw[1]; sola.external_transfer]
        canb = BeyondHulten.external_balance_canary(bf, xa)
        cana = BeyondHulten.external_balance_canary(al, Xa)
        Lfrozen = bf.data.labor_share
        Lcm = BeyondHulten._cost_minimizing_labor(xa[1:N], xa[N+1:2N], xa[2N+1], bf)
        @printf("  %s\n", fin_id)
        @printf("    max|p_BF - p_ALPHA| = %.3e   max|y_BF - y_ALPHA| = %.3e\n",
            maximum(abs, xa[1:N] .- Xa[1:N]), maximum(abs, xa[N+1:2N] .- Xa[N+1:2N]))
        @printf("    |w_BF - w_ALPHA| = %.3e     |F_BF - F_ALPHA| = %.3e\n",
            abs(xa[2N+1] - Xa[2N+1]), abs(xa[2N+2] - Xa[2N+2]))
        @printf("    F_BF = %+.12f   F_ALPHA = %+.12f   B_gov = %+.12f\n",
            canb.external_transfer, cana.external_transfer, canb.programme_financing)
        @printf("    booked_BF = %+.12f   booked_ALPHA = %+.12f   gap_BF = %.2e   gap_ALPHA = %.2e\n",
            canb.financing, cana.financing, canb.diff, cana.diff)
        @printf("    sectoral allocation still frozen: max|L_share - L^cm| = %.4f, sum L^cm = %.12f\n",
            maximum(abs, Lfrozen .- Lcm), sum(Lcm))
    end
    return nothing
end

main()
