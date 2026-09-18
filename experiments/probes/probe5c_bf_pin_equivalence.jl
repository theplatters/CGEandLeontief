# experiments/probes/probe5c_bf_pin_equivalence.jl
#
# Third check (same question, same session): if the eta = 0 row differs from the
# eta = 1 (ALPHA) row only through the pinned external transfer, then solving
# the eta = 0 system with the pin set to ALPHA's solved F must reproduce the
# ALPHA cell. Compares the canonical vectors, the household expenditure E and
# the consumption block.

include(joinpath(@__DIR__, "..", "run.jl"))

using Printf
import NonlinearSolve

const DESIGN = "matrix_5x3_v5"

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

function solve_pin(model, init, c)
    f = (o, x, m) -> problem_bf_pin(o, x, m, c)
    x = copy(init)
    res = NonlinearSolve.solve(NonlinearSolve.NonlinearProblem(f, x, model),
        reltol = 1e-8, abstol = 1e-8, maxiters = 20000)
    x = res.u
    rmax = maximum(abs, begin
        o = similar(x); f(o, x, model); o
    end)
    for _ in 1:4
        rmax <= 1e-12 && break
        res = NonlinearSolve.solve(NonlinearSolve.NonlinearProblem(f, x, model),
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

    println("BF eta = 0 solved with the pin set to ALPHA's solved F, vs the ALPHA cell")
    @printf("%-4s %14s %14s %14s %14s\n", "fin", "max|X_BF-X_A|", "max|E_BF-E_A|",
        "max|cdom_BF-cdom_A|", "resid_BF")
    for fin_id in ("F1", "F2", "F3")
        bf = build_cell_model(design_d["cells"]["matrix_5x3-v5-BF-$fin_id"], design_d, data, ψ, g)
        al = build_cell_model(design_d["cells"]["matrix_5x3-v5-ALPHA-$fin_id"], design_d, data, ψ, g)
        sola = solve(al; init = ref.init_warm)
        XA = [sola.prices_raw; sola.quantities; sola.wages_raw[1]; sola.external_transfer]
        x, rmax = solve_pin(bf, ref.init_warm, XA[2N+2])
        blkA = BeyondHulten._mobile_market_demand(al, XA[1:N], XA[N+1:2N], XA[2N+1];
            external_transfer = XA[2N+2])
        blkB = BeyondHulten._mobile_market_demand(bf, x[1:N], x[N+1:2N], x[2N+1];
            external_transfer = x[2N+2])
        @printf("%-4s %14.3e %14.3e %14.3e %14.2e\n", fin_id,
            maximum(abs, x .- XA), abs(blkB.E - blkA.E),
            maximum(abs, blkB.c_dom .- blkA.c_dom), rmax)
        @printf("     F_BF(pin) = %+.12f   F_ALPHA = %+.12f   sum L^cm_BF = %.12f\n",
            x[2N+2], XA[2N+2], sum(BeyondHulten._cost_minimizing_labor(x[1:N], x[N+1:2N], x[2N+1], bf)))
    end
    return nothing
end

main()
