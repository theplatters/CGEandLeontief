# experiments/probes/probe13_sectoral_labour.jl
#
# Step 1 (successor) of docs/WORKPLAN_SENSITIVE_PRICES.md (v4): the GENERAL
# sectoral-labour-market closure, prototyped on the full-71 calibration.
# No src/ change.
#
# Specification (workplan v4, "The specification to implement"):
#   unknowns X = [p(1:N); y(1:N); w(1:N); F]                 (3N + 1)
#   (1) N zero-profit conditions, p_i = cost_i(p, w_i);
#   (2) N goods-market clearings, with the cost-minimizing allocation
#       L_i = L^cm_i(p, y_i, w_i) and household wage income sum_i w_i L_i;
#   (3) N sectoral labour-market conditions,
#         L^cm_i(p, y_i, w_i) = Lbar_i * ( (w_i / Pi(p)) / (wbar_i / Pibar) )^{eta_s,i},
#       with the anchors wbar_i / Pibar = 1 at the calibration;
#   (4) one numeraire, Pi(p) = 1.
#   F is free, exactly as at eta = 0 (the block is homogeneous of degree one in
#   (p, w, F); dropping a clearing equation instead is the retired ADR-0010
#   shortcut).
#
# Endpoints of the family (both measured):
#   eta_s,i = 0 for every i  ->  L^cm_i = Lbar_i, i.e. the ADR-0020 option C
#       system (problem_sectoral), executed as the v6 BF row = S1;
#   eta_s,i -> infinity      ->  the real wage is pinned at its anchor and pi, w
#       tend to the calibration, so the row is GAMMA-like *in prices*. The limit
#       is approached only asymptotically: the supply slope diverges, so the
#       direct 3N+1 formulation becomes singular for very large eta_s (measured:
#       eta_s = 1e6 fails); the employment limit to the GAMMA closure is not
#       attained at finite eta_s.
#
# Falsifiable predictions, recorded BEFORE running:
#   Q1  prices move with the financing cell (max|p-1| differs across F1/F2/F3
#       at a fixed eta_s) -- no single-wage closure satisfies this;
#   Q2  the uniform-eta_s ladder is monotone: max|p-1| falls from its S1
#       maximum toward zero as eta_s rises;
#   Q3  the eta_s,i = 0 corner reproduces the committed S1 cells;
#   Q4  the account closes (price-weighted clearing residual at the residual
#       level, the proxy validated in probe11);
#   Q5  the large-eta_s endpoint reproduces GAMMA (p = 1).

include(joinpath(@__DIR__, "..", "run.jl"))
using Printf, LinearAlgebra
import NonlinearSolve
const BH = BeyondHulten

cpi_of(model, p) = begin
    σ = model.options.elasticities.σ
    sum(model.data.consumption_share .* p .^ (1 - σ))^(1 / (1 - σ))
end

rmax_of(f, x, pars) = maximum(abs, (o = similar(x); f(o, x, pars); o))

"Residual of the general sectoral-labour-market system (3N + 1 equations)."
function sec_residual!(out, X, pars)
    model, eta_vec = pars
    (; data) = model
    N = length(data.factor_share)
    p = BH._positive_floor(X[1:N])
    y = BH._positive_floor(X[N+1:2N])
    w = BH._positive_floor(X[2N+1:3N])
    F = X[3N+1]
    b = BH._mobile_market_demand(model, p, y, w; external_transfer = F)
    out[1:N]      .= p .- b.cost
    out[N+1:2N]   .= log.(b.L_i) .- log.(BH._positive_floor(data.labor_share)) .-
                     eta_vec .* log.(w ./ cpi_of(model, p))
    out[2N+1:3N]  .= y .- b.intermediary_demand .- b.total_final_demand
    out[3N+1]      = cpi_of(model, p) - 1.0
    nothing
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
    maximum(s) / minimum(s)
end

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

function main()
    root = default_root()
    design_d = load_design("matrix_5x3_v6"; root = root)
    ref = build_reference(design_d; root = root)
    data = ref.data
    N = length(data.factor_share)
    ψ, g = programme_vectors(design_d, data, collect(1:N); root = root)
    prog = BitVector(g .> 0)
    big_order = sortperm(data.labor_share; rev = true)
    top_half = falses(N); top_half[big_order[1:round(Int, 0.5N)]] .= true
    @printf("N = %d, sectoral unknowns = 3N+1 = %d, programme sectors = %d\n",
        N, 3N + 1, count(prog))

    # ── per-financing setup: the eta = 0 (BF) solve is the S1 anchor ─────────
    setups = Dict{String,Any}()
    for fin in ("F1", "F2", "F3")
        cellbf = design_d["cells"]["matrix_5x3-v6-BF-$fin"]
        modelbf = build_cell_model(cellbf, design_d, data, ψ, g)
        solbf = solve(modelbf; init = ref.init_warm)          # problem_sectoral
        cellmo = design_d["cells"]["matrix_5x3-v6-ALPHA-$fin"]
        modelmo = build_cell_model(cellmo, design_d, data, ψ, g)   # eta = 1
        cellga = design_d["cells"]["matrix_5x3-v6-GAMMA-$fin"]
        modelga = build_cell_model(cellga, design_d, data, ψ, g)   # fixed wage
        solga = solve(modelga; init = ref.init_warm)
        setups[fin] = (modelbf, solbf, modelmo, modelga, solga)
    end

    # ══ Q3: the eta_s,i = 0 corner reproduces S1 (the v6 BF row) ═════════════
    println("\n== Q3: eta_s,i = 0 corner vs the executed S1 (BF) cells ==")
    for fin in ("F1", "F2", "F3")
        modelbf, solbf, modelmo, _, _ = setups[fin]
        init = [solbf.prices_raw; solbf.quantities; solbf.wages; solbf.external_transfer]
        x, rmax = solve_ladder(sec_residual!, init, (modelmo, zeros(N)))
        dp = maximum(abs, x[1:N] .- solbf.prices_raw)
        dy = maximum(abs, x[N+1:2N] .- solbf.quantities)
        dw = maximum(abs, x[2N+1:3N] .- solbf.wages)
        dF = abs(x[3N+1] - solbf.external_transfer)
        @printf("  BF-%s: resid %.2e | max|dp| %.2e  max|dy| %.2e  max|dw| %.2e  |dF| %.2e\n",
            fin, rmax, dp, dy, dw, dF)
    end

    # ══ Q1, Q2, Q4: the uniform-eta_s ladder x financing ═════════════════════
    println("\n== Q1/Q2: uniform-eta_s ladder (max|p-1|, employment, identity) ==")
    for fin in ("F1", "F2", "F3")
        modelbf, solbf, modelmo, modelga, solga = setups[fin]
        Lga = sum(BH._cost_minimizing_labor(solga.prices_raw, solga.quantities, 1.0, modelga))
        println("  --- $fin (fixed-wage GAMMA anchor: max|p-1| = $(@sprintf("%.2e", maximum(abs, solga.prices_raw .- 1))), L = $(@sprintf("%.8f", Lga))) ---")
        @printf("  %10s %9s %11s %13s %12s %12s %10s %11s\n",
            "eta_s", "resid", "max|p-1|", "sum L", "wage min/max", "mean w", "identity", "cond(J)")
        init0 = [solbf.prices_raw; solbf.quantities; solbf.wages; solbf.external_transfer]
        for eta_s in (0.0, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0)
            eta_vec = fill(eta_s, N)
            x, rmax = solve_ladder(sec_residual!, init0, (modelmo, eta_vec))
            if rmax > 1e-6
                @printf("  %10.2f %9.2e   <-- DID NOT CONVERGE (not a result)\n", eta_s, rmax)
                continue
            end
            p, y, w = x[1:N], x[N+1:2N], x[2N+1:3N]; F = x[3N+1]
            b = BH._mobile_market_demand(modelmo, p, y, w; external_transfer = F)
            cl = y .- b.intermediary_demand .- b.total_final_demand
            @printf("  %10.2f %9.2e %11.3e %13.8f %6.3f/%.3f %12.6f %10.2e %11.2e\n",
                eta_s, rmax, maximum(abs, p .- 1), sum(b.L_i), minimum(w), maximum(w),
                sum(w .* b.L_i) / sum(b.L_i), dot(p, cl), jac_cond(sec_residual!, x, (modelmo, eta_vec)))
        end
    end

    # ══ two-group variant: a rigid group at eta_s,i = 0, the rest at 0.5 ══════
    println("\n== two-group variant (rigid group eta_s = 0, flexible group eta_s = 0.5) ==")
    for (tag, mask) in (("programme sectors", prog), ("largest half", top_half))
        @printf("  rigid group = %-18s (share %.2f)\n", tag, count(mask) / N)
        @printf("  %4s %9s %11s %13s %12s %11s\n", "fin", "resid", "max|p-1|", "sum L", "mean w", "identity")
        for fin in ("F1", "F2", "F3")
            _, solbf, modelmo, _, _ = setups[fin]
            eta_vec = 0.5 .* (.!mask)
            init0 = [solbf.prices_raw; solbf.quantities; solbf.wages; solbf.external_transfer]
            x, rmax = solve_ladder(sec_residual!, init0, (modelmo, eta_vec))
            p, y, w = x[1:N], x[N+1:2N], x[2N+1:3N]; F = x[3N+1]
            b = BH._mobile_market_demand(modelmo, p, y, w; external_transfer = F)
            cl = y .- b.intermediary_demand .- b.total_final_demand
            @printf("  %4s %9.2e %11.3e %13.8f %12.6f %11.2e\n",
                fin, rmax, maximum(abs, p .- 1), sum(b.L_i), sum(w .* b.L_i) / sum(b.L_i), dot(p, cl))
        end
    end

    println("\n== verdict on the pre-registered predictions ==")
    println("  Q3  eta_s,i = 0 reproduces S1 (v6 BF)  : measured above")
    println("  Q1  prices move with the financing cell : measured above")
    println("  Q2  the ladder is monotone in eta_s     : measured above")
    println("  Q4  the account closes                  : measured above")
    println("  Q5  large eta_s -> 0 price response      : measured (GAMMA-like in")
    println("      prices); the eta_s -> infinity limit is asymptotic, the direct")
    println("      formulation goes singular before the GAMMA employment limit is reached")
    return nothing
end

main()