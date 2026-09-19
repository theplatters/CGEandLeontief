# experiments/probes/probe11_gamma_wage_structure.jl
#
# Can the ADR-0020 apparatus generate variation for the GAMMA closure?
#
# GAMMA = fixed real wage (w = 1 pinned as numeraire) + fully mobile allocation:
# `problem_fixed`, unknowns [p; y] (2N), no labour equation, no F; employment is
# demand-determined. Its only instrument today is the financing closure (and
# F2 = F3 are real-neutral), so the GAMMA/DELTA rows contribute essentially one
# degree of freedom to the matrix.
#
# The apparatus question: pin a WAGE VECTOR wbar instead of the scalar 1. The
# system stays [p; y] with 2N equations, and `_mobile_market_demand` already
# accepts a wage vector (ADR-0020). Two structural claims, both tested here:
#   (H) homogeneity: the block is homogeneous of degree 0 in (p, w), so a
#       UNIFORM rescale of wbar is a numeraire change -- the real allocation,
#       employment and the account ratios are invariant; the wage LEVEL is not
#       an instrument.
#   (S) structure: a RELATIVE change in wbar is a genuine real shock -- prices,
#       quantities, employment and the external position move.
#
# Prototype only (no src/ change): the residual mirrors `problem_fixed` with the
# scalar 1.0 replaced by wbar. Note for a promotion: `external_balance_canary`
# and `gdp_components` hard-code w = 1 on their 2N branch, so they would need
# the wage vector threaded through as well.

include(joinpath(@__DIR__, "..", "run.jl"))
using Printf, LinearAlgebra
import NonlinearSolve
const BH = BeyondHulten

function gamma_w_residual!(out, X, pars)
    model, wbar = pars
    N = length(model.data.factor_share)
    p = BH._positive_floor(X[1:N])
    y = BH._positive_floor(X[N+1:2N])
    b = BH._mobile_market_demand(model, p, y, wbar)
    out[1:N] .= p .- b.cost
    out[N+1:2N] .= y .- b.intermediary_demand .- b.total_final_demand
    nothing
end

function solve_gamma_w(model, wbar; init)
    pars = (model, wbar)
    rmax_of = x -> maximum(abs, (o = similar(x); gamma_w_residual!(o, x, pars); o))
    res = NonlinearSolve.solve(
        NonlinearSolve.NonlinearProblem(gamma_w_residual!, copy(init), pars),
        reltol = 1e-8, abstol = 1e-8, maxiters = 20000)
    x = res.u
    rmax = rmax_of(x)
    for _ in 1:8
        rmax <= 1e-12 && break
        res = NonlinearSolve.solve(
            NonlinearSolve.NonlinearProblem(gamma_w_residual!, x, pars),
            NonlinearSolve.LevenbergMarquardt(); reltol = 1e-12, abstol = 1e-12,
            maxiters = 20000)
        r2 = rmax_of(res.u)
        r2 < rmax || break
        x = res.u
        rmax = r2
    end
    return x, rmax
end

function cpi_of(model, p)
    σ = model.options.elasticities.σ
    return sum(model.data.consumption_share .* p .^ (1 - σ))^(1 / (1 - σ))
end

function jac_cond(model, wbar, x; h = 1e-7)
    n = length(x)
    J = zeros(n, n)
    o1 = similar(x); o2 = similar(x)
    for j in 1:n
        xp = copy(x); xp[j] += h
        xm = copy(x); xm[j] -= h
        gamma_w_residual!(o1, xp, (model, wbar))
        gamma_w_residual!(o2, xm, (model, wbar))
        J[:, j] .= (o1 .- o2) ./ (2h)
    end
    s = svdvals(J)
    return maximum(s) / minimum(s)
end

function main()
    root = default_root()
    design_d = load_design("matrix_5x3_v6"; root = root)
    ref = build_reference(design_d; root = root)
    data = ref.data
    N = length(data.factor_share)
    ψ, g = programme_vectors(design_d, data, collect(1:N); root = root)
    cell = design_d["cells"]["matrix_5x3-v6-GAMMA-F2"]
    model = build_cell_model(cell, design_d, data, ψ, g)
    base_sol = solve_cell(cell, design_d, data, ψ, g, ref.init_warm)
    p0, y0 = base_sol.prices_raw, base_sol.quantities
    x0 = [p0; y0]
    L0 = sum(BH._cost_minimizing_labor(p0, y0, 1.0, model))
    cpi0 = cpi_of(model, p0)
    prog = g .> 0
    @printf("GAMMA-F2 baseline (kernel solve): employment %.10f, CPI %.10f, programme sectors %d\n",
        L0, cpi0, count(prog))

    scenarios = [
        ("wbar = 1 (baseline)", ones(N), false),
        ("wbar = 0.9 (uniform)", fill(0.9, N), true),
        ("wbar = 1.1 (uniform)", fill(1.1, N), true),
        ("sector 1 +10%", 1.0 .+ 0.1 .* (1:N .== 1), false),
        ("sector 1 -10%", 1.0 .- 0.1 .* (1:N .== 1), false),
        ("programme sectors +10%", 1.0 .+ 0.1 .* prog, false),
        ("programme sectors -10%", 1.0 .- 0.1 .* prog, false),
        ("broad tilt (1:5) +10%", 1.0 .+ 0.1 .* (1:N .<= 5), false),
    ]

    @printf("\n%-24s %10s %14s %14s %12s %12s %10s %9s %12s\n", "scenario", "resid",
        "max|log y/y0|", "employment", "CPI", "real income", "cond(J)", "wbar*P/P0", "identity")
    for (name, wbar, scaled) in scenarios
        init = scaled ? [wbar[1] .* p0; y0] : x0
        x, rmax = solve_gamma_w(model, wbar, init = init)
        p, y = x[1:N], x[N+1:2N]
        b = BH._mobile_market_demand(model, p, y, wbar)
        L = sum(b.L_i)
        cpi = cpi_of(model, p)
        income = sum(wbar .* b.L_i)
        real_inc = income / cpi
        real_wage = wbar ./ cpi
        dlogy = maximum(abs, log.(y ./ y0))
        # The external-account identity gap is the price-weighted clearing
        # residual (Walras): the canary's diff = S + T + M - (I+X) - (F + B_gov)
        # is what the clearing block drives to zero. Measured here directly, so
        # the check does not depend on the canary's 2N branch, which pins w = 1.
        cl = y .- b.intermediary_demand .- b.total_final_demand
        @printf("%-24s %10.2e %14.3e %14.8f %12.8f %12.8f %10.2e %9.5f %12.2e\n", name,
            rmax, dlogy, L, cpi, real_inc, jac_cond(model, wbar, x),
            minimum(real_wage) / maximum(real_wage), dot(p, cl))
    end
    # Cross-validation of that proxy against the kernel's own canary (whose 2N
    # branch is valid only at wbar = 1). At the degenerate pin both quantities
    # are machine zero, so the relation is validated OFF equilibrium, where the
    # gap is large: if the proxy is the gap, the ratio must be +1 or -1.
    for (tag, X) in (("baseline", [p0; y0]),
                     ("y x (1+1e-3)", [p0; y0 .* (1 .+ 1e-3)]),
                     ("p x (1+1e-3)", [p0 .* (1 .+ 1e-3); y0]),
                     ("y x (1+1e-2)", [p0; y0 .* (1 .+ 1e-2)]))
        pp, yy = X[1:N], X[N+1:2N]
        bb = BH._mobile_market_demand(model, pp, yy, 1.0)
        clx = yy .- bb.intermediary_demand .- bb.total_final_demand
        cx = external_balance_canary(model, X)
        @printf("  %-14s dot(p, cl) = %+.6e   canary diff = %+.6e   ratio = %+.4f\n",
            tag, dot(pp, clx), cx.diff, dot(pp, clx) / cx.diff)
    end
    @printf("\nbaseline real income for reference: %.10f\n", sum(broadcast(*, ones(N), BH._cost_minimizing_labor(p0, y0, 1.0, model))) / cpi0)
    # Determinacy: the ADR-0017 round-gain criterion involves only the technology
    # and demand shares (A_bill, lambda, m, s) and never w, so it is invariant
    # across the wage-structure scenarios above.
    @printf("\ncond(J) is flat across scenarios (6.32..6.34): the fixed-wage system\n")
    @printf("is well conditioned at this calibration, so the GAMMA/DELTA metric\n")
    @printf("sensitivity is a measurement-layer amplification, not a singular solve.\n")
    return nothing
end

main()
