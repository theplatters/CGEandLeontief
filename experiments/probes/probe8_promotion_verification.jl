# experiments/probes/probe8_promotion_verification.jl
#
# ADR-0020 promotion check (step 2 of the plan).
#
# (i)  The twelve non-BF v5 cells must reproduce their committed manifests to
#      machine precision: the eta = 0 change is confined to `problem_sectoral`,
#      and the shared scalar-wage paths were left arithmetically unchanged.
# (ii) The three BF cells now run the sectoral-wage system through the KERNEL
#      (`solve` -> `problem_sectoral`), and must reproduce the numbers measured
#      by probe7 (F, B_gov, booked, identity gap).
# (iii) `external_balance_canary` and `gdp_components` must work on the new
#      3N+1 vector (the manifests' diagnostics).

include(joinpath(@__DIR__, "..", "run.jl"))

using Printf, TOML

const DESIGN = "matrix_5x3_v5"
const NONBF = ("ALPHA", "BETA", "GAMMA", "DELTA")

function promotion_main()
    root = default_root()
    design_d = load_design(DESIGN; root = root)
    ref = build_reference(design_d; root = root)
    data = ref.data
    N = length(data.factor_share)
    ψ, g = programme_vectors(design_d, data, collect(1:N); root = root)
    @printf("N = %d; eta = 0 canonical vector length = %d (3N+1)\n", N, 3N + 1)

    keys = ("gdp_rel", "consumption_rel", "employment", "wage", "external_transfer",
        "programme_financing", "external_position", "gdp_wedge")
    dkeys = ("canary_s", "canary_ixm", "canary_diff", "gdp_g", "gdp_m_final", "gdp_m_int",
        "gdp_t_int")

    println("\n(i) the twelve non-BF v5 cells vs their committed manifests")
    worst = 0.0
    worst_id = ""
    for labor in NONBF, fin in ("F1", "F2", "F3")
        id = "matrix_5x3-v5-$labor-$fin"
        cell = design_d["cells"][id]
        model = build_cell_model(cell, design_d, data, ψ, g)
        sol = solve_cell(cell, design_d, data, ψ, g, ref.init_warm)
        ev = evaluate_gates(cell, design_d, model, sol, ref.sol)
        man = TOML.parsefile(joinpath(root, "runs", id, "manifest.toml"))
        d = 0.0
        dk = ""
        for k in keys
            v = abs(ev.metrics[k] - man["metrics"][k])
            v > d && (d = v; dk = k)
        end
        for k in dkeys
            v = abs(ev.diagnostics[k] - man["diagnostics"][k])
            v > d && (d = v; dk = k)
        end
        @printf("  %-12s max|metric - manifest| = %.3e (%s)   resid %.2e vs %.2e   %s\n", "$labor-$fin", d,
            dk, ev.gates["residual"]["value"], man["gates"]["residual"]["value"], ev.gates["overall"])
        d > worst && (worst = d; worst_id = "$labor-$fin")
    end
    @printf("  WORST over the twelve: %.3e (%s)  -> %s\n", worst, worst_id,
        worst <= 1e-14 ? "BIT-IDENTICAL (within manifest precision)" : "DIFFERS -- investigate")

    println("\n(ii) the three BF cells through the kernel (sectoral wages, eta = 0)")
    @printf("  %-8s %14s %14s %14s %14s %14s %12s\n", "cell", "resid", "F", "B_gov",
        "booked", "gap", "cond(J)")
    for fin in ("F1", "F2", "F3")
        id = "matrix_5x3-v5-BF-$fin"
        cell = design_d["cells"][id]
        model = build_cell_model(cell, design_d, data, ψ, g)
        sol = solve_cell(cell, design_d, data, ψ, g, ref.init_warm)
        ev = evaluate_gates(cell, design_d, model, sol, ref.sol)
        X = [sol.prices_raw; sol.quantities; sol.wages_raw; sol.external_transfer]
        can = external_balance_canary(model, X)
        comp = gdp_components(model, sol)
        @printf("  %-8s %14.3e %14.8f %14.8f %14.8f %14.2e %12s\n", fin,
            maximum(abs, equilibrium_residuals(model, X)),
            sol.external_transfer, can.programme_financing, can.financing, can.diff,
            ev.gates["overall"])
        @printf("           wage vector: min %.6f max %.6f | wage_bill %.8f | gdp_components wedge %.2e\n",
            minimum(sol.wages_raw), maximum(sol.wages_raw), comp.wage_bill, comp.wedge)
        @printf("           metrics: consumption_rel %+.6f  gdp_rel %+.6f  employment %.6f\n",
            ev.metrics["consumption_rel"], ev.metrics["gdp_rel"], ev.metrics["employment"])
    end
    println("  probe7 reference: F = -0.00581497 / -0.00852841 / -0.02337469,")
    println("                    booked = -0.00581497 / -0.00852841 / -0.00852841,")
    println("                    gap <= 5e-13, B_gov(F3) = +0.01484628")

    println("\n(iii) financing neutrality at eta = 0 through the kernel")
    sols = Dict{String,Any}()
    for fin in ("F2", "F3")
        cell = design_d["cells"]["matrix_5x3-v5-BF-$fin"]
        model = build_cell_model(cell, design_d, data, ψ, g)
        sols[fin] = (sol = solve_cell(cell, design_d, data, ψ, g, ref.init_warm), model = model)
    end
    s2, s3 = sols["F2"].sol, sols["F3"].sol
    X3 = [s3.prices_raw; s3.quantities; s3.wages_raw; s3.external_transfer]
    B3 = external_balance_canary(sols["F3"].model, X3).programme_financing
    @printf("  max|p_F2-p_F3| = %.3e  max|y_F2-y_F3| = %.3e  max|w_F2-w_F3| = %.3e\n",
        maximum(abs, s2.prices_raw .- s3.prices_raw),
        maximum(abs, s2.quantities .- s3.quantities),
        maximum(abs, s2.wages_raw .- s3.wages_raw))
    @printf("  F_F3 - (F_F2 - B_gov) = %.3e   booked_F2 - booked_F3 = %.3e\n",
        s3.external_transfer - (s2.external_transfer - B3),
        s2.external_transfer - s3.external_transfer)
    return nothing
end

promotion_main()
