# experiments/welfare-f1-decomp.jl — why W_a and W_b differ under F1
#
# Under F1 the programme is NOT a public purchase (`additive_demand` returns
# zeros), so `g` is unchanged and the programme lives entirely inside the
# household's consumption via the preference tilt d_i. The question is then why
# the consumption index (W_a) and the c+g bundle index (W_b) disagree.
#
# Two candidate causes, separated here:
#   (i)  the functional form on the SAME consumption basket:
#        CES(lambda = 0) vs W_a' (exact) vs W_a (Tornqvist);
#   (ii) adding g: the lambda ladder from 0 to 1.
#
#   julia --project=. experiments/welfare-f1-decomp.jl

using OpenCGE, TOML, Printf

const AbstractPricing = OpenCGE.AbstractPricingClosure   # harness/kernel name drift

isdefined(Main, :load_design) || include(joinpath(@__DIR__, "run.jl"))

const SIGMA = 0.9

function main()
    root = default_root()
    design = "matrix_5x3_v10"
    dd = load_design(design; root = root)
    ref = build_reference(dd; root = root)
    drops = Vector{Int}(dd["data"]["drops"])
    n_full = length(ref.data.factor_share) + length(drops)
    kept = sort(setdiff(1:n_full, drops))
    ψ, g = programme_vectors(dd, ref.data, kept; root = root)

    # baseline bundle sizes, for the dilution weights
    b0 = final_demand_bundle(ref.sol.model, ref.sol.prices_raw, ref.sol.quantities,
        OpenCGE._bundle_wage(ref.sol))
    println(@sprintf("baseline: C = %.4f  G = %.4f   C/(C+G) = %.4f", sum(b0.c), sum(b0.g),
        sum(b0.c) / sum(b0.c .+ b0.g)))

    for fin in ["F1", "F2", "F3"], lab in ["BF", "ALPHA", "BETA", "GAMMA", "DELTA"]
        id = "matrix_5x3-v10-$lab-$fin" * (lab == "BETA" ? "-etas05" : "")
        cell = design_cell(dd, id)
        sol = solve_cell(cell, dd, ref.data, ψ, g, ref.init_warm)
        wa = 100 * (real_consumption(sol) / real_consumption(ref.sol) - 1)
        wap = 100 * (real_consumption_ev(sol, ref.sol) - 1)
        ces(λ) = 100 * (real_final_demand_ces(sol, ref.sol; lambda = λ, sigma = SIGMA) - 1)
        bs = final_demand_bundle(sol.model, sol.prices_raw, sol.quantities,
            OpenCGE._bundle_wage(sol); external_transfer = sol.external_transfer)
        gs = 100 * (sum(bs.g) / sum(b0.g) - 1)
        @printf("%-3s %-6s  W_a %+7.3f  W_a' %+7.3f | CES c-only %+7.3f | lam=.25 %+7.3f  .5 %+7.3f | W_b %+7.3f | dsum(g) %+.1e\n",
            fin, lab, wa, wap, ces(0.0), ces(0.25), ces(0.5), ces(1.0), gs)
    end
    return 0
end

exit(main())
