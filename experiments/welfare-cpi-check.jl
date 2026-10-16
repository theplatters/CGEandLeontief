# experiments/welfare-cpi-check.jl — is Pi(p) == 1 on every solve path?
#
# `real_consumption_ev` (W_a') is defined as the ratio of consumption
# expenditure `sum(p .* c)`. That is the EXACT money-metric only if the
# consumption price index is the numeraire, i.e. Pi(p) = 1 in the cell, since
# the true money-metric ratio is (E_1/E_0)/Pi(p_1). The BF and BETA rows move
# prices, so the claim has to be checked rather than assumed.
#
#   julia --project=. experiments/welfare-cpi-check.jl

using OpenCGE, Printf

const AbstractPricing = OpenCGE.AbstractPricingClosure

isdefined(Main, :load_design) || include(joinpath(@__DIR__, "run.jl"))

function main()
    root = default_root()
    dd = load_design("matrix_5x3_v10"; root = root)
    ref = build_reference(dd; root = root)
    drops = Vector{Int}(dd["data"]["drops"])
    kept = sort(setdiff(1:(length(ref.data.factor_share) + length(drops)), drops))
    ψ, g = programme_vectors(dd, ref.data, kept; root = root)
    @printf("reference: cpi = %.16f\n", cpi(ref.sol))
    for fin in ["F1", "F2", "F3"], lab in ["BF", "ALPHA", "BETA", "GAMMA", "DELTA"]
        id = "matrix_5x3-v10-$lab-$fin" * (lab == "BETA" ? "-etas05" : "")
        sol = solve_cell(design_cell(dd, id), dd, ref.data, ψ, g, ref.init_warm)
        P = cpi(sol)
        er = sum(sol.prices .* final_demand_bundle(sol.model, sol.prices_raw, sol.quantities,
            OpenCGE._bundle_wage(sol); lambda = 0.0).c) /
             sum(ref.sol.prices .* final_demand_bundle(ref.sol.model, ref.sol.prices_raw,
            ref.sol.quantities, OpenCGE._bundle_wage(ref.sol); lambda = 0.0).c)
        @printf("%-3s %-6s  cpi %+.3e | W_a' raw %+8.4f | /cpi %+8.4f | max|p-1| %.2e\n",
            fin, lab, P - 1, 100 * (er - 1), 100 * (er / P - 1), maximum(abs.(sol.prices_raw .- 1)))
    end
    return 0
end

exit(main())
