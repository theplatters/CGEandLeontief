# experiments/welfare-ownweights.jl — the own-weights bundle index, all 15 cells
#
# The own-weights aggregate evaluates each state's bundle under that state's
# OWN composition weights. By the identity U(q, q/sum(q)) = sum(q) it reduces
# to the ratio of quantity sums -- which coincides with the bundle's VALUE
# ratio only where prices are at their base level. BF and BETA move prices, so
# both variants are printed to show where the specification has to choose.
#
#   julia --project=. experiments/welfare-ownweights.jl

using OpenCGE, Printf

const AbstractPricing = OpenCGE.AbstractPricingClosure
isdefined(Main, :load_design) || include(joinpath(@__DIR__, "run.jl"))

const SIGMA = 0.9

function U(q, w)
    w = w ./ sum(w)
    return (sum(w .^ (1 / SIGMA) .* q .^ ((SIGMA - 1) / SIGMA)))^(SIGMA / (SIGMA - 1))
end

function main()
    root = default_root()
    dd = load_design("matrix_5x3_v10"; root = root)
    ref = build_reference(dd; root = root)
    drops = Vector{Int}(dd["data"]["drops"])
    kept = sort(setdiff(1:(length(ref.data.factor_share) + length(drops)), drops))
    ψ, g = programme_vectors(dd, ref.data, kept; root = root)

    bundle(sol; lam = 1.0) = final_demand_bundle(sol.model, sol.prices_raw,
        sol.quantities, OpenCGE._bundle_wage(sol); lambda = lam,
        external_transfer = sol.external_transfer)
    b0 = bundle(ref.sol)

    @printf("%-4s %-6s %10s %10s %10s %12s\n", "fin", "labour",
        "ownweights", "value ratio", "CES base", "max|p-1|")
    for fin in ["F1", "F2", "F3"], lab in ["BF", "ALPHA", "BETA", "GAMMA", "DELTA"]
        id = "matrix_5x3-v10-$lab-$fin" * (lab == "BETA" ? "-etas05" : "")
        sol = solve_cell(design_cell(dd, id), dd, ref.data, ψ, g, ref.init_warm)
        b = bundle(sol)
        own = sum(b.f) / sum(b0.f) - 1                       # own weights == sum(q)
        vlr = sum(sol.prices .* b.f) / sum(ref.sol.prices .* b0.f) - 1
        ces = real_final_demand_ces(sol, ref.sol; lambda = 1.0, sigma = SIGMA) - 1
        @printf("%-4s %-6s %9.4f%% %9.4f%% %9.4f%% %12.2e\n", fin, lab,
            100 * own, 100 * vlr, 100 * ces, maximum(abs.(sol.prices_raw .- 1)))
    end
    return 0
end

exit(main())
