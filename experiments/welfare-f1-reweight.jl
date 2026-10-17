# experiments/welfare-f1-reweight.jl — would spelling out F1 numbers help?
#
# The proposal: print W_b and W_c in the F1 rows under "updated weightings",
# so the table itself shows they carry no information beyond W_a'.
#
# Both weight vectors are recovered from the solutions rather than from the
# tilt formula: at p = 1 the demand is c_i ∝ alpha_i p_i^(-sigma), so
# alpha_i ∝ c_i p_i^sigma is the preference weight implied by any bundle.
#
#   julia --project=. experiments/welfare-f1-reweight.jl

using OpenCGE, Printf

if !isdefined(OpenCGE, :AbstractPricing)   # the pricing-door line renames the type;
    @eval const AbstractPricing = OpenCGE.AbstractPricingClosure   # the release exports it
end
isdefined(Main, :load_design) || include(joinpath(@__DIR__, "run.jl"))

const SIGMA = 0.9

"""Direct CES utility, weights normalised to sum to one (money-metric units)."""
function U(c, wts)
    w = wts ./ sum(wts)
    return (sum(w .^ (1 / SIGMA) .* c .^ ((SIGMA - 1) / SIGMA)))^(SIGMA / (SIGMA - 1))
end

function main()
    root = default_root()
    dd = load_design("matrix_5x3_v10"; root = root)
    ref = build_reference(dd; root = root)
    drops = Vector{Int}(dd["data"]["drops"])
    kept = sort(setdiff(1:(length(ref.data.factor_share) + length(drops)), drops))
    ψ, g = programme_vectors(dd, ref.data, kept; root = root)

    b0 = final_demand_bundle(ref.sol.model, ref.sol.prices_raw, ref.sol.quantities,
        OpenCGE._bundle_wage(ref.sol))
    ω = b0.c ./ sum(b0.c)

    @printf("%-6s %9s %9s %9s %9s %9s %9s\n",
        "cell", "W_a'", "perstate", "upd/both", "base/both", "W_a", "W_b")
    for lab in ["ALPHA", "GAMMA", "DELTA"]
        sol = solve_cell(design_cell(dd, "matrix_5x3-v10-$lab-F1"), dd, ref.data, ψ, g,
            ref.init_warm)
        p = sol.prices
        bs = final_demand_bundle(sol.model, p, sol.quantities, OpenCGE._bundle_wage(sol);
            lambda = 0.0, external_transfer = sol.external_transfer)
        # restrict to the common support: c_i = 0 in many of the 62 sectors and
        # the CES needs a negative exponent to be finite
        k = (bs.c .> 0) .& (b0.c .> 0)
        c1, c0 = bs.c[k], b0.c[k]
        α = c1 .* p[k] .^ SIGMA
        α = α ./ sum(α)
        wk = ω[k] ./ sum(ω[k])
        wap = real_consumption_ev(sol, ref.sol) - 1
        # the four cells of the 2x2: each bundle valued under each weight vector
        u1a, u1w = U(c1, α), U(c1, wk)   # counterfactual bundle: own / baseline weights
        u0a, u0w = U(c0, α), U(c0, wk)   # baseline bundle:      updated / own weights
        perstate = u1a / u0w - 1         # each state by its own weights
        updboth = u1a / u0a - 1          # updated aggregator, both bundles
        baseboth = u1w / u0w - 1         # baseline aggregator, both bundles
        cross = u1w / u0a - 1            # crossed
        wa = 100 * (real_consumption(sol) / real_consumption(ref.sol) - 1)
        wb = real_final_demand_ces(sol, ref.sol; lambda = 1.0, sigma = SIGMA) - 1
        @printf("%-6s %8.4f%% %8.4f%% %8.4f%% %8.4f%% %8.4f%% %8.4f%%\n",
            lab, 100 * wap, 100 * perstate, 100 * updboth, 100 * baseboth, wa, 100 * wb)
        @printf("        U values: c1@own %.8f  c1@base %.8f  |  c0@upd %.8f  c0@own %.8f  | crossed %+8.4f%%\n",
            u1a, u1w, u0a, u0w, 100 * cross)
        @printf("        max|p - prices_raw| = %.1e ,  sum(bundle weights on g) = %.4f\n",
            maximum(abs.(sol.prices .- sol.prices_raw)), sum(b0.g) / sum(b0.f))
    end
    return 0
end

exit(main())
