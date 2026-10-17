# experiments/consumption-vs-gdp-gap.jl
#
# Why is W_C (household consumption welfare) 0.48 pp ABOVE Delta Y in the
# Leontief/fixed-wage closure under F3, when the programme's spending is public?
#
# Candidate channel: baseline government purchases g^G are tax-financed at a
# rate tau that is solved, so a higher wage bill lowers the rate needed to fund
# a fixed g^G and disposable income rises by more than the wage bill.
# Test: decompose E = (1-tau)*w*sum(L) + R.
#
#   julia --project=<OPEN_CGE_CHECKOUT> experiments/consumption-vs-gdp-gap.jl.jl
#
#   NOT `--project=.`: this script calls kernel functions added in OpenCGE
#   09e2b83 (diag(WELFARE)), which the paper repo's Manifest does not yet pin.
#   Point --project at a checkout that contains them, or the load fails.

using OpenCGE, Printf

if !isdefined(OpenCGE, :AbstractPricing)   # the pricing-door line renames the type;
    @eval const AbstractPricing = OpenCGE.AbstractPricingClosure   # the release exports it
end
isdefined(Main, :load_design) || include(joinpath(@__DIR__, "run.jl"))

function main()
    root = default_root()
    dd = load_design("matrix_5x3_v10"; root = root)
    ref = build_reference(dd; root = root)
    drops = Vector{Int}(dd["data"]["drops"])
    kept = sort(setdiff(1:(length(ref.data.factor_share) + length(drops)), drops))
    ψ, g = programme_vectors(dd, ref.data, kept; root = root)

    probe(sol) = begin
        b = final_demand_bundle(sol.model, sol.prices_raw, sol.quantities,
            OpenCGE._bundle_wage(sol); lambda = 1.0, external_transfer = sol.external_transfer)
        w = sol.wages_raw[1]
        L = sum(b.L_i)
        (w = w, L = L, wage_bill = w * L, E = b.E, R = sol.external_transfer,
         real = real_consumption(sol))
    end
    r0 = probe(ref.sol)
    g0 = 100 * (r0.real)

    for lab in ["GAMMA", "DELTA"]
        for fin in ["F2", "F3"]
            sol = solve_cell(design_cell(dd, "matrix_5x3-v10-$lab-$fin"), dd, ref.data, ψ, g,
                ref.init_warm)
            r = probe(sol)
            man = TOML.parsefile(joinpath(root, "runs", "matrix_5x3-v10-$lab-$fin", "manifest.toml"))
            @printf("\n== %s-%s ==  (manifest: Delta Y %+.2f, W_C %+.2f)\n", lab, fin,
                100 * man["metrics"]["gdp_rel"], 100 * man["metrics"]["consumption_rel"])
            @printf("   wage      %8.6f -> %8.6f   (%+.4f %%)\n", r0.w, r.w, 100 * (r.w / r0.w - 1))
            @printf("   employment%8.6f -> %8.6f   (%+.4f %%)\n", r0.L, r.L, 100 * (r.L / r0.L - 1))
            @printf("   wage bill %8.6f -> %8.6f   (%+.4f %%)\n", r0.wage_bill, r.wage_bill,
                100 * (r.wage_bill / r0.wage_bill - 1))
            @printf("   income E  %8.6f -> %8.6f   (%+.4f %%)\n", r0.E, r.E, 100 * (r.E / r0.E - 1))
            @printf("   R         %8.6f -> %8.6f\n", r0.R, r.R)
            @printf("   W_C                              (%+.4f %%)   <- the reported welfare change\n",
                100 * (r.real / r0.real - 1))
            @printf("   implied (1-tau): %+.4f %%   [income %% minus wage-bill %%]\n",
                (100 * (r.E / r0.E - 1)) - (100 * (r.wage_bill / r0.wage_bill - 1)))
        end
    end
    return 0
end

exit(main())
