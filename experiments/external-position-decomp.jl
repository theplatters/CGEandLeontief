# experiments/external-position-decomp.jl
#
# Why is the booked external position Delta D = F + B_gov NEGATIVE in the
# flexible-wage closures?
#
# Delta D is not a policy choice: it is the position the real allocation
# implies, booked as the endogenous transfer F plus, under F3 only, the
# programme's own financing B_gov = sum(p.*g). This prints the seven
# gdp_components changes (in GDP units, baseline GDP = 1) alongside F, so the
# driver is visible rather than asserted.
#
#   julia --project=. experiments/external-position-decomp.jl

using OpenCGE, Printf

const AbstractPricing = OpenCGE.AbstractPricingClosure
isdefined(Main, :load_design) || include(joinpath(@__DIR__, "run.jl"))

const NAMES = ["C_gross", "G+prog", "I", "X", "-M_final", "-M_int", "-T_int"]

function main()
    root = default_root()
    dd = load_design("matrix_5x3_v10"; root = root)
    ref = build_reference(dd; root = root)
    drops = Vector{Int}(dd["data"]["drops"])
    kept = sort(setdiff(1:(length(ref.data.factor_share) + length(drops)), drops))
    ψ, g = programme_vectors(dd, ref.data, kept; root = root)
    r0 = gdp_components(ref.sol.model, ref.sol)

    for lab in ["BF", "ALPHA", "GAMMA"]
        println("\n== ", lab, " ==")
        @printf("%-9s %8s %8s | %s\n", "fin", "F", "B_gov", join(NAMES, " "))
        @printf("%-9s %8s %8s | %s\n", "", "", "", join(fill("--------", 7), " "))
        for fin in ["F1", "F2", "F3"]
            id = "matrix_5x3-v10-$lab-$fin"
            sol = solve_cell(design_cell(dd, id), dd, ref.data, ψ, g, ref.init_warm)
            r = gdp_components(sol.model, sol)
            d = 100 .* (r.V .- r0.V)
            @printf("%-9s %+8.3f %8.3f | %s\n", fin, 100 * r.external_transfer,
                100 * r.programme_financing, join([@sprintf("%+7.3f", x) for x in d], " "))
            @printf("%-9s booked = %+.2f  (F + B_gov);  sum(components) = %+.3f, wage bill = %+.3f\n",
                "", 100 * r.external_financing, sum(d), 100 * (r.wage_bill - r0.wage_bill))
        end
    end
    return 0
end

exit(main())
