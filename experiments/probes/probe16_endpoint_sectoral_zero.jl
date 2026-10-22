# probe16_endpoint_sectoral_zero.jl -- is no-reallocation the η_s = 0 endpoint
# of the sectoral labour-supply family?
#
# READ-ONLY probe: creates no run dir, touches no src/ file, commits nothing.
#   julia --project=. experiments/probes/probe16_endpoint_sectoral_zero.jl
#
# Question (revision: Section 3 vs Appendix A).
#   Section 3 says the elastic-supply family "spans the two rigid closures
#   without terminating in either"; Appendix A's coincidence table says
#   no-reallocation IS "the endpoint in which the sectoral family terminates,
#   exact by construction". Both cannot be right.
#   The executed cells cannot settle it: every BF cell runs η = 0 AND η_s = 0
#   jointly, and every sectoral cell runs η_s > 0. The missing cell is
#   η = 1 with η_s,i = 0 for all i.
#
# Method: use the v10 design's own calibration, reference continuation and
# warm start, so the numbers are directly comparable to matrix_5x3-v10-BF-*.
# Solve BF (η = 0) and the sectoral η_s = 0 cell per financing regime, then
# compare the solutions component by component.

using OpenCGE, Printf

isdefined(Main, :load_design) || include(joinpath(@__DIR__, "..", "run.jl"))

function main()
    root = default_root()
    dd = load_design("matrix_5x3_v10"; root = root)
    ref = build_reference(dd; root = root)
    drops = Vector{Int}(dd["data"]["drops"])
    kept = sort(setdiff(1:(length(ref.data.factor_share) + length(drops)), drops))
    ψ, g = programme_vectors(dd, ref.data, kept; root = root)
    N = length(ref.data.factor_share)
    println("N = ", N, "  G0 = ", @sprintf("%.6f", sum(g)))

    for fin in ["F1", "F2", "F3"]
        cbf = design_cell(dd, "matrix_5x3-v10-BF-$fin")
        sbf = solve_cell(cbf, dd, ref.data, ψ, g, ref.init_warm)
        ebf = evaluate_gates(cbf, dd, sbf.model, sbf, ref.sol)

        c0 = Dict{String,Any}("labor" => "BETA", "financing" => fin, "eta" => 1.0,
            "eta_s" => 0.0, "theta" => 0.5, "epsilon" => 0.5, "sigma" => 0.9,
            "eta_s_rigid_group" => "none",
            "note" => "endpoint test: sectoral supply at eta_s = 0, eta = 1")
        shocks, _ = resolve_supply_shock(dd, c0, ref.data, ψ)
        finobj = cell_financing(c0["financing"], ψ, g, dd["programme"]["f1_shift"], ref.data)
        m0 = mobile_labor_model(ref.data, shocks, 0.5, 0.5, 0.9, 1.0;
            financing = finobj, eta_s_vec = zeros(N))
        x0 = [sbf.prices_raw; sbf.quantities; sbf.wages_raw; sbf.external_transfer]
        s0 = try
            solve(m0; init = x0)
        catch e
            println("$fin  sectoral eta_s=0  FAILED: ",
                first(split(sprint(showerror, e), "\n")))
            continue
        end
        e0 = evaluate_gates(c0, dd, m0, s0, ref.sol)
        mb, mz = ebf.metrics, e0.metrics
        @printf("\n── %s ────────────────────────────────────────────────────────\n", fin)
        @printf("  BF  (η=0, η_s=0)      gdp %+8.5f  defl %+8.5f  L %+9.5f  W %+9.5f\n",
            100 * mb["gdp_rel"], 100 * (mb["gdp_deflator"] - 1),
            100 * (mb["employment"] - 1), 100 * mb["consumption_rel"])
        @printf("  S0  (η=1, η_s,i=0)    gdp %+8.5f  defl %+8.5f  L %+9.5f  W %+9.5f\n",
            100 * mz["gdp_rel"], 100 * (mz["gdp_deflator"] - 1),
            100 * (mz["employment"] - 1), 100 * mz["consumption_rel"])
        @printf("  |Δp| %.3e  |Δq| %.3e  |Δw| %.3e  |ΔF| %.3e  max|p-1| %.3e\n",
            maximum(abs, s0.prices_raw .- sbf.prices_raw),
            maximum(abs, s0.quantities .- sbf.quantities),
            maximum(abs, s0.wages_raw .- sbf.wages_raw),
            abs(s0.external_transfer - sbf.external_transfer),
            maximum(abs, s0.prices_raw .- 1))
        @printf("  gates S0: %s\n", e0.gate_summary)
    end
    println("\ndone.")
    return 0
end

exit(main())
