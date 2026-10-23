# =====================================================================
# Table 6 --- W_D for the capacity-door cells (GENERATED sidecar)
# ---------------------------------------------------------------------
# Do NOT edit the output by hand. Regenerate with:
#   julia --project=. experiments/make_table6_wd.jl
# then rebuild the table with revision/tables/make_pricing_door.py.
#
# Uses the committed OpenCGE v0.5.0 release pin (ADR-0044).
#
# Job: W_D is not a manifest field, so it is computed here, exactly as in
# Table 5: the quantity index of demand, sum(c + g) / sum(c0 + g0) - 1, at
# lambda = 1, against the design's reference continuation. Every cell Table 6
# carries (the d0, dp05-s0 and dm05-s0 cells of the three capacity designs) is
# re-solved and gated against its runs/<id>/manifest.toml: every metric the
# manifest records must be reproduced within GATE_TOL, or the script fails and
# writes nothing. No run directory is created or modified.
#
# Output: revision/tables/table6_wd.csv (run_id, w_d in percent).
# =====================================================================

using OpenCGE, TOML, Printf

isdefined(Main, :load_design) || include(joinpath(@__DIR__, "run.jl"))

const LAMBDA = 1.0
const GATE_TOL = 1e-6
const DESIGNS = ["pricing_capacity_rows_v2", "pricing_capacity_v3", "pricing_mobile_v3"]
const TABLE6_CELL = r"-(d0|dp05-s0|dm05-s0)$"

"""W_D in percent: the own-weights quantity index of the c + lambda*g bundle."""
function w_d(sol::Solution, ref_sol::Solution)
    f(s) = final_demand_bundle(s.model, s.prices_raw, s.quantities, OpenCGE._bundle_wage(s);
        lambda = LAMBDA, external_transfer = s.external_transfer).f
    return 100 * (sum(f(sol)) / sum(f(ref_sol)) - 1)
end

"""Largest deviation between the re-solved metrics and the manifest's."""
function manifest_gap(got::Dict{String,Any}, man::Dict{String,Any})
    gap = 0.0
    for (k, v) in man["metrics"]
        (v isa Real && haskey(got, k) && got[k] isa Real) || continue
        gap = max(gap, abs(got[k] - v))
    end
    return gap
end

function main()
    root = default_root()
    rows = Tuple{String,Float64}[]
    ok = true
    for name in DESIGNS
        dd = load_design(name; root = root)
        ids = sort([id for id in keys(dd["cells"]) if occursin(TABLE6_CELL, id)])
        @printf("building reference continuation for %s (%d cells) ...\n", name, length(ids))
        ref = build_reference(dd; root = root)
        drops = Vector{Int}(dd["data"]["drops"])
        kept = sort(setdiff(1:(length(ref.data.factor_share) + length(drops)), drops))
        ψ, g = programme_vectors(dd, ref.data, kept; root = root)
        for id in ids
            man = TOML.parsefile(joinpath(root, "runs", id, "manifest.toml"))
            man["provenance"]["opencge"]["repo_rev"] == "v0.5.0" || error("wrong kernel generation: $id")
            cell = design_cell(dd, id)
            sol = solve_cell(cell, dd, ref.data, ψ, g, ref.init_warm)
            ev = evaluate_gates(cell, dd, sol.model, sol, ref.sol)
            gap = manifest_gap(ev.metrics, man)
            pass = man["status"] == "executed" && gap < GATE_TOL
            ok &= pass
            wd = w_d(sol, ref.sol)
            @printf("  %-42s gate %s (%.1e)  W_C %+7.3f  W_D %+7.3f\n", id, pass ? "ok" : "FAIL",
                gap, 100 * man["metrics"]["consumption_rel"], wd)
            push!(rows, (id, wd))
        end
    end
    ok || (println("manifest gate failed; nothing written"); return 1)
    path = joinpath(root, "revision", "tables", "table6_wd.csv")
    open(path, "w") do io
        println(io, "run_id,w_d")
        for (id, wd) in sort(rows)
            @printf(io, "%s,%.12g\n", id, wd)
        end
    end
    println("wrote ", path)
    return 0
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    exit(main())
end
