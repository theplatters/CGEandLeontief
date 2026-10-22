# =====================================================================
# Table 7 --- W_D for the markup-door cells (GENERATED sidecar)
# ---------------------------------------------------------------------
# Do NOT edit the output by hand. Regenerate with:
#   julia --project=<OPEN_CGE_CHECKOUT> experiments/make_table7_wd.jl
# then rebuild the table with revision/tables/make_pricing_door.py.
#
#   NOT `--project=.`: like make_table6_wd.jl, this script calls
#   `final_demand_bundle` (OpenCGE diag(WELFARE), branch welfare-indices =
#   v0.2.1 + the welfare diagnostics), which the paper repo's Manifest does
#   not yet pin. Point --project at a checkout that contains it.
#
# Job: Table 7 reports the same welfare measure as Table 6 (W_D, the own-weights
# quantity index of the c + lambda*g bundle), which no manifest records. Every
# cell Table 7 carries is re-solved here and gated against its
# runs/<id>/manifest.toml: every metric the manifest records must be reproduced
# within GATE_TOL, or the script fails and writes nothing. No run directory is
# created or modified.
#
# Cells (75): the kappa = 0 reference row (the matrix_5x3-v10 no-door twins),
# the retained markup cells (kappa in {0.25, 0.5}) and the kappa = 0.5
# distributed cells (s_pi in {0.5, s}).
#
# Output: revision/tables/table7_wd.csv (run_id, w_d in percent).
# =====================================================================

using OpenCGE, TOML, Printf

isdefined(Main, :load_design) || include(joinpath(@__DIR__, "run.jl"))

const LAMBDA = 1.0
const GATE_TOL = 1e-6
const TABLE7_CELLS = Dict(
    # the no-door row of Table 7 uses the eta_s = 0.5 cell on BETA (the plain
    # matrix_5x3-v10-BETA-F1 run is stale and reproduces ALPHA), as Table 4/6 do.
    "matrix_5x3_v10"     => r"^matrix_5x3-v10-([A-Z]+-F[123]|BETA-F[123]-etas05)$",
    "pricing_markup_v2"  => r"^pricing_markup-v2-[A-Z]+-F[123]-k(025|05)$",
    "pricing_profit_v2"  => r"^pricing_profit-v2-[A-Z]+-F[123]-k05-spi(05|s)$",
)

# The stale plain BETA cells of matrix_5x3_v10 (eta recorded as 0.5 but reproducing
# ALPHA) are replaced by the eta_s = 0.5 cells; never solve them.
const SKIP = r"^matrix_5x3-v10-BETA-F[123]$"

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
    for (name, re) in TABLE7_CELLS
        dd = load_design(name; root = root)
        ids = sort([id for id in keys(dd["cells"]) if occursin(re, id) && !occursin(SKIP, id)])
        @printf("building reference continuation for %s (%d cells) ...\n", name, length(ids))
        ref = build_reference(dd; root = root)
        drops = Vector{Int}(dd["data"]["drops"])
        kept = sort(setdiff(1:(length(ref.data.factor_share) + length(drops)), drops))
        ψ, g = programme_vectors(dd, ref.data, kept; root = root)
        for id in ids
            man = TOML.parsefile(joinpath(root, "runs", id, "manifest.toml"))
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
    path = joinpath(root, "revision", "tables", "table7_wd.csv")
    open(path, "w") do io
        println(io, "run_id,w_d")
        for (id, wd) in sort(rows)
            @printf(io, "%s,%.12g\n", id, wd)
        end
    end
    println("wrote ", path, "  (", length(rows), " cells)")
    return 0
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    exit(main())
end
