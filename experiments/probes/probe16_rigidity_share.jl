# probe16_rigidity_share.jl -- how much does the GROUPING RULE matter?
#
# READ-ONLY probe: no run is created, no src/ file is touched, nothing is
# committed. Run from the repository root:
#   julia --project=. experiments/probes/probe16_rigidity_share.jl
#
# Why: the executed matrix_5x3_v9 cells show that the rigid-GROUP choice moves
# the price response by a factor 2.5 and flips the sign of the employment effect
# (rigid programme sectors: -0.34 %; rigid largest half: +0.02 %; uniform:
# +0.13 %). That is a two-point comparison. This probe maps the whole surface:
# for a FIXED ranking, sweep the rigid SHARE s in {0, 0.25, 0.5, 0.75, 1}; and
# for a fixed share (0.5), sweep three RANKINGS (programme incidence descending,
# baseline employment descending, programme incidence ASCENDING). Together they
# separate the two things a grouping rule does: how MUCH of the economy is rigid,
# and WHICH sectors are.
#
# The rigid group gets eta_s,i = 0 (a vertical supply curve); the flexible group
# gets eta_s = 0.5. s = 0 is the uniform cell, s = 1 is the eta = 0 endpoint.
#
# Reported per cell: max|p-1| at F1 and F2, the F1-F2 difference (the
# demand-sensitivity signature), employment, consumption, the deflator, and the
# wage spread.

using BeyondHulten
using CSV, DataFrames, LinearAlgebra, Printf

const ROOT = joinpath(@__DIR__, "..", "..") |> normpath

data_full = read_data("I-O_DE2019_formatiert.csv"; datadir = ROOT)
data = recalibrate_open(retained_dataset(data_full, Vector{Int}([])); exo_scale = 1.0)
N = length(data.factor_share)

imp = CSV.read(joinpath(ROOT, "cbase2/data_raw/impulses.csv"), DataFrame)
rows = imp[imp.year .== 2024, :]
ψ = Matrix{Float64}(rows[1:1, 3:73])[:][collect(1:N)]
ψ = ψ ./ sum(ψ)
g = (40300.0 / data.gdp_production) .* ψ

shocks = Shocks(ones(N), ones(N), zeros(N))
ref_sol = solve(mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 1.0);
    init = [ones(N); data.λ; 1.0; 0.0])
ref_cons = real_consumption(ref_sol)

function f1_tilt(baseline, ψ, g)
    pos = baseline .> 0
    ψ1 = ψ .* pos; ψ1 = ψ1 ./ sum(ψ1)
    return 1.0 .+ sum(g) .* ψ1 ./ max.(baseline, 1e-12)
end
fin_of(id) = id == "F1" ? PreferenceReallocation(f1_tilt(data.household_baseline, ψ, g)) :
             id == "F2" ? TaxFinanced(g) : error("unknown financing")

rank_prog  = sortperm(ψ; rev = true)
rank_emp   = sortperm(data.labor_share; rev = true)
rank_least = sortperm(ψ)                       # rigid = the LEAST exposed sectors

function rigid_vec(s::Float64, rank::Vector{Int}; η = 0.5)
    v = fill(η, N)
    k = round(Int, s * N)
    k > 0 && (v[rank[1:k]] .= 0.0)
    return v
end

"""
Solve one cell at (F1, F2) and print the two-financing comparison.
"""
function row(tag, esv)
    out = Dict{String,Any}()
    for F in ("F1", "F2")
        m = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 1.0;
            financing = fin_of(F), eta_s_vec = esv)
        sol = try
            solve(m; init = [ones(N); data.λ; ones(N); 0.0])
        catch e
            println(@sprintf("%-34s %s FAILED: %s", tag, F, first(split(sprint(showerror, e), "\n"))))
            return
        end
        p, q, w = sol.prices_raw, sol.quantities, sol.wages_raw
        L = sum(sectoral_labor_demand(p, q, w, m))
        out[F] = (maxp = maximum(abs, p .- 1), L = L,
                  cons = real_consumption(sol) / ref_cons,
                  defl = gdp_deflator(sol, ref_sol),
                  spread = maximum(w) / minimum(w))
    end
    @printf("%-34s max|p-1| F1=%.6f F2=%.6f (d=%+.6f)  L=%.6f  cons=%+.5f  defl=%.6f  spread=%.3f\n",
        tag, out["F1"].maxp, out["F2"].maxp, out["F1"].maxp - out["F2"].maxp,
        out["F2"].L, out["F2"].cons, out["F2"].defl, out["F2"].spread)
end

println("N = ", N, "  programme sectors with psi > 0: ", count(>(0.0), ψ),
    "  G0 = ", @sprintf("%.5f", sum(g)))
println("\n── Rigid SHARE sweep, ranking = programme incidence (descending) ───────")
for s in (0.0, 0.25, 0.5, 0.75, 1.0)
    row("rigid share s=$(s) (programme rank)", rigid_vec(s, rank_prog))
end
println("\n── Rigid SHARE sweep, ranking = baseline employment (descending) ───────")
for s in (0.0, 0.25, 0.5, 0.75, 1.0)
    row("rigid share s=$(s) (employment rank)", rigid_vec(s, rank_emp))
end
println("\n── RANKING sweep at a fixed share (s = 0.5) ────────────────────────────")
row("rigid = programme sectors", rigid_vec(0.5, rank_prog))
row("rigid = largest employers", rigid_vec(0.5, rank_emp))
row("rigid = LEAST exposed sectors", rigid_vec(0.5, rank_least))
println("\n── Reference: uniform eta_s = 0.5 and the eta = 0 endpoint ─────────────")
row("uniform eta_s = 0.5", fill(0.5, N))
row("all rigid (eta_s = 0)", zeros(N))
println("\ndone.")