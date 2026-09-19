# probe18_gamma_bf_eta.jl -- GAMMA carrying the BF allocation elasticity.
#
# READ-ONLY probe: no run, no src/ change, nothing committed.
#   julia --project=. experiments/probes/probe18_gamma_bf_eta.jl
#
# The BF closure's parameter IS the allocation elasticity eta (ADR-0010): the
# BF row of the matrix is eta = 0, and every other row runs eta = 1. The fixed
# wage row (GAMMA) runs eta = 1 -- the cost-minimizing allocation with the wage
# pinned -- so the *combination* "fixed wage + BF allocation rule" is a distinct
# variant that no cell of the matrix covers: the wage pinned AND the allocation
# frozen. It is the double-rigidity corner, and it is admissible in the kernel
# (eta in {0, 1} are the two admitted formulations; the interpolation 0 < eta < 1
# was retired with the B&F reallocation wedge, ADR-0010).
#
# Economics: with w = 1 and L_i = Lbar_i, nothing can absorb the programme --
# neither the wage, nor employment, nor (under a demand-only shock) the price
# level, since zero profit with w = 1 gives p = 1. The programme then moves only
# the *composition* of output and the tax/transfer mix. It is the corner the
# labour axis would sit at if labour were immobile and wages rigid at once, and
# the natural benchmark against which the other GAMMA variants are read.
#
# Reported: max|p-1|, deflator, employment (a datum here), consumption, the
# external position, and -- to show what does move -- the output-composition
# spread max_i |y_i/lambda_i - 1|.

using BeyondHulten
using CSV, DataFrames, LinearAlgebra, Printf

const ROOT = joinpath(@__DIR__, "..", "..") |> normpath

data = recalibrate_open(retained_dataset(read_data("I-O_DE2019_formatiert.csv"; datadir = ROOT),
    Vector{Int}([])); exo_scale = 1.0)
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
             id == "F2" ? TaxFinanced(g) : ExternalDebt(g)

println("N = ", N, "  programme sectors = ", count(>(0.0), ψ))
println("\nrow                       fin  eta  max|p-1|   deflator   L        cons      F        max|dy/lam-1|")

for (tag, η) in (("GAMMA (cost-min alloc)", 1.0), ("GAMMA + BF rule (eta=0)", 0.0))
    for fin_id in ("F1", "F2", "F3")
        m = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, η;
            closure = :fixed, financing = fin_of(fin_id))
        sol = try
            solve(m; init = [ones(N); data.λ])
        catch e
            println(@sprintf("%-24s %-4s %.1f  FAILED: %s", tag, fin_id, η,
                first(split(sprint(showerror, e), "\n"))))
            continue
        end
        p, y = sol.prices_raw, sol.quantities
        L = sum(sectoral_labor_demand(p, y, 1.0, sol.model))
        @printf("%-24s %-4s %.1f  %.6f   %.6f   %.6f  %+.6f  %+.5f  %.6f\n",
            tag, fin_id, η, maximum(abs, p .- 1), gdp_deflator(sol, ref_sol), L,
            real_consumption(sol) / ref_cons - 1, sol.external_transfer,
            maximum(abs, y ./ data.λ .- 1))
    end
end
println("\ndone.")