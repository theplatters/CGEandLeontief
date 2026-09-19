# probe15_sectoral_supply.jl -- the sectoral closure under a SUPPLY shock.
#
# READ-ONLY probe: no run is created, no src/ file is touched, nothing is
# committed. Run from the repository root:
#   julia --project=. experiments/probes/probe15_sectoral_supply.jl
#
# Why: the executed sectoral evidence (matrix_5x3_v9) is demand-only, and under
# demand-only shocks eta_s is NOT identified -- the ladder is a sensitivity band.
# The supply arm is the route to identification (Remark in paper/equivalence.tex
# v3, Section 4). This probe asks three questions.
#
#   Q1 (identification)  Does a supply shock separate the eta_s rungs? Under
#                        A[1] = 1.2 the real wage leaves its anchor, so the
#                        sectoral supply curves should bite.
#   Q2 (interaction)     Does the price response change in size or sign when the
#                        shock is a technology shock rather than a demand
#                        composition? Compare max|p-1| and the deflator against
#                        the executed demand-only cells.
#   Q3 (the canary)      Does the rigid corner still reproduce the eta = 0
#                        endpoint when A != 1 (i.e. is the nesting independent of
#                        the shock)? It should: the proof never uses A = 1.
#
# Design: A = ones(N) with A[1] = 1.2 (+20 % in sector 1, as in
# supply_identification.jl), theta = eps = 0.5, sigma = 0.9, eta = 1.
#   Row set 1  supply shock only (no programme), eta_s in {0, 0.5, 2}
#   Row set 2  supply shock + the 2024 programme under F2, eta_s in {0, 0.5, 2}
#   Row set 3  the canary: sectoral at eta_s,i = 0 vs the eta = 0 endpoint, with
#              and without the shock, F2.

using BeyondHulten
using CSV, DataFrames, LinearAlgebra, Printf

const ROOT = joinpath(@__DIR__, "..", "..") |> normpath

# ── calibration (mirrors experiments/run.jl) ────────────────────────────────
data_full = read_data("I-O_DE2019_formatiert.csv"; datadir = ROOT)
drops = Vector{Int}([])
data_v1 = retained_dataset(data_full, drops)
N = length(data_v1.factor_share)
data = recalibrate_open(data_v1; exo_scale = 1.0)

# ── programme incidence (2024 impulses, design matrix_5x3_v9) ──────────────
imp = CSV.read(joinpath(ROOT, "cbase2/data_raw/impulses.csv"), DataFrame)
rows = imp[imp.year .== 2024, :]
raw = Matrix{Float64}(rows[1:1, 3:73])[:]
v = raw[collect(1:N)]
ψ = v ./ sum(v)
g = (40300.0 / data.gdp_production) .* ψ

shocks_demand = Shocks(ones(N), ones(N), zeros(N))
A = ones(N); A[1] = 1.2
shocks_supply = Shocks(A, ones(N), zeros(N))

ref_sol = solve(mobile_labor_model(data, shocks_demand, 0.5, 0.5, 0.9, 1.0);
    init = [ones(N); data.λ; 1.0; 0.0])
ref_cons = real_consumption(ref_sol)

cpi_of(p, σ) = sum(vec(data.consumption_share) .* p .^ (1 - σ))^(1 / (1 - σ))

println("N = ", N, "  programme G0 = ", @sprintf("%.5f", sum(g)),
    "  A[1] = ", A[1])

"""
Solve one sectoral cell and print its metrics. `init` is the 3N+1 warm start.
"""
function cell_row(tag, shocks, esv, fin; init = nothing)
    m = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 1.0;
        financing = fin, eta_s_vec = esv)
    x0 = init === nothing ? [ones(N); data.λ; ones(N); 0.0] : init
    sol = try
        solve(m; init = x0)
    catch e
        println(@sprintf("%-26s FAILED: %s", tag, first(split(sprint(showerror, e), "\n"))))
        return nothing
    end
    p, q, w = sol.prices_raw, sol.quantities, sol.wages_raw
    X = [p; q; w; sol.external_transfer]
    resid = maximum(abs, equilibrium_residuals(m, X))
    Lsum = sum(sectoral_labor_demand(p, q, w, m))
    cpi = cpi_of(p, 0.9)
    cons = real_consumption(sol) / ref_cons
    @printf("%-26s resid=%.1e L=%.6f w/CPI=[%.5f,%.5f] Lmax/Lmin=%.4f max|p-1|=%.6f defl=%.6f cons=%+.5f F=%+.5f\n",
        tag, resid, Lsum, minimum(w) / cpi, maximum(w) / cpi,
        maximum(sectoral_labor_demand(p, q, w, m)) / minimum(sectoral_labor_demand(p, q, w, m)),
        maximum(abs, p .- 1), gdp_deflator(sol, ref_sol), cons, sol.external_transfer)
    return sol
end

println("\n── Row set 1: SUPPLY shock only (no programme) ─────────────────────────")
for η in (0.0, 0.5, 2.0)
    cell_row("A=1.2 eta_s=$η", shocks_supply, fill(η, N), NoFinancing())
end

println("\n── Row set 2: supply shock + programme (F2) ────────────────────────────")
for η in (0.0, 0.5, 2.0)
    cell_row("A=1.2 + g, eta_s=$η", shocks_supply, fill(η, N), TaxFinanced(g))
end

println("\n── Row set 3: the canary, with and without the shock ──────────────────")
for (nm, shocks) in (("A=1", shocks_demand), ("A[1]=1.2", shocks_supply))
    m0 = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 0.0; financing = TaxFinanced(g))
    s0 = solve(m0; init = [ones(N); data.λ; ones(N); 0.0])
    x0 = [s0.prices_raw; s0.quantities; s0.wages_raw; s0.external_transfer]
    sz = cell_row("  sectoral eta_s=0 ($nm)", shocks, zeros(N), TaxFinanced(g); init = x0)
    if sz !== nothing
        @printf("  canary %-10s max|dp|=%.2e max|dq|=%.2e max|dw|=%.2e |dF|=%.2e\n",
            nm, maximum(abs, sz.prices_raw .- s0.prices_raw),
            maximum(abs, sz.quantities .- s0.quantities),
            maximum(abs, sz.wages_raw .- s0.wages_raw),
            abs(sz.external_transfer - s0.external_transfer))
    end
end

println("\n── Reference points ────────────────────────────────────────────────────")
cell_row("demand-only, eta_s=0.5", shocks_demand, fill(0.5, N), TaxFinanced(g))
cell_row("demand-only, eta_s=0", shocks_demand, zeros(N), TaxFinanced(g))
println("\ndone.")