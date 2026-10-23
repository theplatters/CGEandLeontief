# probe21_sectoral_family_limit.jl -- the elastic end of the sectoral family.
#
# READ-ONLY probe: no run is created, no src/ file is touched, nothing committed.
#   julia --project=. experiments/probes/probe21_sectoral_family_limit.jl
#
# Why: Appendix A and Appendix B assert that the sectoral elastic family
# "approaches" the employment-adjustment closure as eta_s,i -> infinity. This
# probe asks whether that sequence has a limit at all, and where it lands.
#
# Answer (2026-10-23): it has a limit, and the limit is NOT the matrix row.
# Wages and prices converge to the baseline (max|p-1| falls 1e-1 -> 1e-4 over
# eta_s 0.5 -> 100), employment converges to 1.002327 and consumption to
# 0.987117, and the external transfer converges to F* = +0.00086. The
# wage-pinned row of the matrix (F = 0) reads 1.001260 and 0.984667. The gap
# is the transfer's contribution to demand: +0.107 pp of employment and
# +0.245 pp of consumption.
#
# Why the limit is the F-free member: in the mobile system the labour-market
# equation is what determines F (Appendix B, eq. finstrument). At infinite
# elasticity that equation degenerates into "the real wage is at its anchor"
# but remains an equation with a degree of freedom -- and it is the one that
# pins F. The wage-pinned row deletes that equation and pins F = 0 instead to
# keep the system square.
#
# Design: demand-only shocks, theta = eps = 0.5, sigma = 0.9, eta = 1, the
# 2024 programme (G0 = 1.331 % of GDP), uniform sectoral eta_s ladder, F free.

using OpenCGE
using CSV, DataFrames, LinearAlgebra, Printf

const ROOT = joinpath(@__DIR__, "..", "..") |> normpath

# ── calibration (mirrors experiments/run.jl) ────────────────────────────────
data_full = read_data("I-O_DE2019_formatiert.csv"; datadir = ROOT)
data_v1 = retained_dataset(data_full, Vector{Int}([]))
N = length(data_v1.factor_share)
data = recalibrate_open(data_v1; exo_scale = 1.0)

imp = CSV.read(joinpath(ROOT, "cbase2/data_raw/impulses.csv"), DataFrame)
rows = imp[imp.year .== 2024, :]
raw = Matrix{Float64}(rows[1:1, 3:73])[:]
v = raw[collect(1:N)]
ψ = v ./ sum(v)
g = (40300.0 / data.gdp_production) .* ψ

include(joinpath(ROOT, "experiments", "run.jl"))
fin1 = cell_financing("F1", ψ, g, "tilt_g0_over_c0", data)
fin2 = cell_financing("F2", ψ, g, "tilt_g0_over_c0", data)
fin3 = cell_financing("F3", ψ, g, "tilt_g0_over_c0", data)

shocks = Shocks(ones(N), ones(N), zeros(N))
ref_sol = solve(mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 1.0);
    init = [ones(N); data.λ; 1.0; 0.0])
ref_cons = real_consumption(ref_sol)
cpi_of(p, σ) = sum(vec(data.consumption_share) .* p .^ (1 - σ))^(1 / (1 - σ))

function report(tag, m, x0)
    sol = try
        solve(m; init = x0)
    catch e
        println(rpad(tag, 30), "SOLVE FAILED: ", first(split(sprint(showerror, e), "\n")))
        return nothing
    end
    p, q, w = sol.prices_raw, sol.quantities, sol.wages_raw
    L = sum(sectoral_labor_demand(p, q, w, m))
    @printf("%-30s L=%.6f cons=%+.6f max|p-1|=%.3e F=%+.6f w=[%.5f,%.5f]\n",
        tag, L, real_consumption(sol) / ref_cons, maximum(abs, p .- 1),
        sol.external_transfer, minimum(w) / cpi_of(p, 0.9), maximum(w) / cpi_of(p, 0.9))
    return sol
end

println("N=", N, "  G0=", @sprintf("%.5f", sum(g)))
println("\n--- wage-pinned reference (GAMMA, F = 0) ---")
for (nm, fin) in (("GAMMA-F1", fin1), ("GAMMA-F2", fin2), ("GAMMA-F3", fin3))
    report(nm, mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 1.0;
        closure = :fixed, financing = fin), [ones(N); data.λ])
end

println("\n--- sectoral elastic family, uniform eta_s, F free (F2) ---")
x0 = [ones(N); data.λ; ones(N); 0.0]
for η in (0.5, 2.0, 10.0, 100.0, 1000.0)
    s = report("eta_s=$η (F2)",
        mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 1.0;
            financing = fin2, eta_s_vec = fill(η, N)), x0)
    s === nothing || (global x0 = [s.prices_raw; s.quantities; s.wages_raw; s.external_transfer])
end

println("\n--- same ladder, F3 (neutrality along the sequence) ---")
x0 = [ones(N); data.λ; ones(N); 0.0]
for η in (2.0, 100.0)
    s = report("eta_s=$η (F3)",
        mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 1.0;
            financing = fin3, eta_s_vec = fill(η, N)), x0)
    s === nothing || (global x0 = [s.prices_raw; s.quantities; s.wages_raw; s.external_transfer])
end
println("\ndone.")
