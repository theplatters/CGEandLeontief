# probe17_gamma_capacity.jl -- GAMMA (fixed wage) + the capacity/Verdoorn
# externality in the cost hook, the one GAMMA variant with no measurement yet.
#
# READ-ONLY probe: no run, no src/ change, nothing committed.
#   julia --project=. experiments/probes/probe17_gamma_capacity.jl
#
# The reduced form (workplan S5, probe14): A_eff_i = (y_i / lambda_i)^(-delta)
# enters the unit-cost hook. In the fixed-wage system the price block is no
# longer separable from quantities, so the model is a simultaneous system. It is
# solved here by the fixed-point iteration the reduced form implies:
#
#     A <- A .* (y/lambda)^(-delta),  re-solve the 2N fixed-wage system
#
# At the fixed point A is consistent with the solved y, so the returned solution
# IS the solution of the simultaneous system (not an approximation of it). The
# iteration is damped geometrically; convergence is reported, and non-convergence
# is itself a finding (delta < 0 is the Kaldor-Verdoorn sign: higher output
# lowers cost, a positive feedback).
#
# Canary: delta = 0 must reproduce the executed GAMMA cells (matrix_5x3-v9-
# GAMMA-*) exactly, since the externality then vanishes identically.

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

shocks0 = Shocks(ones(N), ones(N), zeros(N))
ref_sol = solve(mobile_labor_model(data, shocks0, 0.5, 0.5, 0.9, 1.0);
    init = [ones(N); data.λ; 1.0; 0.0])
ref_cons = real_consumption(ref_sol)

function f1_tilt(baseline, ψ, g)
    pos = baseline .> 0
    ψ1 = ψ .* pos; ψ1 = ψ1 ./ sum(ψ1)
    return 1.0 .+ sum(g) .* ψ1 ./ max.(baseline, 1e-12)
end
fin_of(id) = id == "F1" ? PreferenceReallocation(f1_tilt(data.household_baseline, ψ, g)) :
             id == "F2" ? TaxFinanced(g) : ExternalDebt(g)

"""
Solve the fixed-wage system with the externality by fixed-point iteration.
Returns `(sol, iterations, last update, max|resid|)` or `(nothing, it, Δ, NaN)`
when it fails to converge.
"""
function gamma_capacity(fin, δ::Float64; maxit = 80, tol = 1e-12, damp = 0.6)
    A = ones(N)
    init = nothing
    lastΔ = NaN
    for it in 1:maxit
        m = mobile_labor_model(data, Shocks(A, ones(N), zeros(N)), 0.5, 0.5, 0.9, 1.0;
            closure = :fixed, financing = fin)
        sol = try
            solve(m; init = init)
        catch e
            return nothing, it, NaN, NaN
        end
        y = sol.quantities
        target = (y ./ data.λ) .^ (-δ)
        Anew = exp.((1 - damp) .* log.(A) .+ damp .* log.(target))
        lastΔ = maximum(abs, Anew .- A)
        init = [sol.prices_raw; y]
        A = Anew
        lastΔ < tol && return sol, it, lastΔ, maximum(abs, equilibrium_residuals(m, init))
    end
    return nothing, maxit, lastΔ, NaN
end

println("N = ", N, "  programme sectors = ", count(>(0.0), ψ),
    "  G0 = ", @sprintf("%.5f", sum(g)))
println("\ncanary: delta = 0 must reproduce the executed GAMMA cells")
println("row     fin   delta  it   max|p-1|   deflator   L        consumption   converged")

for fin_id in ("F1", "F2", "F3")
    for δ in (0.0, 0.5, -0.5)
        fin = fin_of(fin_id)
        sol, it, Δ, resid = gamma_capacity(fin, δ)
        tag = δ == 0.0 ? "GAMMA" : "GAMMA+cap"
        if sol === nothing
            @printf("%-8s %-5s %+5.1f  %3d  --- no convergence (last |dA| = %.2e)\n",
                tag, fin_id, δ, it, Δ)
        else
            p, y = sol.prices_raw, sol.quantities
            L = sum(sectoral_labor_demand(p, y, 1.0, sol.model))
            @printf("%-8s %-5s %+5.1f  %3d  %.6f   %.6f   %.6f  %+.6f    yes (|dA|=%.1e)\n",
                tag, fin_id, δ, it, maximum(abs, p .- 1), gdp_deflator(sol, ref_sol),
                L, real_consumption(sol) / ref_cons - 1, Δ)
        end
    end
end

println("\nexecuted GAMMA cells for comparison (matrix_5x3-v9-GAMMA-*):")
for fin_id in ("F1", "F2", "F3")
    m = mobile_labor_model(data, shocks0, 0.5, 0.5, 0.9, 1.0; closure = :fixed,
        financing = fin_of(fin_id))
    sol = solve(m; init = [ones(N); data.λ])
    p, y = sol.prices_raw, sol.quantities
    L = sum(sectoral_labor_demand(p, y, 1.0, sol.model))
    @printf("GAMMA    %-5s  0.0    -   %.6f   %.6f   %.6f  %+.6f\n",
        fin_id, maximum(abs, p .- 1), gdp_deflator(sol, ref_sol), L,
        real_consumption(sol) / ref_cons - 1)
end
println("\ndone.")