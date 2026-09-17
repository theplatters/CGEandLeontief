# solver_sweep.jl -- init/algorithm sweep for the v3 mobile equilibrium.
# Purpose: establish whether the residual floors are init/path artifacts
# (the same system converged to 2e-10 under one init, stalled at 3.6e-4
# under another) and whether AutoFiniteDiff / LM-primary fix them.
# Run from the repo root:  julia --threads=4 cbase2/scripts/solver_sweep.jl
using CSV, DataFrames, LinearAlgebra, NonlinearSolve

const CB = joinpath(@__DIR__, "..")
for f in ["interface.jl", "solution.jl", "ces.jl", "mobile_labor.jl", "leontief.jl", "util.jl"]
    include(joinpath(CB, "src", "core", f))
end
include(joinpath(CB, "src", "closures.jl"))
include(joinpath(CB, "src", "financing.jl"))
include(joinpath(CB, "src", "calibration.jl"))

function main()
    drops = DATASET_VARIANTS["70s"]
    data_v1 = drop_sectors(read_data(joinpath(CB, "data_raw", "I-O_DE2019_formatiert.csv")), drops)
    N = length(data_v1.factor_share)
    shocks = Shocks(ones(N), ones(N), zeros(N))
    es = recalibrate_open(data_v1, CB; exo_scale=1.0, drops=drops)
    (; Ω_raw, factor_share, consumption_share, import_margin, gov_demand,
     exo_demand, exports_demand, saving_rate) = es
    M0 = Ω_raw' * Diagonal(1.0 .- factor_share)
    Gk = M0 + Diagonal(1.0 .- import_margin) * consumption_share * factor_share' *
         (1 - saving_rate) * (1 - sum(gov_demand))
    ylin = (I - Gk) \ ((1.0 .- import_margin) .* (gov_demand .+ exo_demand) .+ exports_demand)

    for θ in (2.0, 1.0, 0.5)
        mdl = mobile_labor_model(es, shocks, θ, 0.5, 0.9, 0.5)
        inits = (
            default = [ones(N); es.λ; 1.0],
            linear  = [ones(N); ylin; 1.0],
            perturb = [ones(N); es.λ .+ 0.02 .* sin.(1:N); 1.0],
        )
        for (nm, init) in pairs(inits)
            r0 = maximum(abs, equilibrium_residuals(mdl, init))
            prob = NonlinearSolve.NonlinearProblem(problem, init, mdl)
            res = NonlinearSolve.solve(prob, reltol=1e-8, abstol=1e-8, maxiters=20000)
            r = equilibrium_residuals(mdl, res.u)
            println("θ=", θ, " init=", rpad(nm, 8), " Newton: ret=", res.retcode,
                    " resid=", round(maximum(abs, r); digits=10))
            if maximum(abs, r) > 1e-6
                res2 = NonlinearSolve.solve(
                    NonlinearSolve.NonlinearProblem(problem, res.u, mdl),
                    LevenbergMarquardt(); reltol=1e-10, abstol=1e-10, maxiters=20000)
                r2 = equilibrium_residuals(mdl, res2.u)
                println("      init=", rpad(nm, 8), " LM polish: ret=", res2.retcode,
                        " resid=", round(maximum(abs, r2); digits=10))
            end
        end
    end
end
main()
