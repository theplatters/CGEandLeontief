# test_model.jl — Test of the equilibrium solver (MATLAB-faithful calibration)
#
# Run: julia --project=. src/test_model.jl
#
# Uses make_model() which replicates Master_file_3.m lines 328–335:
# init_lambda[1]=1, HtM=0, Ricardian=2, tomorrow=1; phi_htm[labor]=1 (benchmark).

using Printf, LinearAlgebra, Statistics

push!(LOAD_PATH, @__DIR__)
include("io_table.jl")
include("shocks.jl")
include("network.jl")
include("model.jl")

println("="^60)
println("Testing MCP Equilibrium Solver (MATLAB-faithful calibration)")
println("="^60)

# ── 1. Load data ──
DATA_DIR = "data"
N, YEAR = 66, 2015

io = load_io_table(joinpath(DATA_DIR, "IO_data_2018.mat"); N=N, year=YEAR)
shocks = load_shocks(DATA_DIR; N=N)
sf = build_standard_form(io)
D = sf.D

# ── 2. Build model parameters (benchmark regime: complementarity, Fig 2) ──
m0 = make_model(sf; htm_share=0.0, benchmark=true,
                sigma=1.0, epsilon=0.6, eta=0.6, theta1=0.2)

@printf("D=%d, init_lambda[1]=%.6f, [D-2]=%.6f, [D-1]=%.6f, [D]=%.6f\n",
        D, m0.init_lambda[1], m0.init_lambda[D-2], m0.init_lambda[D-1], m0.init_lambda[D])

println("\n── Test 1: Residual at MATLAB-corrected initial point ──")
println("Expected: p = ones(D), λ = init_lambda (MATLAB-adjusted) → residual ≈ 0")
z0 = vcat(m0.init_p, m0.init_lambda)
F_test = zeros(2*D)
equilibrium_residual!(F_test, z0, m0)
resid0 = norm(F_test)
@printf("Residual norm at initial point: %.2e\n", resid0)
if resid0 < 1e-9
    println("✅ Initial point is an equilibrium (residual ≈ 0).")
else
    println("⚠️  Residual = $(resid0) — investigating largest components:")
    idx = sortperm(abs.(F_test), rev=true)[1:5]
    for i in idx
        if i ≤ D
            @printf("  F[%3d] (p)     = %.4e\n", i, F_test[i])
        else
            @printf("  F[%3d] (λ)     = %.4e\n", i-D, F_test[i])
        end
    end
end

println("\n── Test 2: Solve equilibrium at t=0 ──")
try
    p, λ, converged, iters = solve_equilibrium(m0, z0=z0, tol=1e-10, maxiter=200)
    @printf("Converged: %s\n", converged ? "yes" : "no")
    @printf("Solver iterations: %d\n", iters)
    @printf("Max price deviation from 1: %.2e\n", maximum(abs.(p .- 1.0)))
    @printf("p[%d] (numeraire): %.10f\n", D, p[D])
    @printf("λ[1] (nominal GDP): %.6f\n", λ[1])
    if abs(p[D] - 1.0) < 1e-8 && maximum(abs.(p .- 1.0)) < 1e-6
        println("✅ Numeraire holds and all prices ≈ 1 (trivial equilibrium).")
    else
        println("⚠️  Prices deviate from 1 — check init_lambda.")
    end
catch e
    println("❌ Solver failed: $e")
end

println("\n── Test 3: Continuation with small shock (t = 0, 0.01, 0.05) ──")
println("(Uses eq_continuation with shock_type=1 logic: supply + sectoral demand + aggregate demand)")
try
    t_grid = [0.0, 0.01, 0.05]
    result = eq_continuation(m0, shocks.BLS_shock, shocks.PCE_shock,
                              t_grid, show_trace=true)
    @printf("\nResults:\n")
    for ti in 1:length(t_grid)
        @printf("  t=%.4f  RGDP=%.6f  Hulten=%.6f  nomGDP=%.6f  retcode=%d\n",
                t_grid[ti], result["GDP"][ti], result["Hulten"][ti],
                result["nominal_GDP"][ti], result["retcodes"][ti])
    end
    if all(result["retcodes"] .== 0)
        println("✅ Continuation completed — all solves converged.")
    else
        println("⚠️  Some solves did not converge (retcode=1).")
    end
catch e
    println("❌ Continuation failed: $e")
    println("\n── Fallback: Solve t=0 only ──")
    p, λ, converged, iters = solve_equilibrium(m0, z0=z0, tol=1e-10, maxiter=100)
    @printf("Converged: %s, iters=%d, p[D]=%.6f\n", converged ? "yes" : "no", iters, p[D])
end
