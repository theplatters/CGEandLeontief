# test_model.jl — Quick test of the equilibrium solver
#
# Run: julia --project=. src/test_model.jl

using Printf, LinearAlgebra, Statistics

push!(LOAD_PATH, @__DIR__)
include("io_table.jl")
include("shocks.jl")
include("network.jl")
include("model.jl")

println("="^60)
println("Testing MCP Equilibrium Solver")
println("="^60)

# ── 1. Load data ──
DATA_DIR = "data"
N, YEAR = 66, 2015

io = load_io_table(joinpath(DATA_DIR, "IO_data_2018.mat"); N=N, year=YEAR)
shocks = load_shocks(DATA_DIR; N=N)
sf = build_standard_form(io)

# ── 2. Build model parameters ──
D = sf.D

# Elasticities: benchmark regime (complementarity, Fig 2)
theta = make_theta(N; sigma=1.0, epsilon=0.6, eta=0.6, theta1=0.2)

# Cobb-Douglas indicator
cd_vec = make_cobb_douglas(N, D, sf.factor; benchmark=true)

# Shock vectors at t=0 (no shock)
A_init = ones(D)
B_init = ones(D, D)  # all ones — demand shock is applied through row 1 and 5N+3

# Initial prices and lambda (from Domar weights)
init_p = ones(D)
init_lambda = vec(sf.Psi_re[1, :])  # first row of Leontief inverse → column vector

# Ensure positivity for lambda that might be 0
for i in 1:D
    if init_lambda[i] ≤ 0
        init_lambda[i] = 1e-10
    end
end

# HtM phi_htm (all 1 = full insurance)
phi_htm = ones(D)
# Labor gets partial insurance: phi_htm[3N+2:4N+1] = 0.8
phi_htm[(3*N+2):(4*N+1)] .= 0.8

println("\n── Test 1: Equilibrium at t=0 (no shock) ──")
println("Expected: p = ones(D), λ = Psi_re[1,:] (Domar weights)")
println("Residual norm should be near 0.")

m0 = MCPModel(D, N, sf.Omega_re, sf.factor, sf.keynes,
              theta, cd_vec, sf.phi, phi_htm,
              ones(D), ones(D),
              A_init, B_init,
              init_lambda, init_p, sf.chi)

z0 = vcat(init_p, init_lambda)
F_test = zeros(2*D)
equilibrium_residual!(F_test, z0, m0)
resid0 = norm(F_test)
@printf("Residual norm at initial point: %.2e\n", resid0)
if resid0 < 1e-10
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
    p, λ, converged, iters = solve_equilibrium(m0, z0=z0, tol=1e-8, maxiter=200)
    @printf("Converged: %s\n", converged ? "yes" : "no")
    @printf("Solver iterations: %d\n", iters)
    @printf("Max price deviation from 1: %.2e\n", maximum(abs.(p .- 1.0)))
    @printf("p[%d] (numeraire): %.10f\n", D, p[D])
    if abs(p[D] - 1.0) < 1e-8
        println("✅ Numeraire holds (p[D] = 1).")
    end
catch e
    println("❌ Solver failed: $e")
end

println("\n── Test 3: Continuation with tiny shock (t=0.001) ──")
try
    t_grid = [0.0, 0.001, 0.01]
    result = eq_continuation(m0, shocks.BLS_shock, shocks.PCE_shock,
                              t_grid, show_trace=true)
    @printf("\nResults:\n")
    for ti in 1:length(t_grid)
        @printf("  t=%.4f  GDP=%.6f  Hulten=%.6f  nomGDP=%.6f  retcode=%d\n",
                t_grid[ti], result["GDP"][ti], result["Hulten"][ti],
                result["nominal_GDP"][ti], result["retcodes"][ti])
    end
    println("✅ Continuation completed.")
catch e
    println("❌ Continuation failed: $e")
    # Try just t=0
    println("\n── Fallback: Solve t=0 only ──")
    p, λ, converged, iters = solve_equilibrium(m0, z0=z0, tol=1e-8, maxiter=100)
    @printf("Converged: %s, iters=%d, p[D]=%.6f\n", converged ? "yes" : "no", iters, p[D])
end