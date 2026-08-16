# oil_shock.jl — Historical oil-shock calibration (Oil_Shock.m)
#
# Parameters from Oil_Shock.m:
#   Year = 1971    (not 1980 — different IO structure)
#   ε = 0.02, θ = 0.25, σ = 0.25
#   Shock: A[7] = 0.9  (10% negative TFP shock to oil)
#   GDP: real_gdp_mc formula (same as MATLAB line)
#   Fixed-labor solver (Simulation_Derivs.m → solve_bf)

function run_oil_shock(data;
        A_shock_val::Float64 = 0.9,
        ε::Float64 = 0.02,
        θ::Float64 = 0.25,
        σ::Float64 = 0.25,
        sector::Int = 7)

    A = ones(data.N)
    A[sector] = A_shock_val

    sol = solve_bf(A, data.Ω, data.α, data.β, data.L, ε, θ, σ)
    if !sol.converged
        @warn "Oil-shock solver did not converge (residual: $(sol.residual_norm))"
    end

    # GDP = real_gdp_mc formula (matching MATLAB line)
    gdp = real_gdp_mc(sol, A, data.α, data.L, ε)

    # New Domar weight at shocked equilibrium
    λ_new = sol.p .* sol.y ./ gdp

    # Hulten benchmark
    λ = (I - Diagonal(1 .- data.α) * data.Ω)' \ data.β
    hulten = exp(dot(λ, log.(A)))

    return (gdp = gdp,
            log_gdp_change = log(gdp),
            hulten = hulten,
            amplification = log(gdp) / log(hulten),
            λ_new = λ_new,
            p = sol.p,
            y = sol.y,
            converged = sol.converged)
end