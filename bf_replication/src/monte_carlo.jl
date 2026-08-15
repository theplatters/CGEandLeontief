
# monte_carlo.jl -- Monte-Carlo robustness check (R3), ports GDP_Simulation_88sectorKLEMS.m.
# Reproduces the B&F (2019) Monte Carlo (Table I / Fig. 6):
#   * draw independent sectoral TFP shocks A = exp(z),  z ~ N(-0.5 * diag(Cov), Cov)
#     with Cov DIAGONAL (independent sectors) -- this is exactly the MATLAB line
#         A = exp(mvnrnd(-1/2*diag(Sigma), diag(diag(Sigma))))'
#     i.e. Cov = Diagonal(diag(Sigma_yearly)) for the paper's annual benchmark, or
#     Diagonal(diag(Sigma_4year)) for the JK 4-year variant.
#   * solve the equilibrium with solve_bf (numerical Newton, verified against the
#     pedagogical solver to 1e-13 on the oil-shock calibration).
#   * record the REAL (CES-welfare) GDP exactly as the MATLAB does (line 123):
#         C = sum( L .* p .* A.^((eps-1)/eps) .* alpha.^(1/eps) .* y.^(1/eps) .* (1 ./ L).^(1/eps) )
#     NOTE: this is NOT nominal GDP (w'L); the paper reports moments of log(C).
#   * keep draws with log(C) in (-0.4, 0.3) (the MATLAB "correct" filter) and
#     report mean / std / skewness / excess kurtosis of log(C).
using Random

"""
    real_gdp_mc(sol, A, α, L, ε)

Real (CES-welfare) GDP -- exact port of the MATLAB MC line 123:
    C = (p .* A.^((ε-1)/ε) .* α.^(1/ε) .* y.^(1/ε) .* (1 ./ L).^(1/ε))' * L
"""
function real_gdp_mc(sol::BFSolution, A::Vector{Float64}, α::Vector{Float64},
                     L::Vector{Float64}, ε::Float64)
    core = (A .^ ((ε - 1) / ε)) .* (α .^ (1 / ε)) .* (sol.y .^ (1 / ε)) .* ((1 ./ L) .^ (1 / ε))
    return sum(L .* sol.p .* core)
end

# Sample moments of log GDP, matching the MATLAB biased estimators:
#   skewness = m3 / m2^(3/2),  kurtosis = m4 / m2^2  ->  excess = kurtosis - 3.
function moments_loggdp(x::Vector{Float64})
    n = length(x)
    n == 0 && return (mean = NaN, std = NaN, skewness = NaN, excess_kurtosis = NaN)
    μ  = mean(x)
    d  = x .- μ
    m2 = sum(d .^ 2) / n
    m3 = sum(d .^ 3) / n
    m4 = sum(d .^ 4) / n
    return (mean = μ, std = sqrt(m2),
            skewness = m3 / m2^1.5,
            excess_kurtosis = m4 / m2^2 - 3)
end

function run_monte_carlo(data, Cov::AbstractMatrix{Float64};
        trials::Int = 1000,
        jacobian::Symbol = :numerical,
        seed::Int = 1234)
    N = data.N
    ε, θ, σ = 0.5, 0.001, 0.9   # paper's Monte-Carlo calibration (lines 90-92)
    Random.seed!(seed)          # reproducible draws (draw_logA uses the global RNG)

    loggdp = Float64[]
    n_converged = 0
    n_correct   = 0

    for k in 1:trials
        # independent sectoral shocks (diagonal Cov), mean-corrected E[A] = 1
        A_shock = exp.(draw_logA(Cov; n = 1)[:, 1])

        sol = solve_bf(A_shock, data.Ω, data.α, data.β, data.L, ε, θ, σ;
                       jacobian = jacobian)
        if sol.converged
            n_converged += 1
            g = real_gdp_mc(sol, A_shock, data.α, data.L, ε)
            if g > 0 && -0.4 < log(g) < 0.3   # MATLAB "correct" filter
                push!(loggdp, log(g))
                n_correct += 1
            end
        end
    end

    return (n_converged = n_converged, n_correct = n_correct,
            moments = moments_loggdp(loggdp), log_gdp = loggdp)
end
