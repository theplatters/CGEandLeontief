
# shocks.jl -- TFP shock generation for the Monte Carlo (R3).
#
# Faithful port of the MATLAB data pipeline in GDP_Simulation_88sectorKLEMS.m /
# GDP_Simulation_88sectorKLEMS_JK.m:
#
#   * stfp (sectoral TFP growth, N x T) -> covariance matrix
#   * 4-year cumulative covariance (Sigma_4year) exactly as the JK variant:
#         cum_stfp      = cumsum(log(1 + stfp), dims=2)
#         cum_stfp_4yr  = cum_stfp[:, 1:4:end]
#         Sigma_4year   = cov(diff(cum_stfp_4yr, dims=2)')
#   * draws:  A = exp(z),  z ~ N( mu, Cov ),  mu = -1/2 * diag(Cov)
#     so that E[A] = 1. The covariance used is the DIAGONAL (independent
#     sectoral shocks), matching the active line in the MATLAB code; set
#     correlated=true to draw from the full covariance instead.

using LinearAlgebra, Statistics

# Local covariance: M is d x n (d variables, n observations) -> d x d.
function mycov(M::Matrix{Float64})
    n = size(M, 2)
    X = M .- mean(M, dims = 2)
    return (X * X') ./ (n - 1)
end

# Empirical covariances from the cleaned stfp matrix (N x T).
# Returns (Sigma_yearly, Sigma_4year). Matches MATLAB cov(stfp') and the
# 4-year construction above.
function empirical_covariances(stfp::Matrix{Float64})
    Sigma_yearly = mycov(stfp)                       # N x N
    cum_stfp     = cumsum(log.(1 .+ stfp); dims = 2)
    cum_stfp_4yr = cum_stfp[:, 1:4:end]
    Sigma_4year  = mycov(diff(cum_stfp_4yr; dims = 2))
    return Sigma_yearly, Sigma_4year
end

# Multivariate normal sampler (port of mvnrnd.m: Cholesky with an
# eigendecomposition fallback). Returns an n x d matrix of draws (rows).
function mvnrnd(mu::Vector{Float64}, Sigma::Matrix{Float64}; n::Int = 1)
    d = length(mu)
    try
        U = cholesky(Sigma).U
        return randn(n, d) * U .+ mu'
    catch
        F = eigen(Sigma)
        if minimum(F.values) < 0
            error("Sigma must be positive semi-definite.")
        end
        U = Diagonal(sqrt.(max.(F.values, 0.0))) * F.vectors'
        return randn(n, d) * U .+ mu'
    end
end

# Draw log-A vectors. Each COLUMN is one draw (N x n).
# Default: independent sectoral shocks (diagonal covariance) as in MATLAB.
function draw_logA(Sigma::Matrix{Float64}; n::Int = 1, correlated::Bool = false)
    mu  = -0.5 .* diag(Sigma)
    Cov = correlated ? Sigma : Diagonal(diag(Sigma))
    return mvnrnd(mu, Matrix(Cov); n = n)'   # N x n
end
