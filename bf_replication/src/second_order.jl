# second_order.jl — Second-order approximation (R6)
#
# Faithful port of Second_Order_Simulation.m.
#
# Two cases:
#   No-reallocation (fixed labor):
#     Hessian via finite differences of nominal GDP w.r.t. A_i
#     Formula: log GDP ≈ λ'·log(A) + ½·log(A)'·Hess·log(A)
#
#   Full-reallocation (mobile labor, w=1):
#     Hessian via finite differences of Domar weights λ = p·y/C
#     Same formula, different Hessian.
#
using LinearAlgebra, Random

# ---------------------------------------------------------------------------
# Baseline Domar weights
# ---------------------------------------------------------------------------
function _domar_weights(Ω, α, β)
    N = length(α)
    return (I - Diagonal(1 .- α) * Ω)' \ β
end

# ---------------------------------------------------------------------------
# Part 1: No-reallocation (fixed-labor) Hessian
#
# For each sector i, perturb A_i by h, solve fixed-labor equilibrium,
# compute GDP_nominal = w'L = gdp_nominal, return (GDP_pert - GDP_base) / h.
# ---------------------------------------------------------------------------
function second_order_hessian_norealloc(Ω, α, β, L, ε, θ, σ; h::Float64=1e-6)
    N = length(α)
    A0 = ones(N)
    gdp0 = gdp_nominal(solve_bf(A0, Ω, α, β, L, ε, θ, σ))

    H = zeros(N, N)
    for i in 1:N
        A_i = copy(A0); A_i[i] += h
        sol = solve_bf(A_i, Ω, α, β, L, ε, θ, σ)
        gdp_i = gdp_nominal(sol)
        d_gdp = (gdp_i - gdp0) / h
        # MATLAB: Hess(:,i) = scalar → all rows of column i get the same value
        H[:, i] .= d_gdp
    end
    return H
end

# ---------------------------------------------------------------------------
# Part 2: Full-reallocation (mobile labor) Hessian
#
# For each sector i, perturb A_i by h, call solve_bf_realloc to get prices,
# then compute real GDP C = (β'·p^(1-σ))^(1/(σ-1)), quantities y, and
# new Domar weights λ_new = p·y / C.
# Hess[:,i] = (λ_new - λ_base) / h  — change in Domar weights in response
# to a TFP perturbation.
# ---------------------------------------------------------------------------
function second_order_hessian_realloc(Ω, α, β, ε, θ, σ; h::Float64=1e-5)
    N = length(α)
    λ0 = _domar_weights(Ω, α, β)

    H = zeros(N, N)
    for i in 1:N
        A_i = ones(N); A_i[i] += h
        sol = solve_bf_realloc(A_i, Ω, α, β, ε, θ, σ)
        p = sol.p
        C = sol.real_gdp
        # q — intermediate price index
        q = (Ω * (p .^ (1 - θ))) .^ (1 / (1 - θ))
        # y from the market-clearing system:
        # y = (I - M') \ final_demand  where M = diag(p)^ε diag(A)^(ε-1) diag(q)^(θ-ε) diag(1-α) Ω diag(p)^(-θ)
        M = Diagonal(p .^ ε) * Diagonal(A_i .^ (ε - 1)) * Diagonal(q .^ (θ - ε)) *
            Diagonal(1 .- α) * Ω * Diagonal(p .^ (-θ))
        final_demand = (β' * Diagonal(p .^ (-σ))) * C
        # In MATLAB: y = ((beta'*diag(p)^(-sigma)*C)*inv(I - M)')'
        # This is: y = (I - M') \ final_demand
        y = (I - M') \ vec(final_demand)
        # New Domar weights
        λ_new = p .* y ./ C
        H[:, i] = (λ_new .- λ0) ./ h
    end
    return H
end

# ---------------------------------------------------------------------------
# Second-order approximation MC
#   log GDP ≈ 0 + λ'·log(A) + ½·log(A)'·H·log(A)
#
# Draws A from the same diagonal-normal distribution as the main MC and
# evaluates the approximation formula (no solver per draw — fast).
# ---------------------------------------------------------------------------
function second_order_mc(λ::Vector{Float64}, H::Matrix{Float64},
                          Cov::Matrix{Float64};
                          trials::Int=1000, seed::Int=12345)
    Random.seed!(seed)
    loggdp = zeros(trials)
    for k in 1:trials
        logA = draw_logA(Cov; n=1)[:, 1]  # N-vector of log shocks
        first_order = dot(λ, logA)
        second_order = 0.5 * dot(logA, H * logA)
        loggdp[k] = first_order + second_order
    end
    return moments_loggdp(loggdp)
end

# Convenience wrapper: load data once, compute Hessian, run MC.
function run_second_order_norealloc(data, Cov::Matrix{Float64};
        trials::Int=1000, seed::Int=12345,
        h::Float64=1e-6)
    ε, θ, σ = 0.5, 0.001, 0.9
    λ = _domar_weights(data.Ω, data.α, data.β)
    H = second_order_hessian_norealloc(data.Ω, data.α, data.β, data.L, ε, θ, σ; h=h)
    mc = second_order_mc(λ, H, Cov; trials=trials, seed=seed)
    return (λ=λ, H=H, mc=mc)
end

function run_second_order_realloc(data, Cov::Matrix{Float64};
        trials::Int=1000, seed::Int=12345,
        h::Float64=1e-5)
    ε, θ, σ = 0.5, 0.001, 0.9
    λ = _domar_weights(data.Ω, data.α, data.β)
    H = second_order_hessian_realloc(data.Ω, data.α, data.β, ε, θ, σ; h=h)
    mc = second_order_mc(λ, H, Cov; trials=trials, seed=seed)
    return (λ=λ, H=H, mc=mc)
end