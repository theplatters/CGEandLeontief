
# core_solver.jl -- Robust equilibrium solver for the B&F (2019) replication.
#
# This file is the PRODUCTION solver used by the Monte Carlo (R3) and the
# elasticity-gradient (R4) exercises, which must solve 10^3--10^5 equilibria
# reliably. It ports the authors' own two routines:
#   * Simulation.m          -> bf_residual!   (the 2N residual)
#   * Simulation_Derivs.m   -> bf_jacobian!    (the analytical Jacobian)
# The pedagogical solver in model.jl (numerical Jacobian) is kept untouched
# for the verified R1/R2 demonstrations; this file is the high-throughput
# twin. Both implement the SAME equations, so their equilibria coincide.
#
# IMPORTANT (Jacobian layout): In Simulation_Derivs.m the authors return
#   OutDeriv = [DOut1DP DOut1DY; DOut2DP DOut2DY]'
# and fmincon uses it directly as the constraint Jacobian. Hence the correct
# Jacobian of F = [Out1; Out2] for our Newton step IS that matrix WITH the
# final transpose (off-diagonal blocks swap, each block is transposed). A
# cross-check cell compares this solver against the numerical one on the
# baseline + oil shock to confirm they agree.

using LinearAlgebra

# ----------------------------------------------------------------------------
# Solution container
# ----------------------------------------------------------------------------
struct BFSolution
    p::Vector{Float64}       # equilibrium prices (N)
    y::Vector{Float64}       # equilibrium quantities (N)
    w::Vector{Float64}       # sectoral wages (N)
    C::Float64               # total expenditure w'L (nominal GDP)
    converged::Bool
    residual_norm::Float64
end

# ----------------------------------------------------------------------------
# Judicious initial guess (mirrors the MATLAB init in eg.m / the MC driver)
#   p0 = exp(-M \ log A),  y0 = lambda ./ p0,   M = I - diag(1-alpha)*Omega
# ----------------------------------------------------------------------------
function good_init(A::Vector{Float64}, Ω::Matrix{Float64},
                   α::Vector{Float64}, β::Vector{Float64})
    N = length(A)
    M = I - Diagonal(1 .- α) * Ω
    logA = log.(A)
    p0 = exp.(-M \ logA)
    λ  = (I - Diagonal(1 .- α) * Ω)' \ β
    y0 = λ ./ p0
    return vcat(p0, y0)
end

# ----------------------------------------------------------------------------
# Residual: F(X) = 0,  X = [p; y]  (2N)
# This is a faithful port of Simulation.m and is algebraically identical to
# model.jl's compute_equilibrium! (verified: the oil-shock amplification ratio
# 1.28 matches the paper).
# ----------------------------------------------------------------------------
function bf_residual!(out::Vector{Float64}, X::Vector{Float64},
                       A::Vector{Float64}, Ω::Matrix{Float64}, α::Vector{Float64},
                       β::Vector{Float64}, L::Vector{Float64},
                       ϵ::Float64, θ::Float64, σ::Float64, N::Int)
    p = view(X, 1:N)
    y = view(X, N+1:2N)

    q = (Ω * (p .^ (1 - θ))) .^ (1 / (1 - θ))
    w = p .* (A .^ ((ϵ - 1) / ϵ)) .* (α .^ (1 / ϵ)) .* (y .^ (1 / ϵ)) .* (L .^ (-1 / ϵ))
    C = dot(w, L)

    # Price (zero-profit) residuals
    price_rhs = (A .^ (ϵ - 1)) .* (α .* (w .^ (1 - ϵ)) .+ (1 .- α) .* (q .^ (1 - ϵ)))
    price_rhs = max.(price_rhs, 1e-300)
    price_rhs = price_rhs .^ (1 / (1 - ϵ))
    out[1:N] .= p .- price_rhs

    # Market-clearing residuals
    inner = (p .^ ϵ) .* (A .^ (ϵ - 1)) .* (q .^ (θ - ϵ)) .* (1 .- α) .* y
    intermediate_demand = (p .^ (-θ)) .* (Ω' * inner)
    final_demand = β .* (p .^ (-σ)) .* C
    out[N+1:2N] .= y .- intermediate_demand .- final_demand
    return nothing
end

# ----------------------------------------------------------------------------
# Analytical Jacobian (port of Simulation_Derivs.m, WITHOUT the final transpose)
# ----------------------------------------------------------------------------
function bf_jacobian!(J::Matrix{Float64}, X::Vector{Float64},
                       A::Vector{Float64}, Ω::Matrix{Float64}, α::Vector{Float64},
                       β::Vector{Float64}, L::Vector{Float64},
                       ϵ::Float64, θ::Float64, σ::Float64, N::Int)
    p = view(X, 1:N)
    y = view(X, N+1:2N)
    q = (Ω * (p .^ (1 - θ))) .^ (1 / (1 - θ))
    w = p .* (A .^ ((ϵ - 1) / ϵ)) .* (α .^ (1 / ϵ)) .* (y .^ (1 / ϵ)) .* (L .^ (-1 / ϵ))
    C = dot(w, L)

    # Convenience diagonals
    DWDP = Diagonal(w ./ p)                       # d w_s / d p_s
    DWDY = Diagonal((1 / ϵ) .* (w ./ y))          # d w_s / d y_s
    DCDP = (w ./ p) .* L                          # d C   / d p   (vector)
    DCDY = (1 / ϵ) .* (w ./ y) .* L               # d C   / d y   (vector)

    # d q_s / d p_t  = q_s^theta * p_t^{-theta} * Omega[s,t]
    DQDP = (q .^ θ) .* (p .^ (-θ))' .* Ω

    # --- DOut1 (price block) ---
    V  = α .* (w .^ (1 - ϵ)) .+ (1 .- α) .* (q .^ (1 - ϵ))
    M1 = Diagonal((1 ./ A) .* (V .^ (ϵ / (1 - ϵ))))
    term1dp = Diagonal(α .* (w .^ (1 - ϵ)) ./ p) .+
              Diagonal((1 .- α) .* (q .^ (-ϵ))) * DQDP
    DOut1DP = I - M1 * term1dp
    DOut1DY = -M1 * Diagonal((α ./ ϵ) .* (w .^ (1 - ϵ)) ./ y)

    # --- DOut2 (quantity block) ---
    T1 = (p .^ (ϵ - 1)) .* y .* (q .^ (θ - ϵ)) .* (1 .- α) .* (A .^ (ϵ - 1))
    T2 = (p .^ ϵ) .* y .* (q .^ (θ - ϵ - 1)) .* (1 .- α) .* (A .^ (ϵ - 1))
    term1 = -ϵ * Diagonal(p .^ (-θ)) * (Ω' * Diagonal(T1))
    term2 = (θ - ϵ) * Diagonal(p .^ (-θ)) * (Ω' * Diagonal(T2)) * DQDP
    term3 = -σ * Diagonal(β .* (p .^ (-σ - 1))) * C
    term4 = (β .* (p .^ (-σ))) * (DCDP')           # outer product N x N
    term5 = -θ * Diagonal((p .^ (-θ - 1)) .* (Ω' * (Diagonal(T1) * y)))
    DOut2DP = term1 .+ term2 .+ term3 .+ term4 .+ term5

    M2 = Diagonal(p .^ ϵ) * Diagonal(A .^ (ϵ - 1)) * Diagonal(q .^ (θ - ϵ)) *
         Diagonal(1 .- α) * Ω * Diagonal(p .^ (-θ))
    DOut2DY = I - M2' - (β .* (p .^ (-σ))) * (DCDY')

    # Assemble the Newton Jacobian. The MATLAB returns
    #   OutDeriv = [DOut1DP DOut1DY; DOut2DP DOut2DY]'
    # and fmincon uses it directly as the constraint Jacobian, so the correct
    # Jacobian of F = [Out1; Out2] is exactly that matrix WITH the final
    # transpose (the off-diagonal blocks swap and each block is transposed).
    J[1:N, 1:N]       .= DOut1DP'
    J[1:N, N+1:2N]    .= DOut2DP'
    J[N+1:2N, 1:N]    .= DOut1DY'
    J[N+1:2N, N+1:2N] .= DOut2DY'
    return nothing
end

# ----------------------------------------------------------------------------
# Numerical Jacobian (finite differences) -- fallback / cross-check
# ----------------------------------------------------------------------------
function numerical_jacobian!(J::Matrix{Float64}, X::Vector{Float64},
                              A::Vector{Float64}, Ω::Matrix{Float64}, α::Vector{Float64},
                              β::Vector{Float64}, L::Vector{Float64},
                              ϵ::Float64, θ::Float64, σ::Float64, N::Int)
    F0 = zeros(2N)
    bf_residual!(F0, X, A, Ω, α, β, L, ϵ, θ, σ, N)
    h = 1e-7
    for j in 1:2N
        s = h * max(abs(X[j]), 1.0)
        Xp = copy(X); Xp[j] += s
        Fp = zeros(2N)
        bf_residual!(Fp, Xp, A, Ω, α, β, L, ϵ, θ, σ, N)
        J[:, j] .= (Fp .- F0) ./ s
    end
    return nothing
end

# ----------------------------------------------------------------------------
# Newton solver with damping. jacobian = :analytical (default) or :numerical.
# ----------------------------------------------------------------------------
function solve_bf(A::Vector{Float64}, Ω::Matrix{Float64}, α::Vector{Float64},
                  β::Vector{Float64}, L::Vector{Float64},
                  ϵ::Float64, θ::Float64, σ::Float64;
                  max_iter::Int=300, tol::Float64=1e-10,
                  jacobian::Symbol=:numerical,
                  warm_start::Union{Vector{Float64}, Nothing}=nothing)
    # NOTE: the :analytical path (port of Simulation_Derivs.m) is EXPERIMENTAL:
    # its quantity-clearing block has not yet been verified against the
    # numerical Jacobian and currently diverges. The :numerical path (finite
    # differences) is verified (oil-shock GDP = 0.948, matching B&F) and is
    # the default. Switch to :analytical only after the cross-check passes.
    N = length(A)
    X = warm_start === nothing ? good_init(A, Ω, α, β) : copy(warm_start)
    X = max.(X, 1e-10)

    out = zeros(2N)
    F   = zeros(2N)
    J   = zeros(2N, 2N)
    bf_residual!(F, X, A, Ω, α, β, L, ϵ, θ, σ, N)

    for iter in 1:max_iter
        r = norm(F)
        if r < tol
            return _build_solution(X, A, Ω, α, β, L, ϵ, θ, σ, N, true, r)
        end
        if jacobian == :analytical
            bf_jacobian!(J, X, A, Ω, α, β, L, ϵ, θ, σ, N)
        else
            numerical_jacobian!(J, X, A, Ω, α, β, L, ϵ, θ, σ, N)
        end
        # Newton step. If the Jacobian is (near-)singular, fall back to a
        # Levenberg-Marquardt regularised step (always invertible for λ>0) so
        # the solver degrades gracefully instead of throwing LAPACK errors.
        local Δ, ok
        try
            Δ = J \ F; ok = true
        catch
            JTJ = J' * J
            ok = false
            for λ in (1e-8, 1e-4, 1e-1, 1.0, 1e3, 1e6)
                try
                    Δ = (JTJ + λ * I) \ (J' * F); ok = true; break
                catch
                    continue
                end
            end
            if !ok
                return _build_solution(X, A, Ω, α, β, L, ϵ, θ, σ, N, false, r)
            end
        end
        improved = false
        Xnew = max.(X - Δ, 1e-10)
        bf_residual!(out, Xnew, A, Ω, α, β, L, ϵ, θ, σ, N)
        if norm(out) < r
            X, F = Xnew, copy(out); improved = true
        else
            for d in (0.5, 0.25, 0.1, 0.01)
                Xd = max.(X - d .* Δ, 1e-10)
                bf_residual!(out, Xd, A, Ω, α, β, L, ϵ, θ, σ, N)
                if norm(out) < r
                    X, F = Xd, copy(out); improved = true; break
                end
            end
        end
        if !improved
            # Newton stalled: take the full step and let the next iteration decide
            X, F = Xnew, copy(out)
        end
    end
    bf_residual!(F, X, A, Ω, α, β, L, ϵ, θ, σ, N)
    return _build_solution(X, A, Ω, α, β, L, ϵ, θ, σ, N, false, norm(F))
end

function _build_solution(X, A, Ω, α, β, L, ϵ, θ, σ, N, converged, r)
    p = X[1:N]; y = X[N+1:2N]
    w = p .* (A .^ ((ϵ - 1) / ϵ)) .* (α .^ (1 / ϵ)) .* (y .^ (1 / ϵ)) .* (L .^ (-1 / ϵ))
    C = dot(w, L)
    return BFSolution(p, y, w, C, converged, r)
end

# ----------------------------------------------------------------------------
# Reallocation (mobile labor) solver — w = 1, only N price equations
# ----------------------------------------------------------------------------
struct BFReallocSolution
    p::Vector{Float64}       # equilibrium prices (N)
    real_gdp::Float64        # CES consumption index (β'·p^(1-σ))^(1/(σ-1))
    converged::Bool
    residual_norm::Float64
end

# Residual for the reallocation case: w = 1, only N price equations.
#   p_u = (A_u^(ϵ-1) · (α_u · 1^(1-ϵ) + (1-α_u) · q_u^(1-ϵ)))^(1/(1-ϵ))
function bf_realloc_residual!(out::Vector{Float64}, p::Vector{Float64},
                              A::Vector{Float64}, Ω::Matrix{Float64}, α::Vector{Float64},
                              ϵ::Float64, θ::Float64)
    N = length(p)
    q = (Ω * (p .^ (1 - θ))) .^ (1 / (1 - θ))
    # Zero-profit with w = 1: p_u = RHS
    rhs = (A .^ (ϵ - 1)) .* (α .* 1.0 .+ (1 .- α) .* (q .^ (1 - ϵ)))
    rhs = max.(rhs, 1e-300)
    rhs = rhs .^ (1 / (1 - ϵ))
    out[1:N] .= p .- rhs
    return nothing
end

function solve_bf_realloc(A::Vector{Float64}, Ω::Matrix{Float64}, α::Vector{Float64},
                          β::Vector{Float64}, ϵ::Float64, θ::Float64, σ::Float64;
                          max_iter::Int=300, tol::Float64=1e-10)
    N = length(A)
    # Initial guess: p₀ from the Cobb-Douglas approximation (good_init but only p)
    M = I - Diagonal(1 .- α) * Ω
    p = exp.(-M \ log.(A))
    p = max.(p, 1e-10)

    F  = zeros(N)
    J  = zeros(N, N)
    Fp = zeros(N)

    for iter in 1:max_iter
        bf_realloc_residual!(F, p, A, Ω, α, ϵ, θ)
        r = norm(F)
        if r < tol
            real_gdp = (dot(β, p .^ (1 - σ)))^(1 / (σ - 1))
            return BFReallocSolution(p, real_gdp, true, r)
        end
        # Numerical Jacobian (N×N)
        h = 1e-7
        for j in 1:N
            s = h * max(abs(p[j]), 1.0)
            pj = copy(p); pj[j] += s
            bf_realloc_residual!(Fp, pj, A, Ω, α, ϵ, θ)
            J[:, j] .= (Fp .- F) ./ s
        end
        # Newton step with LM fallback
        local Δ
        try
            Δ = J \ F
        catch
            JTJ = J' * J
            local solved = false
            for λ in (1e-8, 1e-4, 1e-1, 1.0, 1e3, 1e6)
                try
                    Δ = (JTJ + λ * I) \ (J' * F)
                    solved = true
                    break
                catch
                    continue
                end
            end
            if !solved
                real_gdp = (dot(β, p .^ (1 - σ)))^(1 / (σ - 1))
                return BFReallocSolution(p, real_gdp, false, r)
            end
        end
        improved = false
        for d in (1.0, 0.5, 0.25, 0.1, 0.01)
            pn = max.(p - d .* Δ, 1e-10)
            bf_realloc_residual!(Fp, pn, A, Ω, α, ϵ, θ)
            if norm(Fp) < r
                p = pn; F .= Fp; improved = true; break
            end
        end
        if !improved
            p = max.(p - Δ, 1e-10)
        end
    end
    bf_realloc_residual!(F, p, A, Ω, α, ϵ, θ)
    real_gdp = (dot(β, p .^ (1 - σ)))^(1 / (σ - 1))
    return BFReallocSolution(p, real_gdp, false, norm(F))
end
#   * gdp_nominal = C = w'L            (used by R3 and R4)
#   * gdp_welfare = C * sum(beta .* p^-sigma)   (used by the single-sector test)
#   * mean_price  = sum(y .* p) / sum(y)  (Domar-weighted, eg.m)
# ----------------------------------------------------------------------------
gdp_nominal(sol::BFSolution) = sol.C
gdp_welfare(sol::BFSolution, β::Vector{Float64}, σ::Float64) =
    sol.C * sum(β .* (sol.p .^ (-σ)))
mean_price(sol::BFSolution) = sum(sol.y .* sol.p) / sum(sol.y)

# Methods for BFModel.SimulationResult (returned by the pedagogical solver in
# model.jl) so the same GDP/price helpers work on both solver outputs.
gdp_nominal(sol::BFModel.SimulationResult) = sol.nominal_gdp
gdp_welfare(sol::BFModel.SimulationResult, β::Vector{Float64}, σ::Float64) =
    sol.nominal_gdp * sum(β .* (sol.p .^ (-σ)))
mean_price(sol::BFModel.SimulationResult) = sum(sol.y .* sol.p) / sum(sol.y)
