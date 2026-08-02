"""
    model.jl — Equilibrium system for B&F (2019) replication

Ports the MATLAB `Simulation.m` function to Julia.

The B&F model is a one-factor (labor) CES economy with N sectors.
The equilibrium is defined by a 2N system of equations in:
  - Prices p (N equations: zero-profit / price normalization)
  - Quantities y (N equations: market clearing)

Key equations (matching MATLAB notation):

  Price index for intermediates:
    q_u = (Σ_s Ω_us · p_s^(1-θ))^(1/(1-θ))

  Sectoral wage (from labor FOC):
    w_u = p_u · A_u^((ϵ-1)/ϵ) · α_u^(1/ϵ) · y_u^(1/ϵ) · L_u^(-1/ϵ)

  Total consumption expenditure:
    C = w' · L

  Zero-profit condition (price equations):
    p_u = (A_u^(ϵ-1) · (α_u · w_u^(1-ϵ) + (1-α_u) · q_u^(1-ϵ)))^(1/(1-ϵ))

  Market clearing (quantity equations):
    y_u = Σ_s Ω_su · p_s^ϵ · A_s^(ϵ-1) · q_s^(θ-ϵ) · (1-α_s) · y_s · p_u^(-θ)
          + β_u · p_u^(-σ) · C

  CPI (consumption price index):
    CPI = (Σ_u β_u · p_u^(1-σ))^(1/(1-σ))

Note: In the MATLAB code, CPI is NOT explicitly computed inside Simulation.m.
The consumption demand term β_u · p_u^(-σ) · C uses total expenditure C = w'·L
without a CPI normalization. The Julia BeyondHulten package adds CPI explicitly.
"""
module BFModel

using LinearAlgebra
using Statistics

export BFParameters, compute_equilibrium, compute_gdp, compute_cpi,
       compute_mean_price, SimulationResult

"""
    BFParameters

Parameters for the B&F (2019) CES equilibrium.
"""
struct BFParameters
    A::Vector{Float64}       # productivity vector (N)
    Ω::Matrix{Float64}      # cost share matrix (N×N)
    α::Vector{Float64}      # factor shares (N)
    β::Vector{Float64}      # consumption shares (N)
    L::Vector{Float64}      # labor allocation (N)
    ϵ::Float64              # elasticity: intermediates vs labor
    θ::Float64              # elasticity: between intermediates
    σ::Float64              # elasticity: consumption substitution
    N::Int                  # number of sectors
end

function BFParameters(A, Ω, α, β, L, ϵ, θ, σ)
    N = length(α)
    @assert length(A) == N "A must have N elements"
    @assert size(Ω) == (N, N) "Ω must be N×N"
    @assert length(β) == N "β must have N elements"
    @assert length(L) == N "L must have N elements"
    BFParameters(A, Ω, α, β, L, ϵ, θ, σ, N)
end

"""
    SimulationResult

Result of a simulation: equilibrium prices, quantities, and derived quantities.
"""
struct SimulationResult
    p::Vector{Float64}       # equilibrium prices (N)
    y::Vector{Float64}       # equilibrium quantities (N)
    w::Vector{Float64}       # sectoral wages (N)
    q::Vector{Float64}       # intermediate price indices (N)
    C::Float64               # total consumption expenditure
    cpi::Float64             # consumption price index
    nominal_gdp::Float64     # w'·L (nominal GDP = total wage bill)
    real_gdp::Float64        # (β'·p^(1-σ))^(1/(σ-1)) — CES real consumption
    converged::Bool          # whether solver converged
    residual_norm::Float64   # ||F(x)||
end

"""
    compute_equilibrium!(out, X, params)

Compute the 2N residual vector for the equilibrium system.

  out[1:N]   = price residuals (zero-profit conditions)
  out[N+1:2N] = quantity residuals (market clearing)

This is a direct port of Simulation.m.
"""
function compute_equilibrium!(out::Vector{Float64}, X::Vector{Float64},
                               params::BFParameters)
    (; A, Ω, α, β, L, ϵ, θ, σ, N) = params

    # Split unknown vector
    p = view(X, 1:N)
    y = view(X, N+1:2*N)

    # Intermediate price index: q_u = (Σ_s Ω_us · p_s^(1-θ))^(1/(1-θ))
    q = (Ω * (p .^ (1 - θ))) .^ (1 / (1 - θ))

    # Sectoral wage: w_u = p_u · A_u^((ϵ-1)/ϵ) · α_u^(1/ϵ) · y_u^(1/ϵ) · L_u^(-1/ϵ)
    w = p .* (A .^ ((ϵ - 1) / ϵ)) .* (α .^ (1 / ϵ)) .* (y .^ (1 / ϵ)) .* (L .^ (-1 / ϵ))

    # Total consumption expenditure
    C = dot(w, L)

    # Price equations: p_u = (A_u^(ϵ-1) · (α_u·w_u^(1-ϵ) + (1-α_u)·q_u^(1-ϵ)))^(1/(1-ϵ))
    # Residual: p_u - RHS = 0
    price_rhs = (A .^ (ϵ - 1)) .* (α .* (w .^ (1 - ϵ)) + (1 .- α) .* (q .^ (1 - ϵ)))
    # Handle potential negative values under the root
    price_rhs = max.(price_rhs, 1e-300)
    price_rhs = price_rhs .^ (1 / (1 - ϵ))

    out[1:N] .= p .- price_rhs

    # Market clearing: y_u = intermediates_demand + final_demand
    # intermediates_demand = Σ_s Ω_su · p_s^ϵ · A_s^(ϵ-1) · q_s^(θ-ϵ) · (1-α_s) · y_s · p_u^(-θ)
    # In matrix form: diag(p)^(-θ) · Ω' · diag(p^ϵ · A^(ϵ-1) · q^(θ-ϵ) · (1-α) · y)

    # MATLAB: y' * diag(p)^epsilon * diag(A)^(epsilon-1) * diag(q)^(theta-epsilon)
    #         * diag(1-alpha) * Omega * diag(p)^(-theta)
    # Note: MATLAB's Omega is the cost share matrix (row-normalized IO)
    # The market clearing equation uses Ω' (transpose)

    # intermediate_demand vector:
    # id_u = p_u^(-θ) · Σ_s Ω_su · p_s^ϵ · A_s^(ϵ-1) · q_s^(θ-ϵ) · (1-α_s) · y_s
    inner = (p .^ ϵ) .* (A .^ (ϵ - 1)) .* (q .^ (θ - ϵ)) .* (1 .- α) .* y
    intermediate_demand = (p .^ (-θ)) .* (Ω' * inner)

    # Final demand: β_u · p_u^(-σ) · C
    final_demand = β .* (p .^ (-σ)) .* C

    out[N+1:2N] .= y .- intermediate_demand .- final_demand

    nothing
end

"""
    compute_equilibrium(params; init) -> SimulationResult

Solve the equilibrium system using NonlinearSolve.
"""
function compute_equilibrium(params::BFParameters; init::Vector{Float64}=Float64[])
    (; A, Ω, α, β, L, ϵ, θ, σ, N) = params

    # Default initial guess: prices=1, quantities=λ (Domar weights)
    if isempty(init)
        # MATLAB: init = [ones(N,1); lambda]
        # lambda = (I - diag(1-alpha)*Omega)' \ beta
        λ = (I - diagm(0 => 1 .- α) * Ω)' \ β
        init = [ones(N); λ]
    end

    # Smart initial guess from MATLAB:
    # init = [exp(-(I-diag(1-alpha)*Omega) \ log(A)); lambda ./ exp(...)]
    # This uses the Cobb-Douglas approximation as starting point
    if all(A .≈ 1.0)
        # Baseline: just use unit prices and Domar weights
    else
        try
            log_A = log.(A)
            p_init = exp.(-((I - diagm(0 => 1 .- α) * Ω) \ log_A))
            λ = (I - diagm(0 => 1 .- α) * Ω)' \ β
            y_init = λ ./ p_init
            init = [p_init; y_init]
        catch
            # Fall back to default init if matrix is singular
        end
    end

    # Ensure init is positive
    init = max.(init, 1e-10)

    # Define the residual function
    F!(out, x, p) = compute_equilibrium!(out, x, p)

    # Use simple Newton's method (avoids NonlinearSolve dependency issues)
    result = _newton_solve(F!, init, params)

    return result
end

"""
    _newton_solve(F!, x0, params; max_iter, tol)

Simple Newton-Raphson solver with numerical Jacobian.
This avoids the need for NonlinearSolve.jl and matches the
fmincon approach in the MATLAB code (which uses finite differences).
"""
function _newton_solve(F!::Function, x0::Vector{Float64}, params;
                        max_iter::Int=200, tol::Float64=1e-10)
    n = length(x0)
    x = copy(x0)
    F = zeros(n)
    J = zeros(n, n)
    x_pert = zeros(n)

    for iter in 1:max_iter
        F!(F, x, params)
        norm_F = norm(F)

        if norm_F < tol
            # Converged! Compute result
            return _build_result(x, params, true, norm_F)
        end

        # Numerical Jacobian (forward differences)
        h = 1e-7
        for j in 1:n
            x_pert .= x
            x_pert[j] += h * max(abs(x[j]), 1.0)
            F_pert = zeros(n)
            F!(F_pert, x_pert, params)
            J[:, j] .= (F_pert .- F) / (h * max(abs(x[j]), 1.0))
        end

        # Newton step: x_new = x - J⁻¹ F
        try
            Δx = J \ F
            # Damping factor for stability
            α_damp = 1.0
            x_new = x - α_damp .* Δx

            # Enforce positivity (prices and quantities must be positive)
            x_new = max.(x_new, 1e-10)

            # Check if step is productive
            F_new = zeros(n)
            F!(F_new, x_new, params)
            if norm(F_new) < norm_F
                x = x_new
            else
                # Try damped step
                for damping in [0.5, 0.25, 0.1, 0.01]
                    x_damped = x - damping .* Δx
                    x_damped = max.(x_damped, 1e-10)
                    F_damped = zeros(n)
                    F!(F_damped, x_damped, params)
                    if norm(F_damped) < norm_F
                        x = x_damped
                        break
                    end
                end
                # If nothing worked, take full step anyway
                if norm(F_new) >= norm_F
                    x = max.(x - Δx, 1e-10)
                end
            end
        catch e
            # Jacobian is singular — try pseudo-inverse or small perturbation
            @warn "Jacobian singular at iter $iter, perturbing"
            x = x .* (1 .+ 0.001 .* randn(n))
            x = max.(x, 1e-10)
        end
    end

    # Did not converge
    F!(F, x, params)
    norm_F = norm(F)
    @warn "Solver did not converge. Final residual: $norm_F"
    return _build_result(x, params, false, norm_F)
end

"""
    _build_result(x, params, converged, residual_norm)

Build a SimulationResult from the solution vector.
"""
function _build_result(x::Vector{Float64}, params::BFParameters,
                        converged::Bool, residual_norm::Float64)
    (; A, Ω, α, β, L, ϵ, θ, σ, N) = params

    p = x[1:N]
    y = x[N+1:2N]

    q = (Ω * (p .^ (1 - θ))) .^ (1 / (1 - θ))
    w = p .* (A .^ ((ϵ - 1) / ϵ)) .* (α .^ (1 / ϵ)) .* (y .^ (1 / ϵ)) .* (L .^ (-1 / ϵ))
    C = dot(w, L)

    cpi_val = compute_cpi(p, β, σ)
    nominal_gdp = C  # w'·L
    real_gdp = (dot(β, p .^ (1 - σ)))^(1 / (σ - 1))  # = (β'·p^(1-σ))^(1/(σ-1))

    # Note: in MATLAB, GDP for the no-reallocation case is:
    #   GDP = (w' * L) = nominal GDP
    # For the reallocation case:
    #   GDP = (β' * p^(1-σ))^(1/(σ-1)) = real GDP (CES consumption index)
    # Both are relative to baseline (=1 when A=1).

    return SimulationResult(p, y, w, q, C, cpi_val, nominal_gdp, real_gdp,
                             converged, residual_norm)
end

"""
    compute_cpi(p, β, σ) -> Float64

Consumer Price Index: CPI = (Σ_u β_u · p_u^(1-σ))^(1/(1-σ))
"""
function compute_cpi(p::Vector{Float64}, β::Vector{Float64}, σ::Float64)
    (dot(β, p .^ (1 - σ)))^(1 / (1 - σ))
end

"""
    compute_gdp(result::SimulationResult; type=:nominal) -> Float64

GDP computation matching the MATLAB code.

  - type=:nominal → w'·L (total wage bill, MATLAB's GDP for fixed labor)
  - type=:real → (β'·p^(1-σ))^(1/(σ-1)) (CES consumption index, MATLAB's GDP for mobile labor)
  - type=:welfare → Σ(w·L·p^(-σ)·β) (welfare-based GDP from eg.m)
"""
function compute_gdp(result::SimulationResult; type::Symbol=:nominal)
    if type == :nominal
        return result.nominal_gdp
    elseif type == :real
        return result.real_gdp
    elseif type == :welfare
        # From eg.m: sum(wages' * L * p.^(-sigma) .* beta)
        (; p, w, C) = result
        # This requires access to β, σ, L — but they're not in SimulationResult
        # Would need to pass params. For now, use stored values.
        error("Welfare GDP requires parameters. Use compute_gdp(result, params) instead.")
    end
end

"""
    compute_mean_price(result::SimulationResult) -> Float64

Mean (Domar-weighted) price: Σ(y·p) / Σ(y)

This matches the MATLAB eg.m function:
  mean_prices = sum(Soln(N+1:end) .* Soln(1:N)) / sum(Soln(N+1:end))
"""
function compute_mean_price(result::SimulationResult)
    dot(result.y, result.p) / sum(result.y)
end

end # module
