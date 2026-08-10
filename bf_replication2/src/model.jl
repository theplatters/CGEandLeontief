# src/model.jl — Equilibrium system (MCP via Fischer–Burmeister + NLsolve)
#
# Ports Standard_form_Covid_HtM3_logCES.mod into a square system of
# nonlinear equations. The sticky-labor complementarity is replaced by
# the smooth Fischer–Burmeister function φ(x,y) = x + y - √(x² + y²).
#
# Solved with NLsolve (finite-difference Jacobian).

using LinearAlgebra, Printf, Statistics
using NLsolve

export MCPModel, make_theta, equilibrium_residual!, solve_equilibrium, eq_continuation

"""
    MCPModel

Parameters for the equilibrium system.

# Fields
- `D::Int`                   — total dimension = 5N + 4
- `N::Int`                   — number of sectors
- `Omega_re::Matrix`         — standard-form IO matrix (D×D)
- `factor::Vector{Int}`      — factor type codes
- `keynes::Vector{Int}`      — price rigidity codes
- `theta::Vector{Float64}`   — elasticities per sector (D)
- `cobb_douglas::Vector{Int}` — 0=CES, 1=CD, 2=no price eq
- `phi::Vector{Float64}`     — factor supply elasticity (D)
- `phi_htm::Vector{Float64}` — HtM social insurance (D)
- `mu::Vector{Float64}`      — markup (D, all 1)
- `in_mu::Vector{Float64}`   — inverse markup (D, all 1)
- `A::Vector{Float64}`       — supply productivity shock (D)
- `B::Matrix{Float64}`       — demand shock matrix (D×D)
- `init_lambda::Vector{Float64}` — initial inverse demands (Domar weights)
- `init_p::Vector{Float64}`  — initial prices (all 1)
- `chi::Vector{Float64}`     — factor ownership
"""
struct MCPModel
    D::Int
    N::Int
    Omega_re::Matrix{Float64}
    factor::Vector{Int}
    keynes::Vector{Int}
    theta::Vector{Float64}
    cobb_douglas::Vector{Int}
    phi::Vector{Float64}
    phi_htm::Vector{Float64}
    mu::Vector{Float64}
    in_mu::Vector{Float64}
    A::Vector{Float64}
    B::Matrix{Float64}
    init_lambda::Vector{Float64}
    init_p::Vector{Float64}
    chi::Vector{Float64}
end

"""
    make_theta(N; sigma=1.0, epsilon=0.6, eta=0.6, theta1=0.2, factor_ela=0.2, rho=1.0) -> Vector{Float64}

Build the D-vector of sector-specific elasticities, matching Master_file_3.m.

# Layout
- θ[1]              = sigma          (consumption)
- θ[2 .. N+1]       = epsilon        (VA vs intermediates)
- θ[N+2 .. 2N+1]    = eta            (across VA)
- θ[2N+2 .. 3N+1]   = theta1         (across intermediates)
- θ[3N+2 .. 5N+1]   = factor_ela     (labor + capital)
- θ[5N+2]           = 0.95           (HtM)
- θ[5N+3]           = rho            (Ricardian)
- θ[5N+4]           = 0.95           (consumption tomorrow)
"""
function make_theta(N; sigma=1.0, epsilon=0.6, eta=0.6, theta1=0.2, factor_ela=0.2, rho=1.0)
    D = 5 * N + 4
    θ = zeros(D)
    i = 1
    θ[i] = sigma; i += 1
    θ[i:(i+N-1)] .= epsilon; i += N
    θ[i:(i+N-1)] .= eta; i += N
    θ[i:(i+N-1)] .= theta1; i += N
    θ[i:(i+2N-1)] .= factor_ela; i += 2N
    θ[i] = 0.95; i += 1
    θ[i] = rho; i += 1
    θ[i] = 0.95
    return θ
end

"""
    make_cobb_douglas(N, D, factor; benchmark=true) -> Vector{Int}

Build the cobb_douglas indicator vector.

# Benchmark regime (complementarity, Fig 2)
- Goods (factor=1): CES (cobb_douglas=0)
- Row 1 (consumption): CD (cobb_douglas=1)
- Row 5N+3 (Ricardian): CD (cobb_douglas=1)
- Factors and consumers: 2 (no price equation)

# CD regime (Fig 3)
- All goods: CD (cobb_douglas=1)
- Factors: 2 (no price equation)
"""
function make_cobb_douglas(N, D, factor; benchmark=true)
    cd_vec = 2 * ones(Int, D)
    if benchmark
        # Goods (factor=1) are CES
        for k in 1:D
            if factor[k] == 1
                cd_vec[k] = 0
            end
        end
        cd_vec[1] = 1              # consumption today = CD
        cd_vec[5*N+3] = 1          # Ricardian consumer = CD
    else
        # CD regime: everything with factor>0 is CD
        for k in 1:D
            if factor[k] == 1
                cd_vec[k] = 1
            end
        end
        cd_vec[1] = 1
        cd_vec[5*N+3] = 1
        # Factors remain 2
    end
    return cd_vec
end

"""
    equilibrium_residual!(F, z, m::MCPModel)

Fill the residual vector F (length 2D) for the equilibrium system.
Variables z = [p[1:D]; λ[1:D]].

Residual block layout:
- F[1:D]        — price equations (CES/CD for goods; factor clearing for factors;
                   FB for sticky factors)
- F[D+1:2D]     — lambda equations (propagation for goods/factors;
                   HtM/Ricardian budget constraints)
- Special: p[D] - 1 = 0 (numeraire at F[D])
"""
function equilibrium_residual!(F, z, m::MCPModel)
    D = m.D
    p = @view z[1:D]
    λ = @view z[D+1:2D]

    # Clamp to avoid numerical issues
    # p must be strictly positive (prices)
    for i in 1:D
        if p[i] ≤ 1e-15
            p[i] = 1e-15
        end
    end
    for i in 1:D
        if λ[i] ≤ 1e-15
            λ[i] = 1e-15
        end
    end

    Ω = m.Omega_re
    θ = m.theta
    cd = m.cobb_douglas
    μ = m.mu
    in_μ = m.in_mu
    φ = m.phi
    φ_htm = m.phi_htm
    A = m.A
    B = m.B
    λ̄ = m.init_lambda
    f = m.factor
    kyn = m.keynes

    # ──────────────────────────────────────────────
    # 1. Price equations (F[1:D])
    # ──────────────────────────────────────────────
    for k in 1:D
        if kyn[k] == 0
            # Factor clearing — flexible (capital, tomorrow)
            # λ_k = p_k · (p_k / p_1)^φ_k · λ̄_k · A_k
            floor_mult = (p[k] / p[1])^φ[k]
            F[k] = λ[k] - p[k] * floor_mult * λ̄[k] * A[k]

        elseif kyn[k] == -1
            # Factor clearing — sticky labor, FB reformulation
            # x = p - floor ≥ 0
            # y = A·λ̄ - λ/p ≥ 0
            # φ(x,y) = x + y - √(x² + y²) = 0
            floor_val = ((λ[k] / p[k]) / (A[k] * λ̄[k]))^φ[k]
            if φ[k] == 0
                floor_val = 1.0
            end
            x = p[k] - floor_val
            y = A[k] * λ̄[k] - λ[k] / p[k]
            # FB function
            F[k] = x + y - sqrt(x^2 + y^2)

        elseif cd[k] == 1
            # Cobb-Douglas price equation
            # log(p_k) = -log(A_k) + Σ_j B_kj · Ω_kj · log(p_j)
            sum_val = 0.0
            for j in 1:D
                if Ω[k, j] > 0
                    sum_val += B[k, j] * Ω[k, j] * log(p[j])
                end
            end
            F[k] = log(p[k]) + log(A[k]) - sum_val

        elseif cd[k] == 0
            # CES price equation
            # p_k = (μ/in_μ)/A_k · (Σ_j B_kj^θ · Ω_kj · p_j^(1-θ))^(1/(1-θ))
            θk = θ[k]
            if abs(θk - 1.0) < 1e-12
                # Near-unit elasticity: fall back to CD-like behavior
                sum_val = 0.0
                for j in 1:D
                    if Ω[k, j] > 0
                        sum_val += B[k, j] * Ω[k, j] * log(p[j])
                    end
                end
                F[k] = log(p[k]) + log(A[k]) - sum_val
            else
                sum_val = 0.0
                for j in 1:D
                    if Ω[k, j] > 0
                        sum_val += (B[k, j]^θk) * Ω[k, j] * (p[j]^(1 - θk))
                    end
                end
                scale = (μ[k] / in_μ[k]) / A[k]
                F[k] = p[k] - scale * (sum_val)^(1 / (1 - θk))
            end

        else
            # cd[k] == 2: no price equation (e.g. HtM row)
            # Leave this variable free
            F[k] = 0.0
        end
    end

    # ──────────────────────────────────────────────
    # 2. Lambda equations (F[D+1:2D])
    # ──────────────────────────────────────────────
    for k in 1:D
        if f[k] < 2
            # Standard lambda propagation
            # λ_k = Σ_{j: f[j]>0} λ_j · μ_j⁻¹ · B_jk^θ · Ω_jk · (p_k/p_j)^(1-θ) · (μ_j/in_μ_j)^(1-θ) · A_j^(θ-1)
            sum_val = 0.0
            for j in 1:D
                if f[j] > 0
                    θj = θ[j]
                    term = λ[j] / μ[j]
                    if Ω[j, k] > 0
                        term *= (B[j, k]^θj) * Ω[j, k]
                        term *= (p[k] / p[j])^(1 - θj)
                        term *= (μ[j] / in_μ[j])^(1 - θj)
                        term *= A[j]^(θj - 1)
                        sum_val += term
                    end
                end
            end
            F[D + k] = λ[k] - sum_val

        elseif f[k] == 3
            # HtM consumer lambda (row 5N+2)
            # λ_htm = (Σ_{j: f[j]=0} λ̄_j · (1-φ_htm_j) · (1 - λ_j/p_j/λ̄_j)) · λ_last
            income_sum = 0.0
            for j in 1:D
                if f[j] == 0
                    income_sum += λ̄[j] * (1 - φ_htm[j]) * (1 - λ[j] / p[j] / λ̄[j])
                end
            end
            F[D + k] = λ[k] - income_sum * λ[D]

        elseif f[k] == 2
            # Ricardian consumer lambda (row 5N+3)
            # λ_ric = -income_sum · λ_last + Σ_{j: f[j]=0} λ_j
            income_sum = 0.0
            factor_lambda_sum = 0.0
            for j in 1:D
                if f[j] == 0
                    income_sum += λ̄[j] * (1 - φ_htm[j]) * (1 - λ[j] / p[j] / λ̄[j])
                    factor_lambda_sum += λ[j]
                end
            end
            F[D + k] = λ[k] - (-income_sum * λ[D] + factor_lambda_sum)
        end
    end

    # ──────────────────────────────────────────────
    # 3. Numeraire: p[D] = 1
    # ──────────────────────────────────────────────
    # F[D] was set by a price equation or factor clearing; override to numeraire
    F[D] = p[D] - 1.0

    return nothing
end

"""
    solve_equilibrium(m::MCPModel; z0=nothing, tol=1e-10, maxiter=1000) -> (p, λ, converged, iterations)

Solve the equilibrium system using NLsolve.

Returns (p vector, λ vector, converged::Bool, iterations::Int).
Throws on fatal solver errors.
"""
function solve_equilibrium(m::MCPModel; z0=nothing, tol=1e-10, maxiter=1000)
    D = m.D
    if z0 === nothing
        z0 = vcat(m.init_p, m.init_lambda)
    end

    function residual_wrapper!(F, z)
        equilibrium_residual!(F, z, m)
    end

    result = nlsolve(residual_wrapper!, z0;
                     method=:trust_region,
                     ftol=tol,
                     iterations=maxiter,
                     store_trace=false,
                     show_trace=false)

    z_sol = result.zero
    p = z_sol[1:D]
    λ = z_sol[D+1:2D]
    converged = result.f_converged || result.x_converged
    return p, λ, converged, result.iterations
end

"""
    eq_continuation(m_base::MCPModel, shock_A, shock_B, shock_B_agg, t_grid;
                    show_trace=false) -> Dict

Solve the equilibrium for each t in t_grid, using the previous solution
as initializer for the next (continuation method). Returns a Dict of arrays.
"""
function eq_continuation(m_base::MCPModel,
                          shock_A, shock_B,
                          t_grid::Vector{Float64};
                          show_trace=false)

    D = m_base.D
    N = m_base.N
    n_t = length(t_grid)

    # Storage
    GDP_stack = zeros(n_t)
    nom_GDP_stack = zeros(n_t)
    Hulten_stack = zeros(n_t)
    prices_stack = zeros(D, n_t)
    lambda_stack = zeros(D, n_t)
    retcodes = zeros(Int, n_t)

    # Precompute Psi_re (invariant across t)
    Psi_re = inv(I - m_base.Omega_re)

    z0 = vcat(m_base.init_p, m_base.init_lambda)

    for (ti, t) in enumerate(t_grid)
        # Build shock vectors for this t
        # (Same logic as Master_file_3.m shock_type=1: supply + demand)
        A_shock = copy(m_base.A)
        B_shock = copy(m_base.B)

        # Labor supply shock: A[3N+2:4N+1] = (1 - t * BLS_shock)
        A_shock[(3*N+2):(4*N+1)] .= 1.0 .- t .* shock_A

        # Sectoral demand: B[1, 2:N+1] = (1-t*0.66) + t*0.66*(1+PCE_shock)
        B_shock[1, 2:(N+1)] .= (1 - t * 0.66) .+ t * 0.66 .* (1.0 .+ shock_B)
        # Renormalize
        B_shock[1, :] ./= sum(m_base.Omega_re[1, :] .* B_shock[1, :]')

        # Aggregate demand shock: B[5N+3, D] = 1 + min(ti-1,1) * 0.105
        B_shock[5*N+3, D] = 1.0 + min(ti - 1, 1) * 0.105
        B_shock[5*N+3, :] ./= sum(m_base.Omega_re[5*N+3, :] .* B_shock[5*N+3, :]')

        # Build the model for this t
        m_t = MCPModel(
            m_base.D, m_base.N,
            m_base.Omega_re,
            m_base.factor, m_base.keynes,
            m_base.theta, m_base.cobb_douglas,
            m_base.phi, m_base.phi_htm,
            m_base.mu, m_base.in_mu,
            A_shock, B_shock,
            m_base.init_lambda, m_base.init_p,
            m_base.chi
        )

        p, λ, converged, iters = solve_equilibrium(m_t; z0=z0)

        # Store
        prices_stack[:, ti] .= p
        lambda_stack[:, ti] .= λ
        GDP_stack[ti] = λ[1] / p[1]
        nom_GDP_stack[ti] = λ[1]
        retcodes[ti] = converged ? 0 : 1

        # Hulten benchmark: exp(Ψ_re[1,:] · log(A))
        Hulten_stack[ti] = exp(sum(Psi_re[1, :] .* log.(A_shock)))

        # Use current solution as next initializer
        z0 = vcat(p, λ)

        if show_trace
            @printf("t=%.3f  GDP=%.6f  nomGDP=%.6f  conv=%d\n",
                    t, GDP_stack[ti], nom_GDP_stack[ti], converged)
        end
    end

    return Dict(
        "t_grid" => t_grid,
        "GDP" => GDP_stack,
        "nominal_GDP" => nom_GDP_stack,
        "Hulten" => Hulten_stack,
        "prices" => prices_stack,
        "lambda" => lambda_stack,
        "retcodes" => retcodes
    )
end