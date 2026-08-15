
# elasticity_gradient.jl -- Elasticity-gradient sweep (R4), ports elasticity_gradient.m
# and eg.m.
#
# Reproduces the "which elasticity matters" figure (B&F 2019, Fig. 7):
#   * for each sector in `sectors`, apply an idiosyncratic TFP shock (A[sec] = shock_val)
#   * sweep ONE elasticity over range(0.015, 0.99) while holding the other two
#     at their baseline (0.015), recording:
#        - the REAL (CES-welfare) GDP, exactly as eg.m:
#              GDP = (shock.^((eps-1)/eps) .* alpha.^(1/eps) .* y.^(1/eps)
#                     .* (1 ./ L).^(1/eps))' * L
#          i.e. real_gdp_mc(sol, shock, alpha, L, eps) -- note the swept eps is
#          used for the epsilon sweep, baseline eps otherwise.
#        - the Domar-weighted mean price, eg.m:
#              mean_prices = sum(beta .* p.^(-sigma))   (denominator = sum(beta) = 1)
#          using the swept sigma for the sigma sweep, baseline sigma otherwise.
#   * output: 3 x 100 matrices for GDP and mean price, per swept elasticity.
#
# NOTE: the MATLAB uses baseline elasticities eps_orig = theta_orig = sigma_orig
# = 0.015 (NOT the 0.5/0.001/0.9 production calibration) for this figure.
#
# To make the Newton solver track the solution across the (near-degenerate)
# extremes of the sweep, we use CONTINUATION: each grid point is warm-started
# from the previous one's solution. This mirrors how fmincon would track the
# root and avoids the diverging-step failures seen when starting cold.

function _eg_gdp(sol::BFSolution, shock::Vector{Float64}, α::Vector{Float64},
                 L::Vector{Float64}, ε::Float64)
    return real_gdp_mc(sol, shock, α, L, ε)
end

function _eg_mp(sol::BFSolution, β::Vector{Float64}, σ::Float64)
    return sum(β .* (sol.p .^ (-σ)))
end

function run_elasticity_gradient(data;
        sectors::Vector{Int} = [10, 20, 23],
        shock_val::Float64 = 1.1,
        grid = range(0.015, 0.99; length = 100),
        eps_orig::Float64 = 0.015,
        theta_orig::Float64 = 0.015,
        sigma_orig::Float64 = 0.015,
        jacobian::Symbol = :numerical)
    N   = data.N
    g   = collect(grid)
    ng  = length(g)
    ns  = length(sectors)

    gdp_sigma = zeros(ns, ng); gdp_theta = zeros(ns, ng); gdp_eps = zeros(ns, ng)
    mp_sigma  = zeros(ns, ng); mp_theta  = zeros(ns, ng); mp_eps  = zeros(ns, ng)

    for (i, sec) in enumerate(sectors)
        shock = ones(N); shock[sec] = shock_val

        # sweep sigma (eps, theta at baseline) -> GDP uses baseline eps, mp uses swept sigma
        prev_X = nothing
        for (j, sigma) in enumerate(g)
            sol = solve_bf(shock, data.Ω, data.α, data.β, data.L, eps_orig, theta_orig, sigma;
                           jacobian = jacobian, warm_start = prev_X)
            if sol.converged
                gdp_sigma[i, j] = _eg_gdp(sol, shock, data.α, data.L, eps_orig)
                mp_sigma[i, j]  = _eg_mp(sol, data.β, sigma)
                prev_X = vcat(sol.p, sol.y)
            else
                gdp_sigma[i, j] = 1.0; mp_sigma[i, j] = 1.0
            end
        end
        # sweep theta (eps, sigma at baseline) -> GDP uses baseline eps, mp uses baseline sigma
        prev_X = nothing
        for (j, theta) in enumerate(g)
            sol = solve_bf(shock, data.Ω, data.α, data.β, data.L, eps_orig, theta, sigma_orig;
                           jacobian = jacobian, warm_start = prev_X)
            if sol.converged
                gdp_theta[i, j] = _eg_gdp(sol, shock, data.α, data.L, eps_orig)
                mp_theta[i, j]  = _eg_mp(sol, data.β, sigma_orig)
                prev_X = vcat(sol.p, sol.y)
            else
                gdp_theta[i, j] = 1.0; mp_theta[i, j] = 1.0
            end
        end
        # sweep epsilon (theta, sigma at baseline) -> GDP uses swept eps, mp uses baseline sigma
        prev_X = nothing
        for (j, eps) in enumerate(g)
            sol = solve_bf(shock, data.Ω, data.α, data.β, data.L, eps, theta_orig, sigma_orig;
                           jacobian = jacobian, warm_start = prev_X)
            if sol.converged
                gdp_eps[i, j] = _eg_gdp(sol, shock, data.α, data.L, eps)
                mp_eps[i, j]  = _eg_mp(sol, data.β, sigma_orig)
                prev_X = vcat(sol.p, sol.y)
            else
                gdp_eps[i, j] = 1.0; mp_eps[i, j] = 1.0
            end
        end
    end

    return (gdp_eps = gdp_eps, gdp_theta = gdp_theta, gdp_sigma = gdp_sigma,
            mp_eps = mp_eps, mp_theta = mp_theta, mp_sigma = mp_sigma,
            grid = g, sectors = sectors,
            eps_orig = eps_orig, theta_orig = theta_orig, sigma_orig = sigma_orig)
end
