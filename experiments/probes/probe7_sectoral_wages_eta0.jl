# experiments/probes/probe7_sectoral_wages_eta0.jl
#
# Option C of ADR-0020: sector-specific wages at the eta = 0 (BF) endpoint.
#
# Spec under test. At eta = 0 labour cannot reallocate, so the sectoral
# allocation stays frozen at `data.labor_share` and each sector's wage is set
# by that sector's own marginal product at the frozen allocation. Unknowns and
# equations (N = 71):
#
#   X = [p(1:N); y(1:N); w(1:N); F]          3N + 1 unknowns
#   1. zero-profit, per sector:  p_i = cost_i(p, w_i)            (N)
#   2. sectoral FOC at the frozen L_i:
#        log L^cm_i(p_i, y_i, w_i) = log labor_share_i           (N)
#   3. all-N clearing, household income = sum_i w_i L_i          (N)
#   4. numeraire: CPI = 1                                        (1)
#
# Why F is kept: the block {zero-profit, FOC, clearing} is homogeneous of
# degree 1 in (p, w, F), so 3N equations determine 3N - 1 effective unknowns;
# with the FOC replacing the single aggregate labour equation of the mobile
# system the block has one equation more than it has directions to pin, and the
# demand block needs one free scalar to be consistent with the supply side. F
# is that scalar (it enters E after tax, exactly as in ADR-0019). The external
# account then closes by construction: at a FOC solution the identity gap is
# `-(sum_i w_i L^cm_i - sum_i w_i L_i)` = 0, because the frozen allocation is
# cost-minimizing at these wages.
#
# The probe solves the three BF cells (F1/F2/F3), checks the residual, the
# clearing block, the FOC gap, the account closure, F's identification (FD
# Jacobian column + condition number) and multi-start invariance, and reports
# the sectoral wage dispersion plus a comparison with the pinned-F = 0 row of
# the v5 generation.

include(joinpath(@__DIR__, "..", "run.jl"))

using Printf, LinearAlgebra
import NonlinearSolve

const DESIGN = "matrix_5x3_v5"

"Cost-minimizing labour demand with a SECTORAL wage vector (mirrors
`_cost_minimizing_labor`, which hard-codes `log(w)` for a scalar wage)."
function lcm_vec(p, y, w, model)
    (; data, options, shocks) = model
    ϵ = options.elasticities.ϵ
    p, y, w, A, α = BeyondHulten._positive_floor.((p, y, w, shocks.supply_shock, data.factor_share))
    log_demand = ϵ .* (log.(p) .+ ((ϵ - 1) / ϵ) .* log.(A) .+ (1 / ϵ) .* log.(α) .+
        (1 / ϵ) .* log.(y) .- log.(w))
    exp.(clamp.(log_demand, log(floatmin(Float64)), log(floatmax(Float64))))
end

"""
Demand block with a SECTORAL wage vector and the frozen eta = 0 allocation.
Mirrors `_mobile_market_demand`, with two substitutions: `L_i` is the frozen
`data.labor_share` (not `sectoral_labor_demand`), and household wage income is
`sum_i w_i L_i` (each sector's workers earn that sector's wage). No legacy
manna (the matrix designs pass zero manna; asserted below).
"""
function bf_sw_blocks(model, p::AbstractVector, y::AbstractVector,
        w::AbstractVector, F::Real)
    (; data, options, shocks) = model
    N = length(data.factor_share)
    (; consumption_share, Ω_raw, factor_share, labor_share) = data
    (; θ, ϵ, σ) = options.elasticities
    ip = BeyondHulten._intermediate_price(Ω_raw, p, θ)
    fin = model.financing
    ds_eff = BeyondHulten.preference_weights(fin, shocks.demand_shock)
    L = labor_share
    income = dot(w, L)                                  # sectoral wage bill
    program = fin isa Union{TaxFinanced, ExternalDebt} ? fin.g : zeros(N)
    tau = dot(p, data.gov_demand .+ (fin isa TaxFinanced ? program : zeros(N))) / income
    E = (1 - tau) * income + F
    agg = sum(consumption_share .* ds_eff .* p .^ (1 - σ))
    c_dom = (1 .- data.saving_rate) .* (1 .- data.import_margin) .*
            (consumption_share .* ds_eff) .* E .* p .^ (-σ) ./ agg
    inter = p .^ (-θ) .* (Ω_raw' * (p .^ ϵ .* shocks.supply_shock .^ (ϵ - 1) .*
        ip .^ (θ - ϵ) .* (data.A_bill ./ data.λ) .* y))
    final = c_dom .+ (1 .- data.import_margin) .* program .+
            (1 .- data.import_margin) .* (data.gov_demand .+ data.exo_demand) .+
            data.exports_demand
    cost = BeyondHulten._ces_unit_cost(shocks.supply_shock, factor_share, w, ip, ϵ)
    Lcm = lcm_vec(p, y, w, model)
    (; inter, final, c_dom, cost, E, income, tau, Lcm, L, program,
        cpi = sum(consumption_share .* p .^ (1 - σ))^(1 / (1 - σ)))
end

"Residual of the sectoral-wage eta = 0 system (3N+1 equations)."
function problem_bf_sw(out::Vector, X::Vector, model)
    N = length(model.data.factor_share)
    p = BeyondHulten._positive_floor(X[1:N])
    y = BeyondHulten._positive_floor(X[N+1:2N])
    w = BeyondHulten._positive_floor(X[2N+1:3N])
    F = X[3N+1]
    b = bf_sw_blocks(model, p, y, w, F)
    out[1:N] .= p .- b.cost
    out[N+1:2N] .= log.(b.Lcm) .- log.(b.L)
    out[2N+1:3N] .= y .- b.inter .- b.final
    out[3N+1] = b.cpi - 1.0
    nothing
end

"Full external-account decomposition at the sectoral-wage solution."
function sw_account(model, X)
    N = length(model.data.factor_share)
    p = X[1:N]; y = X[N+1:2N]; w = X[2N+1:3N]; F = X[3N+1]
    b = bf_sw_blocks(model, p, y, w, F)
    (; data, options, shocks) = model
    (; θ, ϵ) = options.elasticities
    m = data.import_margin
    # same valuation as external_balance_canary, with the sectoral-wage blocks
    M_cons = dot(p .* (m ./ max.(1 .- m, eps(Float64))), b.c_dom)
    M_inj = dot(p .* m, b.program .+ data.gov_demand .+ data.exo_demand)
    k_bill = p .^ ϵ .* shocks.supply_shock .^ (ϵ - 1) .*
             BeyondHulten._intermediate_price(data.Ω_raw, p, θ) .^ (1 - ϵ)
    M_intl = dot(k_bill .* (data.M_int ./ data.λ), y)
    T_intl = dot(k_bill .* (data.T_int ./ data.λ), y)
    S = data.saving_rate * b.E
    IX = dot(p, data.exo_demand .+ data.exports_demand)
    B_gov = model.financing isa ExternalDebt ? dot(p, b.program) : 0.0
    M = M_cons + M_inj + M_intl
    resource = S + T_intl + M - IX
    (; S, T = T_intl, M, IX, F, B_gov, resource, booked = F + B_gov,
        gap = resource - (F + B_gov),
        labour_wedge = dot(w, b.Lcm) - dot(w, b.L),   # predicted -gap
        income = b.income, E = b.E, Lcm = b.Lcm, L = b.L, w = w, cpi = b.cpi,
        resid = maximum(abs, begin
            o = similar(X); problem_bf_sw(o, X, model); o
        end))
end

"Solve with the ADR-0015 ladder (Newton + residual-gated LM polish)."
function solve_sw(model, init)
    f = problem_bf_sw
    x = copy(init)
    res = NonlinearSolve.solve(NonlinearSolve.NonlinearProblem(f, x, model),
        reltol = 1e-8, abstol = 1e-8, maxiters = 40000)
    x = res.u
    rmax = maximum(abs, begin
        o = similar(x); f(o, x, model); o
    end)
    for _ in 1:6
        rmax <= 1e-11 && break
        res = NonlinearSolve.solve(NonlinearSolve.NonlinearProblem(f, x, model),
            NonlinearSolve.LevenbergMarquardt(); reltol = 1e-13, abstol = 1e-13, maxiters = 40000)
        x_new = res.u
        r_new = maximum(abs, begin
            o = similar(x_new); f(o, x_new, model); o
        end)
        r_new < rmax || break
        x, rmax = x_new, r_new
    end
    x, rmax
end

function reconstruct_main()
    root = default_root()
    design_d = load_design(DESIGN; root = root)
    ref = build_reference(design_d; root = root)
    data = ref.data
    N = length(data.factor_share)
    ψ, g = programme_vectors(design_d, data, collect(1:N); root = root)
    @printf("N = %d, unknowns = 3N+1 = %d, equations = %d\n", N, 3N + 1, 3N + 1)

    for fin_id in ("F1", "F2", "F3")
        cell = design_d["cells"]["matrix_5x3-v5-BF-$fin_id"]
        model = build_cell_model(cell, design_d, data, ψ, g)
        @assert all(iszero, model.shocks.autonomous_demand) &&
                all(iszero, model.shocks.investment_shock) "probe assumes zero legacy manna"
        # warm start: the pinned-F = 0 BF solution, wage vector at one
        bf_pinned = build_cell_model(cell, design_d, data, ψ, g)
        sol0 = solve(bf_pinned; init = ref.init_warm)
        init = [sol0.prices_raw; sol0.quantities; ones(N); 0.0]
        x, rmax = solve_sw(model, init)
        a = sw_account(model, x)

        println("\n=== BF-$fin_id (sectoral wages, eta = 0) ===")
        @printf("  max|resid| = %.3e   (equations 3N+1)\n", rmax)
        @printf("  F = %+.10f   B_gov = %+.10f   booked = %+.10f\n", a.F, a.B_gov, a.booked)
        @printf("  resource = S + T + M - (I+X) = %+.10f   gap = %+.3e\n", a.resource, a.gap)
        @printf("  labour wedge  sum w L^cm - sum w L = %+.3e  (gap = -wedge)\n", a.labour_wedge)
        @printf("  FOC gap max|log L^cm - log L| = %.3e\n",
            maximum(abs, log.(a.Lcm) .- log.(a.L)))
        @printf("  income sum w_i L_i = %.10f   E = %.10f   CPI = %.12f\n", a.income, a.E, a.cpi)
        @printf("  wages: min = %.6f  max = %.6f  ratio = %.1f  sum w_i = %.4f\n",
            minimum(a.w), maximum(a.w), maximum(a.w) / minimum(a.w), sum(a.w))
        top = sortperm(a.w; rev = true)[1:3]
        bot = sortperm(a.w)[1:3]
        for i in top
            @printf("    high wage sector %2d: w = %10.6f  L = %.6f  y = %.6f\n", i, a.w[i], a.L[i], x[N+i])
        end
        for i in bot
            @printf("    low  wage sector %2d: w = %10.6f  L = %.6f  y = %.6f\n", i, a.w[i], a.L[i], x[N+i])
        end

        # identification of F: finite-difference Jacobian at the solution
        J = zeros(3N + 1, 3N + 1)
        h = 1e-7
        for j in 1:(3N+1)
            xp = copy(x); xp[j] += h
            xm = copy(x); xm[j] -= h
            op = similar(x); om = similar(x)
            problem_bf_sw(op, xp, model); problem_bf_sw(om, xm, model)
            J[:, j] .= (op .- om) ./ (2h)
        end
        sv = svdvals(J)
        @printf("  Jacobian: cond = %.3e   F-column norm = %.3e   sigma_min/sigma_max = %.3e\n",
            cond(J), norm(J[:, end]), minimum(sv) / maximum(sv))

        # multi-start invariance
        x2, r2 = solve_sw(model, [ones(N); data.λ; fill(1.2, N); 0.01])
        @printf("  multi-start: max|dX| = %.3e (resid %.1e), dF = %.3e\n",
            maximum(abs, x2 .- x), r2, abs(x2[end] - x[end]))
    end
    return nothing
end

reconstruct_main()
