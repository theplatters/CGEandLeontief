# ═══════════════════════════════════════════════════════════════════════════════
# Mobile Labor with Geometric Intersectoral Reallocation η
# ═══════════════════════════════════════════════════════════════════════════════
#
# This file implements the key model extension: replacing sector-specific
# (immobile) labor with geometric intersectoral reallocation parameter η.
#
#   η = 0     → immobile baseline allocation
#   η = 1     → fully cost-minimizing allocation
#   η outside [0, 1] → extrapolation of the same geometric rule
#
# The key change vs. the base CES model:
#   - The unknown vector gains a scalar wage w: X = [p(1:N); y(1:N); w]
#   - Labor allocation is endogenous: L_i = (∂Y_i/∂L_i = w) → solved from FOC
#   - Labor market clearing: total employment equals the fixed labor bar
#
# Author: calculato (AI research assistant)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    MobileLaborCESElasticities

Elasticity parameters for the mobile-labor CES model.

- `θ` : elasticity of substitution between intermediate goods
- `ϵ`  : elasticity of substitution between labor and intermediate composite
- `σ`  : elasticity of substitution in consumption
- `η`  : intersectoral reallocation parameter; it is not a labor-supply elasticity
"""
struct MobileLaborCESElasticities <: AbstractElasticities
    θ::Float64
    ϵ::Float64
    σ::Float64
    η::Float64
    # BETA closure (Stage 1.2): elasticity of TOTAL labour supply along the
    # real wage; 0.0 = vertical supply (ALPHA). Used only by the :beta closure.
    eta_s::Float64
end
# Parent-compatible 4-arg constructor (no elastic labour supply).
MobileLaborCESElasticities(θ::Real, ϵ::Real, σ::Real, η::Real) =
    MobileLaborCESElasticities(Float64(θ), Float64(ϵ), Float64(σ), Float64(η), 0.0)

"""
    MobileLaborCES

Model type for CES with geometric intersectoral labor reallocation and an explicit wage-regime closure.

- `elasticities` : a `MobileLaborCESElasticities` struct
- `labor_bar`    : fixed total employment capacity L̄ (defaults to Σ labor_share)
"""
struct MobileLaborCES <: ModelType
    elasticities::MobileLaborCESElasticities
    labor_bar::Float64
    closure::Symbol   # :mobile = flexible wage (ALPHA); :fixed = sticky wage (GAMMA);
                      # :beta = flexible wage with elastic total labour supply (BETA)
    function MobileLaborCES(elasticities::MobileLaborCESElasticities, labor_bar::Float64, closure::Symbol)
        closure in (:mobile, :fixed, :beta) || throw(ArgumentError("unsupported MobileLaborCES closure $closure; use :mobile, :fixed, or :beta"))
        new(elasticities, labor_bar, closure)
    end
end

_closure_symbol(closure::Symbol) = closure
_closure_symbol(::FlexibleWageClosure) = :mobile
_closure_symbol(::FixedWageClosure) = :fixed
_closure_symbol(closure) = throw(ArgumentError(
    "unsupported MobileLaborCES closure $closure; use :mobile, :fixed, FlexibleWageClosure(), or FixedWageClosure()"))

MobileLaborCES(e::MobileLaborCESElasticities, lb::Real, closure::AbstractLaborClosure) =
    MobileLaborCES(e, Float64(lb), _closure_symbol(closure))
MobileLaborCES(e::MobileLaborCESElasticities, lb::Real; closure=:mobile) =
    MobileLaborCES(e, Float64(lb), _closure_symbol(closure))
MobileLaborCES(e::MobileLaborCESElasticities, lb::Real, closure::Symbol) = MobileLaborCES(e, Float64(lb), closure)

function MobileLaborCES(elasticities::MobileLaborCESElasticities, data::Data)
    labor_bar = sum(data.labor_share)
    MobileLaborCES(elasticities, labor_bar, :mobile)
end

labor_closure(options::MobileLaborCES) =
    options.closure == :mobile ? FlexibleWageClosure() :
    options.closure == :beta   ? ElasticLaborClosure(options.elasticities.eta_s) :
                                 FixedWageClosure()

# ── Labor-market equation hook ──
# The total-labor-market residual of the flexible-wage systems (problem(), 2N+1).
# ALPHA (FlexibleWageClosure): vertical supply at the bar L̄.
# BETA  (ElasticLaborClosure, defined in cbase2/src/closures.jl): elastic supply
#       L^s = L̄ · (w / w0)^{η_s} along the real wage (CPI = 1 ⇒ w is the real wage).
labor_market_residual(::FlexibleWageClosure, model::Model{MobileLaborCES}, L_sum::Real, w::Real) =
    L_sum - model.options.labor_bar

# ── Unit cost with the Cobb-Douglas limit guard (Stage 1.5, Milestone C) ──
# CES unit cost:  (A^(ϵ-1) · (fs·w^(1-ϵ) + (1-fs)·ip^(1-ϵ)))^(1/(1-ϵ)).
# At ϵ = 1 exactly the formula degenerates (x^Inf with x = 1 ± float error);
# the analytic Cobb-Douglas limit is  w^fs · ip^(1-fs) / A  and is used whenever
# |1 - ϵ| is below machine-relevant tolerance, keeping the CD special case
# sign-safe and continuous.
const _CD_LIMIT_TOL = 1e-9

function _ces_unit_cost(A_eff, fs, w, ip, ϵ)
    if abs(1 - ϵ) < _CD_LIMIT_TOL
        return (w .^ fs) .* (ip .^ (1 .- fs)) ./ A_eff
    end
    ((A_eff .^ (ϵ - 1)) .* (fs .* w .^ (1 - ϵ) .+ (1 .- fs) .* ip .^ (1 - ϵ))) .^ (1 / (1 - ϵ))
end

# ─────────────────────────────────────────────────────────────────────────────────
# Sectoral labor demand from the wage (marginal product condition)
# ─────────────────────────────────────────────────────────────────────────────────
# In the base CES model, the wage equation (from FOC of labor) is:
#   w_i = p_i · A_i^((ε-1)/ε) · α_i^(1/ε) · y_i^(1/ε) · L_i^(-1/ε)
#
# With mobile labor, w is economy-wide (scalar). Inverting for L_i:
#   L_i = [ p_i · A_i^((ε-1)/ε) · α_i^(1/ε) · y_i^(1/ε) / w ]^ε
#
# This gives sectoral labor demand as a function of (p, y, w).
# ─────────────────────────────────────────────────────────────────────────────────

const _ETA_MAX = 50.0
const _ETA_SCALE_INDETERMINACY_TOL = 1.1e-6
const _ALLOCATION_LOG_LIMIT = 50.0
const _WEDGE_MIN = eps(Float64)

function _checked_eta(η)
    isfinite(η) || throw(ArgumentError("η must be finite"))
    abs(η) <= _ETA_MAX || throw(DomainError(η, "η must have magnitude at most 50"))
    Float64(η)
end

_positive_floor(x) = max.(x, eps(Float64))
_positive_floor(x::Real) = max(x, eps(Float64))

"""
	_intermediate_price(Ω_raw, p, θ)

Intermediate-goods price index. CES for θ ≠ 1; the exact Cobb-Douglas
limit ip_u = prod_j p_j^{Ω[u,j]} at θ = 1 (the bounded production-core
choice of Notebook 03b: with genuine sector self-loops -- up to Omega_ii
= 0.57 in the data -- the CES index at theta < 1 is self-referencing and
the zero-profit price system has exploding spiral branches; the
Cobb-Douglas index is log-linear with a unique positive root).

Both branches are valid unit-expenditure functions only when every row of
Ω_raw is a probability vector (Σ_j Ω_raw[u,j] = 1). `conditional_input_shares`
enforces this for the source table and `retained_dataset` re-asserts it after
any sector rebuild; without the renormalization the CES branch silently
rescales with a deficit-dependent factor (the 2026-09-17 θ-discontinuity).
"""
function _intermediate_price(Ω_raw::AbstractMatrix, p::AbstractVector, θ::Real)
	isapprox(θ, 1.0; rtol = 0, atol = 1e-6) && return exp.(Ω_raw * log.(p))
	return (Ω_raw * p .^ (1 - θ)) .^ (1 / (1 - θ))
end

"""Cost-minimizing labor demand, evaluated in log space for stable η extrapolation."""
function _cost_minimizing_labor(p, y, w, model::Model{MobileLaborCES})
    (; data, options, shocks) = model
    (; ϵ) = options.elasticities
    (; factor_share) = data
    p, y, w, A, α = _positive_floor.((p, y, w, shocks.supply_shock, factor_share))
    log_demand = ϵ .* (log.(p) .+ ((ϵ - 1) / ϵ) .* log.(A) .+ (1 / ϵ) .* log.(α) .+
        (1 / ϵ) .* log.(y) .- log(w))
    exp.(clamp.(log_demand, log(floatmin(Float64)), log(floatmax(Float64))))
end

"""Interpolate fixed and cost-minimizing allocations geometrically in log space."""
function _interpolated_labor(p, y, w, model::Model{MobileLaborCES})
    η = _checked_eta(model.options.elasticities.η)
    fixed = _positive_floor(model.data.labor_share)
    optimum = _positive_floor(_cost_minimizing_labor(p, y, w, model))
    log_labor = (1 - η) .* log.(fixed) .+ η .* log.(optimum)
    exp.(clamp.(log_labor, log(floatmin(Float64)), log(floatmax(Float64))))
end

"""Return the allocative-efficiency factor (one at the optimum, never above one).

The log allocation ratio is limited to ±50 to keep trial solver iterates finite.
"""
function _allocation_efficiency_wedge(labor, optimum, factor_share, ϵ)
    ϵ > 0 || throw(DomainError(ϵ, "ϵ must be positive for the allocation wedge"))
    ratio = clamp.(log.(_positive_floor(labor) ./ _positive_floor(optimum)),
                   -_ALLOCATION_LOG_LIMIT, _ALLOCATION_LOG_LIMIT)
    curvature = abs(1 - ϵ) / ϵ
    penalty = 0.5 .* _positive_floor(factor_share) .* max.(1 .- factor_share, 0) .* curvature .* ratio.^2
    # Keep the wedge strictly positive even when the clamped log ratio makes
    # the exponential underflow. `max.` remains elementwise AD-compatible.
    max.(exp.(-penalty), _WEDGE_MIN)
end

function _allocation_efficiency_wedge(p, y, w, labor, model::Model{MobileLaborCES})
    (; factor_share) = model.data
    ϵ = model.options.elasticities.ϵ
    optimum = _cost_minimizing_labor(p, y, w, model)
    _allocation_efficiency_wedge(labor, optimum, factor_share, ϵ)
end

"""
    sectoral_labor_demand(p, y, w, model::Model{MobileLaborCES})

Compute sectoral labor demand L_i given prices, quantities, and the economy-wide wage.
Derived by inverting the marginal-product-of-labor condition.
"""
function sectoral_labor_demand(p, y, w, model::Model{MobileLaborCES})
    _interpolated_labor(p, y, w, model)
end

"""
    economy_wide_wage(p, y, labor, model::Model{MobileLaborCES})

Compute the sector-implied wage from the base CES FOC (used for initialization).
When labor is mobile, all sectors should imply the same wage w.
"""
function economy_wide_wage(p, y, labor, model::Model{MobileLaborCES})
    (; data, options, shocks) = model
    (; ϵ) = options.elasticities
    (; supply_shock) = shocks
    (; factor_share) = data

    # w_i = p_i · A_i^((ε-1)/ε) · α_i^(1/ε) · y_i^(1/ε) · L_i^(-1/ε)
    w_vec = p .* (supply_shock .^ ((ϵ - 1) / ϵ)) .* (factor_share .^ (1 / ϵ)) .* (y .^ (1 / ϵ)) .* (labor .^ (-1 / ϵ))
    # Return the GDP-weighted average as a single wage
    return sum(w_vec .* labor) / sum(labor)
end

# ═══════════════════════════════════════════════════════════════════════════════
# The equilibrium problem (2N+1 equations, 2N+1 unknowns)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    problem(out, X, model::Model{MobileLaborCES})

The equilibrium system for the mobile-labor CES model.

Unknowns: X = [p(1:N); y(1:N); w]  — 2N+1 elements
Equations (2N+1):
  1. Zero-profit equations  (N):     p_i = cost_i(p, w)   for all i=1..N
  2. Market clearing        (N-1):   y_i = intermediary_demand_i + final_demand_i
                                      for i=1..N-1
  3. Labor market clearing  (1):     Σ L_i(p,y,w) = L̄
  4. Numeraire              (1):     CPI = 1  (Σ β_i · p_i^(1-σ))^(1/(1-σ)) = 1)

Note (2026-09-17, user decision: back to the sanctioned form). With the v3
homogeneous budget the system is scale-homogeneous in (p, w): enforcing ALL N
clearing equations leaves solutions on RAYS (measured: r_zp(μp, y, μw) =
μ·r_zp, r_rest invariant, to 9.3e-15), so the all-N/no-numeraire variant is
ray-degenerate with a singular Jacobian (deflated σ_min analysis,
process_comments.md). The numeraire CPI = 1 is the GAUGE FIX that selects the
ray member; the enforced system is square with a generically nonsingular
Jacobian. The N-th market is not independent: the ray member with CPI = 1
clears all N markets, which is asserted post-solve as the canary
(verify_v3.jl). The price level is pinned by the numeraire, not by the real
injections (which are scale-invariant).

Economic note: η changes the geometric allocation between baseline and
cost-minimizing sectoral labor demand. It is not a labor-supply elasticity.
"""
function problem(out::Vector, X::Vector, model::Model{MobileLaborCES})
    (; data, options, shocks) = model
    N = length(data.factor_share)

    # FORMULATION (2026-09-17, user decision: back to the sanctioned
    # N−1+CPI form). 2N+1 unknowns: p1..pN, y1..yN, w. The numeraire CPI = 1
    # is the gauge fix of the scale-homogeneous system (rays in (p, w); see
    # the docstring above); the N-th clearing equation is implied at the
    # selected ray member and asserted post-solve as the canary.
    p = _positive_floor(X[1:N])
    y = _positive_floor(X[N+1:2N])
    w = max(X[2N+1], 1e-10)  # scalar wage, keep positive

    (; supply_shock, demand_shock) = shocks
    (; consumption_share, Ω_raw, factor_share, labor_share) = data
    (; θ, ϵ, σ, η) = options.elasticities

    # ── Intermediate goods price index ──
    intermediate_price = _intermediate_price(Ω_raw, p, θ)

    # ── CPI (consumption price index) ──
    cpi = sum(consumption_share .* p .^ (1 - σ))^(1 / (1 - σ))

    # ── Sectoral labor demand (from inverted FOC) ──
    L_i = sectoral_labor_demand(p, y, w, model)

    # ── Final demand ──
    # Budget-consistent CES demand over the household's (possibly financed)
    # expenditure. `fin = model.financing` (Foundation II):
    #   F1  composes E within its budget via preference weights (no additive
    #       demand; the normalizer keeps Σ p_i c_i = E exactly);
    #   F2  E = w·ΣL − T(p) with T = Σ p_i g_i the lump-sum tax;
    #   F3  E = w·ΣL, externally financed.
    # The unfinanced autonomous/investment manna of the legacy pipeline is
    # retired; additive demand enters only through `additive_demand(fin, N)`.
    fin = model.financing
    ds_eff = preference_weights(fin, demand_shock)
    L_sum = sum(L_i)
    total_income = w * L_sum
    # NOTE: no positivity guard here — the residual function must tolerate the
    # solver's exploration of negative-income trial points (the legacy code
    # did). The E > 0 check belongs to the post-solve validation in the
    # notebooks (headline assertions verify E = w*L - T > 0 at equilibrium).
    E = household_expenditure(fin, model, total_income, p, L_sum)
    agg = sum(consumption_share .* ds_eff .* p .^ (1 - σ))
    # v3 open economy: the household consumes (1-s)E gross (saving s·E leaks);
    # only the DOMESTIC content circulates -- the import content of household,
    # government, investment and programme demand is supplied by the external
    # account. Exports are exogenous domestic sales (no margin). The leakages
    # s·E + M now balance the injections I + X -- the identity S = I + X - M
    # holds at every equilibrium (validated in the notebooks).
    c_dom = (1 .- data.saving_rate) .* (1 .- data.import_margin) .*
            (consumption_share .* ds_eff) .* E .* p .^ (-σ) ./ agg
    total_final_demand = c_dom .+
                         (1 .- data.import_margin) .* additive_demand(fin, N) .+
                         (1 .- data.import_margin) .* (data.gov_demand .+ data.exo_demand) .+
                         data.exports_demand

    # ── Intermediary demand ──
    intermediary_demand = p .^ (-θ) .* (Ω_raw' * (p .^ ϵ .* supply_shock .^ (ϵ - 1) .* intermediate_price .^ (θ - ϵ) .* (1 .- factor_share) .* y))

    # ── B&F (2019) labor-reallocation wedge ──
    # L_opt is the cost-minimizing (mobile) labor demand; L_base is the
    # baseline (immobile, fixed) allocation. When labor cannot reallocate toward
    # its shocked optimum (eta < 1), there is a second-order allocative
    # inefficiency that raises effective cost / lowers TFP. The wedge vanishes
    # for eta=1 (mobile) and at baseline (L_opt == L_base).
    alloc_wedge = _allocation_efficiency_wedge(p, y, w, L_i, model)

    # ── Cost function ──
    # 2026-09-17 (user decision): the allocative-efficiency wedge is BYPASSED.
    # The storyline only needs the extreme closures η ∈ {0, 1}; the wedge's
    # second-order penalty and its accounting rent in the middle ground were
    # speculative complications. The geometric interpolation of labour demand
    # (η in the labour block) is unaffected; _allocation_efficiency_wedge is
    # retained as dormant code. Consequence: at η < 1 the zero-profit block
    # prices at cost-minimizing cost while actual employment is the
    # interpolated L_i — a documented one-wage-device artifact, NOT enforced
    # by any accounting identity (the wedge rent no longer exists as a
    # separate leak; it is simply not modelled).
    cost = _ces_unit_cost(supply_shock, factor_share, w, intermediate_price, ϵ)

    # ── Equation 1: Zero-profit for ALL N sectors ──
    out[1:N] .= p .- cost

    # ── Equation 2: Market clearing for sectors 1..N-1 ──
    # The numeraire CPI = 1 fixes the gauge of the scale-homogeneous system;
    # the ray member it selects clears ALL N markets (verified post-solve as
    # the canary in verify_v3.jl). Dropping the N-th equation here keeps the
    # enforced system square with a generically nonsingular Jacobian.
    out[N+1:2N-1] .= y[1:N-1] .- intermediary_demand[1:N-1] .- total_final_demand[1:N-1]

    # ── Equation 3: Labour market (flexible-wage system: ALPHA / BETA) ──
    out[2N] = labor_market_residual(labor_closure(options), model, sum(L_i), w)

    # ── Equation 4: Numeraire constraint -- CPI = 1 (gauge fix) ──
    out[2N+1] = cpi - 1.0

    nothing
end

"""Return the exact residual vector for either mobile-labor closure.

For the `:fixed` closure the canonical vector is the REDUCED form (2N-1:
p2..pN, y1..yN); a full (p, y) or (p, y, w) vector is also accepted and
reduced internally (p1 is pinned by construction).
"""
function equilibrium_residuals(model::Model{MobileLaborCES}, X::AbstractVector)
    N = length(model.data.factor_share)
    fixed = labor_closure(model.options) isa FixedWageClosure
    Xv = collect(X)
    if fixed
        # Fixed: canonical FULL vector (2N: p1..pN, y1..yN; w = 1 pinned).
        if length(Xv) == 2N
            xr = Xv
        elseif length(Xv) == 2N + 1
            xr = Xv[1:2N]   # drop the wage component of a mobile vector
        else
            throw(DimensionMismatch("fixed closure expects a $(2N)-element vector"))
        end
        out = similar(xr)
        problem_fixed(out, xr, model)
    else
        # Mobile: canonical FULL vector (2N+1: p1..pN, y1..yN, w).
        length(Xv) == 2N + 1 || throw(DimensionMismatch(
            "mobile closure expects a $(2N+1)-element vector"))
        out = similar(Xv)
        problem(out, Xv, model)
    end
    out
end

function _equilibrium_residuals(model::Model{MobileLaborCES}, sol::Solution)
    X = labor_closure(model) isa FixedWageClosure ?
        [sol.prices_raw; sol.quantities] :
        [sol.prices_raw; sol.quantities; sol.wages_raw[1]]
    equilibrium_residuals(model, X)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Solver
# ═══════════════════════════════════════════════════════════════════════════════

"""
    problem_fixed(out, X, model)

Specialized equilibrium system for the :fixed (sticky-wage) closure.
w = 1.0 is hard-coded; unknowns are p(1:N) and y(1:N) only (2N).
Equations:
  1. N zero-profit conditions:  p_i = cost_i(p, w=1.0)
  2. N-1 market-clearing:       y_i = intermed_i + final_i  (drop last sector)
  3. 1 numeraire:               CPI = 1  (replaces the redundant market equation)
Total: 2N equations, 2N unknowns.
Employment `L_i` is computed post-solve and is not constrained to `labor_bar`.
The sticky-wage counterfactual holds the production cost at its direct CES
form. The mobile-only reallocation-efficiency correction is not added because
this fixed closure has no labor-market/accounting equation to support it.
"""
function problem_fixed(out::Vector, X::Vector, model::Model{MobileLaborCES})
    (; data, options, shocks) = model
    N = length(data.factor_share)
    w = 1.0  # sticky wage

    # FIXED-WAGE FORMULATION (2N unknowns: p1..pN, y1..yN; w = 1 pinned as the
    # sticky-wage numeraire). ALL N zero-profit and ALL N clearing equations
    # are enforced -- with the v3 homogeneous budget there is no Walras
    # redundancy at fixed w, and no sector's zero-profit may be dropped.
    p = _positive_floor(X[1:N])
    y = _positive_floor(X[N+1:2N])

    (; supply_shock, demand_shock) = shocks
    (; consumption_share, Ω_raw, factor_share, labor_share) = data
    (; θ, ϵ, σ, η) = options.elasticities

    # Intermediate goods price index
    intermediate_price = _intermediate_price(Ω_raw, p, θ)

    # CPI
    cpi = sum(consumption_share .* p .^ (1 - σ))^(1 / (1 - σ))

    # Sectoral labor demand at w=1.0
    L_i = sectoral_labor_demand(p, y, w, model)

    # Final demand (budget-consistent, financed — same hook as `problem`).
    fin = model.financing
    ds_eff = preference_weights(fin, demand_shock)
    L_sum = sum(L_i)
    total_income = w * L_sum
    # NOTE: no positivity guard here — the residual function must tolerate the
    # solver's exploration of negative-income trial points (the legacy code
    # did). The E > 0 check belongs to the post-solve validation in the
    # notebooks (headline assertions verify E = w*L - T > 0 at equilibrium).
    E = household_expenditure(fin, model, total_income, p, L_sum)
    agg = sum(consumption_share .* ds_eff .* p .^ (1 - σ))
    # v3 open economy: the household consumes (1-s)E gross (saving s·E leaks);
    # only the DOMESTIC content circulates -- the import content of household,
    # government, investment and programme demand is supplied by the external
    # account. Exports are exogenous domestic sales (no margin). The leakages
    # s·E + M now balance the injections I + X -- the identity S = I + X - M
    # holds at every equilibrium (validated in the notebooks).
    c_dom = (1 .- data.saving_rate) .* (1 .- data.import_margin) .*
            (consumption_share .* ds_eff) .* E .* p .^ (-σ) ./ agg
    total_final_demand = c_dom .+
                         (1 .- data.import_margin) .* additive_demand(fin, N) .+
                         (1 .- data.import_margin) .* (data.gov_demand .+ data.exo_demand) .+
                         data.exports_demand

    # Intermediary demand
    intermediary_demand = p .^ (-θ) .* (Ω_raw' * (p .^ ϵ .* supply_shock .^ (ϵ - 1) .* intermediate_price .^ (θ - ϵ) .* (1 .- factor_share) .* y))

    # Direct CES cost at the sticky wage (CD-limit-safe); the reallocation
    # correction is supported only by the mobile closure's labor-market
    # accounting.
    cost = _ces_unit_cost(supply_shock, factor_share, w, intermediate_price, ϵ)

    # 1. Zero-profit for ALL N sectors (p1 = 1 is NOT pinned here: the wage
    #    w = 1 is the numeraire of the sticky-wage regime, and every sector's
    #    zero-profit must hold)
    out[1:N] .= p .- cost

    # 2. Market clearing for ALL N sectors. At pinned w there is no Walras
    #    redundancy under the v3 homogeneous budget: the saving leak s·E
    #    balances the exogenous injections I + X (identity S = I + X - M).
    out[N+1:2N] .= y .- intermediary_demand .- total_final_demand

    nothing
end


"""
    _solve_fixed(model; init)

Solve the sticky-wage (w=1.0) system. Returns a `Solution` with equilibrium
prices, quantities, wage=1.0, consumption, and Tornqvist real GDP computed
from consumption (B&F metric). Employment (`sum(L_i)`) is a post-solve outcome.
"""
function _solve_fixed(model::Model{MobileLaborCES}; init=nothing)
    (; data, options, shocks) = model
    N = length(data.factor_share)

    η = options.elasticities.η
    # The fixed-wage η ≈ 1 system is homogeneous and scale-indeterminate
    # UNLESS an additive demand anchor breaks homogeneity. Only F2/F3
    # (additive public bundles) anchor scale; F1 is purely compositional
    # and does not.
    isapprox(η, 1.0; rtol=0, atol=_ETA_SCALE_INDETERMINACY_TOL) &&
        !has_additive_anchor(model.financing) &&
        throw(ArgumentError("fixed-wage η=1 has a homogeneous, scale-indeterminate equilibrium; add an additive public demand anchor (TaxFinanced / ExternalDebt), or use another η"))

    if init === nothing
        init = [ones(N); data.λ]                # p = 1, y = λ
    elseif length(init) == 2N + 1
        # Full (p, y, w) init from a mobile solution: drop the wage component
        init = init[1:2N]
    end

    # Avoid asking the nonlinear solver to differentiate an already exact
    # baseline, which some solver versions report as stalled.
    x = if maximum(abs, equilibrium_residuals(model, init)) <= 1e-12
        Float64.(init)
    else
        ProbN = NonlinearSolve.NonlinearProblem(problem_fixed, init, model)
        res = NonlinearSolve.solve(ProbN, reltol=1e-6, abstol=1e-6, maxiters=20000)
        # Quality gate = the ACTUAL residual, never the retcode. Bounded LM polish.
        x = res.u
        rmax = maximum(abs, equilibrium_residuals(model, x))
        for _ in 1:3
            rmax <= 1e-6 && break
            res = NonlinearSolve.solve(
                NonlinearSolve.NonlinearProblem(problem_fixed, x, model),
                LevenbergMarquardt(); reltol=1e-8, abstol=1e-8, maxiters=20000)
            x = res.u
            rmax = maximum(abs, equilibrium_residuals(model, x))
        end
        if rmax > 1e-6
            error("MobileLaborCES._solve_fixed did not converge: retcode = $(res.retcode), max|resid| = $rmax")
        end
        x
    end

    p = x[1:N]
    q = max.(x[N+1:2N], 0.0)
    w = 1.0  # sticky wage

    (; θ, ϵ, σ, η) = options.elasticities

    # Sectoral labor demand at equilibrium
    L_i = sectoral_labor_demand(p, q, w, model)

    # Wages vector (all 1.0 — sticky)
    wages = fill(w, N)

    # Consumption (budget-consistent, financed — must match `problem_fixed`)
    numeraire = (data.consumption_share' * p .^ (1 - σ))^(1 / (1 - σ))
    fin = model.financing
    ds_eff = preference_weights(fin, shocks.demand_shock)
    L_sum = sum(L_i)
    total_income = w * L_sum
    E = household_expenditure(fin, model, total_income, p, L_sum)
    agg = sum(data.consumption_share .* ds_eff .* p .^ (1 - σ))
    consumption = (1 .- data.saving_rate) .* (data.consumption_share .* ds_eff .* E .* p .^ (-σ)) ./ agg  # gross household consumption (saving sE leaks)

    # Real GDP: consumption Tornqvist (B&F metric). v3 base = the calibrated
    # baseline household block c0_gross (data.household_baseline) -- the exact
    # v3 baseline demand, so the index is 1 at the baseline by construction.
    # Intermediate homotopy/beta rungs may have E < 0 (negative consumption);
    # report NaN there -- those points are rejected by the callers anyway.
    base_consumption = data.household_baseline
    real_gdp_index = all(>=(0), consumption) ?
        tornqvist_quantity_index(p, consumption, ones(N), base_consumption) : NaN
    nominal_gdp = w * sum(L_i)

    return Solution(p, q, wages, consumption, numeraire, real_gdp_index, nominal_gdp, model)
end


"""
    solve(model::Model{MobileLaborCES}; init)

Solve the mobile-labor CES model. Returns a `Solution` with the equilibrium
prices, quantities, wage (as a vector of the same wage in all sectors),
consumption, and GDP measures.

The initial guess defaults to baseline prices=1, quantities=λ, and wage=1.
"""
function solve(model::Model{MobileLaborCES};
    init = nothing
)
    (; data, options, shocks) = model
    N = length(data.factor_share)

    # ── Sticky-wage path (:fixed closure) ──
    # w = 1.0 is hard-coded; solve a 2N system (prices + quantities only).
    # Employment is computed post-solve.
    if labor_closure(options) isa FixedWageClosure
        return _solve_fixed(model; init=init)
    end

    # ── Standard :mobile path (FULL 2N+1 system: p1..pN, y1..yN, w) ──
    if init === nothing
        # Default initialization: p=1, y=λ, wage=1
        init = [ones(N); data.λ; 1.0]
    end

    # As in the fixed-wage path, avoid asking the nonlinear solver to
    # differentiate an already exact baseline (some versions report this as
    # stalled). Non-exact guesses still retain retcode checks.
    x = if maximum(abs, equilibrium_residuals(model, init)) <= 1e-12
        Float64.(init)
    else
        ProbN = NonlinearSolve.NonlinearProblem(problem, init, model)
        res = NonlinearSolve.solve(ProbN, reltol=1e-6, abstol=1e-6, maxiters=20000)
        # Quality gate = the ACTUAL residual, never the retcode (NonlinearSolve
        # reports Stalled on slow final convergence). Bounded LM polish (up to
        # 3 attempts) from the last point; verify before accepting. Mobile gate:
        # 1e-5 -- at the near-singular labour-equation direction the FD-Newton
        # floor is ~3e-6 (sum L off by 0.0003 percent, economically nil); the
        # budget identities remain EXACT and the DELTA equivalence gate stays
        # at machine precision.
        x = res.u
        rmax = maximum(abs, equilibrium_residuals(model, x))
        for _ in 1:3
            rmax <= 1e-5 && break
            res = NonlinearSolve.solve(
                NonlinearSolve.NonlinearProblem(problem, x, model),
                LevenbergMarquardt(); reltol=1e-8, abstol=1e-8, maxiters=20000)
            x = res.u
            rmax = maximum(abs, equilibrium_residuals(model, x))
        end
        if rmax > 1e-5
            error("MobileLaborCES.solve did not converge: retcode = $(res.retcode), max|resid| = $rmax")
        end
        x
    end

    # Full unknowns -> Solution fields (clamp solver dust below zero; the
    # Tornqvist index requires nonnegative quantities)
    p = x[1:N]
    q = max.(x[N+1:2N], 0.0)
    w = x[2N+1]

    (; θ, ϵ, σ, η) = options.elasticities

    # Sectoral labor demand at equilibrium
    L_i = sectoral_labor_demand(p, q, w, model)

    # Wages vector (all equal to w — mobile labor)
    wages = fill(w, N)

    # Consumption — must match the budget-consistent, financed demand used
    # inside `problem` so the reported allocation is the actual equilibrium
    # allocation.
    numeraire = (data.consumption_share' * p .^ (1 - σ))^(1 / (1 - σ))
    fin = model.financing
    ds_eff = preference_weights(fin, shocks.demand_shock)
    L_sum = sum(L_i)
    total_income = w * L_sum
    E = household_expenditure(fin, model, total_income, p, L_sum)
    agg = sum(data.consumption_share .* ds_eff .* p .^ (1 - σ))
    consumption = (1 .- data.saving_rate) .* (data.consumption_share .* ds_eff .* E .* p .^ (-σ)) ./ agg  # gross household consumption (saving sE leaks)

    # Real GDP: Tornqvist (Divisia) quantity index of FINAL CONSUMPTION (value-added)
    # — not gross output. In B&F (2019), real GDP is a Divisia index of real final
    # demand / value added. Using gross output confounds intermediate flows with
    # welfare-relevant final output and dilutes the reallocation bridge. The
    # consumption (final-demand) Tornqvist isolates the welfare-relevant change.
    # v3 base = the calibrated baseline household block (data.household_baseline),
    # so the index is 1 at the v3 baseline by construction.
    base_consumption = data.household_baseline
    real_gdp_index = all(>=(0), consumption) ?
        tornqvist_quantity_index(p, consumption, ones(N), base_consumption) : NaN
    nominal_gdp = w * sum(L_i)

    return Solution(p, q, wages, consumption, numeraire, real_gdp_index, nominal_gdp, model)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Utility: construct a model from data + elasticities
# ═══════════════════════════════════════════════════════════════════════════════

"""
    mobile_labor_model(data, shocks, θ, ϵ, σ, η; labor_bar=nothing)

Convenience constructor for a MobileLaborCES model.
`closure=:mobile` uses a flexible, market-clearing wage; `closure=:fixed` uses a
sticky wage of one and unconstrained employment demand.
"""
function mobile_labor_model(data::Data, shocks::Shocks, θ::Float64, ϵ::Float64, σ::Float64, η::Float64; labor_bar::Union{Real, Nothing}=nothing, closure=:mobile, financing::Union{AbstractFinancing, Nothing}=nothing, eta_s::Union{Real, Nothing}=nothing)
    # BETA: an explicit supply elasticity forces the :beta closure.
    if eta_s !== nothing
        closure in (:mobile, :beta) || throw(ArgumentError(
            "eta_s applies to the :beta closure (got closure = $closure)"))
        closure = :beta
    end
    el = eta_s === nothing ? MobileLaborCESElasticities(θ, ϵ, σ, η) :
                             MobileLaborCESElasticities(θ, ϵ, σ, η, Float64(eta_s))
    closure_symbol = _closure_symbol(closure)
    labor_bar !== nothing && closure_symbol == :fixed && throw(ArgumentError(
        "fixed closure treats employment as an outcome; labor_bar is not used and must not be supplied"))
    lb = labor_bar === nothing ? sum(data.labor_share) : Float64(labor_bar)
    fin = financing === nothing ? NoFinancing() : financing
    model = Model(data, shocks, MobileLaborCES(el, lb, closure_symbol), fin)
    return model
end

mobile_labor_model(data::Data, shocks::Shocks, θ::Real, ϵ::Real, σ::Real, η::Real;
    labor_bar=nothing, closure=:mobile, financing=nothing, eta_s=nothing) = mobile_labor_model(data, shocks, Float64(θ), Float64(ϵ), Float64(σ), Float64(η);
    labor_bar=labor_bar === nothing ? nothing : Float64(labor_bar),
    closure=_closure_symbol(closure),
    financing=financing,
    eta_s=eta_s)
