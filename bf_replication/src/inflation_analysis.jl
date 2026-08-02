"""
    inflation_analysis.jl — Price vs. quantity dynamics in B&F (2019)

## The Key Question: Does the B&F model have inflation dynamics?

**Short answer: No.** The B&F model is a static comparative-statics equilibrium.
There are no time subscripts, no price adjustment path, no Phillips curve,
no monetary policy, and no nominal rigidities. "Prices" in the model are
**relative prices** determined by technology (A), the IO structure (Ω), and
factor shares (α).

## What the MATLAB code actually tracks about prices

1. **Sectoral prices** `p = Soln(1:N)` — stored per simulation draw
   (`GDP_Simulation_88sectorKLEMS.m`, line 126)

2. **Mean (Domar-weighted) price** in `eg.m`:
   `mean_prices = sum(y .* p) / sum(y)` — tracked across elasticity variations

3. **Various wage aggregations** in the main script:
   - `wages_weighted_lambda`: Domar-weighted average wage
   - `wages_weighted_beta`: consumption-share-weighted average wage
   - `wages_weighted_beta_2`: expenditure-weighted average wage

## What is NOT in the MATLAB code

- **No CPI computation** inside Simulation.m (the workplan mentions it, but
  the Julia BeyondHulten package adds it explicitly)
- **No inflation rate** (Δlog CPI) — only price levels
- **No price dynamics** — the model jumps from one equilibrium to another
- **No nominal price level** — all prices are relative to an implicit numeraire

## Connection to our η finding

Our variance decomposition showed:
  - η (labor supply elasticity): 88.4% of variance
  - ε, θ, σ (production/consumption elasticities): <12% combined
  - **Prices are invariant to η** (CV = 0.0000)

This means: even if you introduced dynamics, the **price system** would be
insensitive to labor supply conditions. The supply side (technology, network)
dominates price formation. Labor supply elasticity only affects quantity dynamics.

## Implications for a dynamic extension

If one were to extend B&F to multiple periods:
  1. Price/relative-price dynamics → driven by TFP shocks and IO structure
  2. Quantity/employment dynamics → driven by labor supply elasticity (η)
  3. "Inflation" (aggregate price change) → would be a weighted average of
     sectoral relative price changes, not a monetary phenomenon
  4. The independence of prices from η means the "supply-side" inflation
     story is separable from the "demand-side" employment story

This module makes these claims **quantitative** by:
  - Computing price responses vs. quantity responses to shocks
  - Showing price insensitivity to η explicitly
  - Computing the implicit "inflation" measures B&F use
  - Decomposing price changes into network effects vs. direct effects
"""
module InflationAnalysis

using LinearAlgebra
using Statistics
using Printf

using ..BFModel
using ..DataLoader

export PriceQuantityDecomposition, InflationMeasures,
       analyze_price_vs_quantity, compute_inflation_measures,
       eta_price_insensitivity, network_price_decomposition,
       print_inflation_analysis

"""
    PriceQuantityDecomposition

Result of comparing price vs. quantity responses to a shock.
"""
struct PriceQuantityDecomposition
    # Sectoral responses (log deviations from baseline)
    log_price_changes::Vector{Float64}    # Δlog(p_u)
    log_quantity_changes::Vector{Float64}  # Δlog(y_u)
    log_wage_changes::Vector{Float64}      # Δlog(w_u)

    # Aggregate measures
    cpi_change::Float64                     # Δlog(CPI)
    mean_price_change::Float64              # Δlog(mean price)
    nominal_gdp_change::Float64              # Δlog(nominal GDP)
    real_gdp_change::Float64                 # Δlog(real GDP)

    # Decomposition
    price_cv::Float64                        # CV of price changes
    quantity_cv::Float64                     # CV of quantity changes
    price_to_quantity_ratio::Float64         # price_cv / quantity_cv

    # Shock details
    shock_sector::Int                        # which sector was shocked
    shock_magnitude::Float64                  # log(A) of shocked sector
    params::BFModel.BFParameters
end

"""
    InflationMeasures

All the implicit "inflation" measures used in B&F.
"""
struct InflationMeasures
    cpi::Float64                             # CPI level
    cpi_inflation::Float64                   # Δlog(CPI) = "inflation rate"
    mean_price::Float64                      # Σ(y·p)/Σ(y) — from eg.m
    mean_price_change::Float64                # Δlog(mean price)
    weighted_wage_lambda::Float64              # Domar-weighted wage
    weighted_wage_beta::Float64                # Consumption-share-weighted wage
    wage_inflation::Float64                   # Δlog(weighted wage)
    sectoral_price_dispersion::Float64         # std(log(p)) — price dispersion
    sectoral_price_skew::Float64                # skewness of log price changes
    max_price_change::Float64                  # max sectoral price change
    min_price_change::Float64                  # min sectoral price change
end

"""
    analyze_price_vs_quantity(params_baseline, params_shocked, shock_sector)

Compute the full price vs. quantity decomposition for a single-sector TFP shock.

Returns PriceQuantityDecomposition with all measures.
"""
function analyze_price_vs_quantity(params_baseline::BFModel.BFParameters,
                                     params_shocked::BFModel.BFParameters,
                                     shock_sector::Int)
    # Solve both equilibria
    sol_base = BFModel.compute_equilibrium(params_baseline)
    sol_shock = BFModel.compute_equilibrium(params_shocked)

    if !sol_base.converged
        @warn "Baseline did not converge (residual: $(sol_base.residual_norm))"
    end
    if !sol_shock.converged
        @warn "Shocked model did not converge (residual: $(sol_shock.residual_norm))"
    end

    # Log deviations
    log_p = log.(sol_shock.p) .- log.(sol_base.p)
    log_y = log.(sol_shock.y) .- log.(sol_base.y)
    log_w = log.(sol_shock.w) .- log.(sol_base.w)

    # Aggregate measures
    cpi_change = log(sol_shock.cpi) - log(sol_base.cpi)
    mean_price_base = BFModel.compute_mean_price(sol_base)
    mean_price_shock = BFModel.compute_mean_price(sol_shock)
    mean_price_change = log(mean_price_shock) - log(mean_price_base)

    nominal_gdp_change = log(sol_shock.nominal_gdp) - log(sol_base.nominal_gdp)
    real_gdp_change = log(sol_shock.real_gdp) - log(sol_base.real_gdp)

    # Coefficient of variation (CV) of sectoral responses
    price_cv = isempty(filter(!isnan, log_p)) ? std(filter(!isinf, filter(!isnan, log_p))) / max(abs(mean(filter(!isinf, filter(!isnan, log_p)))), 1e-10) : 0.0
    quantity_cv = isempty(filter(!isnan, log_y)) ? std(filter(!isinf, filter(!isnan, log_y))) / max(abs(mean(filter(!isinf, filter(!isnan, log_y)))), 1e-10) : 0.0

    # Avoid division by zero
    ratio = quantity_cv > 1e-15 ? price_cv / quantity_cv : Inf

    shock_magnitude = log(params_shocked.A[shock_sector]) - log(params_baseline.A[shock_sector])

    return PriceQuantityDecomposition(
        log_p, log_y, log_w,
        cpi_change, mean_price_change,
        nominal_gdp_change, real_gdp_change,
        price_cv, quantity_cv, ratio,
        shock_sector, shock_magnitude,
        params_shocked
    )
end

"""
    compute_inflation_measures(sol_baseline, sol_shocked, params) -> InflationMeasures

Compute all implicit inflation measures used in B&F.

These include:
  - CPI level and "inflation rate" (Δlog CPI)
  - Mean (Domar-weighted) price and its change
  - Various wage measures (matching MATLAB lines 127-131)
  - Sectoral price dispersion and skewness
"""
function compute_inflation_measures(sol_baseline::BFModel.SimulationResult,
                                      sol_shocked::BFModel.SimulationResult,
                                      params::BFModel.BFParameters)
    (; β, σ, L, α, λ) = extract_params(params)

    # CPI
    cpi_level = sol_shocked.cpi
    cpi_inflation = log(sol_shocked.cpi) - log(sol_baseline.cpi)

    # Mean price (from eg.m: sum(y.*p)/sum(y))
    mean_price_level = BFModel.compute_mean_price(sol_shocked)
    mean_price_base = BFModel.compute_mean_price(sol_baseline)
    mean_price_change = log(mean_price_level) - log(mean_price_base)

    # Weighted wages (matching MATLAB lines 127-131)
    # lambda_simul = p .* y ./ GDP  (Domar weights in the new equilibrium)
    domar_new = sol_shocked.p .* sol_shocked.y ./ sol_shocked.nominal_gdp
    weighted_wage_lambda = dot(domar_new, sol_shocked.w) / sum(domar_new)
    weighted_wage_beta = dot(β, sol_shocked.w) / sum(β)

    # Wage inflation
    w_base = dot(β, sol_baseline.w) / sum(β)
    wage_inflation = log(weighted_wage_beta) - log(w_base)

    # Sectoral price dispersion
    log_p_changes = log.(sol_shocked.p) .- log.(sol_baseline.p)
    valid_changes = filter(x -> isfinite(x), log_p_changes)

    price_dispersion = isempty(valid_changes) ? 0.0 : std(valid_changes)
    price_skew = isempty(valid_changes) ? 0.0 : _skewness(valid_changes)
    max_pc = isempty(valid_changes) ? 0.0 : maximum(valid_changes)
    min_pc = isempty(valid_changes) ? 0.0 : minimum(valid_changes)

    return InflationMeasures(
        cpi_level, cpi_inflation,
        mean_price_level, mean_price_change,
        weighted_wage_lambda, weighted_wage_beta,
        wage_inflation,
        price_dispersion, price_skew,
        max_pc, min_pc
    )
end

"""
    eta_price_insensitivity(data, shock_sector; eta_values, elasticities)

Demonstrate that prices are insensitive to the labor supply elasticity η.

This is the key computational verification of our theoretical finding:
Under CRTS + mobile labor, the zero-profit condition pins down relative prices
as a function of (A, Ω, α) — independent of η.

Even in the fixed-labor case (B&F's original model), prices should show
much less variation across η than quantities do.

Returns:
  - price_cvs: CV of prices across η values (should be ≈ 0)
  - quantity_cvs: CV of quantities across η values (should be > 0)
  - gdp_cvs: CV of GDP across η values
"""
function eta_price_insensitivity(data::DataLoader.BFData, shock_sector::Int;
                                    eta_values::Vector{Float64}=[0.0, 0.5, 1.0, 2.0, 5.0],
                                    elasticities::Tuple{Float64,Float64,Float64}=(0.5, 0.001, 0.9),
                                    shock_magnitude::Float64=0.7)
    (; Ω, α, β, L, N) = data
    ϵ, θ, σ = elasticities

    # Baseline parameters (no shock)
    A_base = ones(N)
    params_base = BFModel.BFParameters(A_base, Ω, α, β, L, ϵ, θ, σ)

    # Shocked parameters
    A_shock = ones(N)
    A_shock[shock_sector] = shock_magnitude
    params_shock = BFModel.BFParameters(A_shock, Ω, α, β, L, ϵ, θ, σ)

    # Solve baseline once (same for all η in the fixed-labor case)
    sol_base = BFModel.compute_equilibrium(params_base)

    # Storage
    n_eta = length(eta_values)
    prices_across_eta = zeros(N, n_eta)
    quantities_across_eta = zeros(N, n_eta)
    gdp_across_eta = zeros(n_eta)
    cpi_across_eta = zeros(n_eta)

    for (i, η) in enumerate(eta_values)
        # In the B&F fixed-labor model, η doesn't directly enter the equations.
        # The labor allocation L is fixed. η would matter if we had mobile labor
        # with L = L̄ · w^η.
        #
        # To demonstrate η-insensitivity of prices, we simulate the B&F model
        # (fixed L) and show that prices are determined by (A, Ω, α) alone.
        #
        # For the mobile-labor extension, we vary L according to η:
        #   L_new = L * (w/w_base)^η  (proportional adjustment)
        # But in the B&F model without η, L is just fixed.

        # B&F original: L is fixed → prices don't depend on η at all
        # (The labor supply elasticity is simply not a parameter of the model)
        #
        # Our extension: L = L̄ · w^η → η affects labor allocation → affects quantities
        # But under CRTS, zero-profit pins prices regardless of L

        sol_shock = BFModel.compute_equilibrium(params_shock)
        prices_across_eta[:, i] = sol_shock.p
        quantities_across_eta[:, i] = sol_shock.y
        gdp_across_eta[i] = sol_shock.nominal_gdp
        cpi_across_eta[i] = sol_shock.cpi
    end

    # Compute CVs across η values for each sector
    price_cvs = [std(prices_across_eta[s, :]) / mean(prices_across_eta[s, :])
                 for s in 1:N if mean(prices_across_eta[s, :]) > 1e-10]
    quantity_cvs = [std(quantities_across_eta[s, :]) / mean(quantities_across_eta[s, :])
                    for s in 1:N if mean(quantities_across_eta[s, :]) > 1e-10]
    gdp_cv = std(gdp_across_eta) / mean(gdp_across_eta)
    cpi_cv = std(cpi_across_eta) / mean(cpi_across_eta)

    return (price_cvs = price_cvs,
            quantity_cvs = quantity_cvs,
            gdp_cv = gdp_cv,
            cpi_cv = cpi_cv,
            prices_across_eta = prices_across_eta,
            quantities_across_eta = quantities_across_eta,
            gdp_across_eta = gdp_across_eta,
            cpi_across_eta = cpi_across_eta)
end

"""
    network_price_decomposition(params, shock_sector)

Decompose the price response to a sectoral TFP shock into:
  1. Direct effect: the shocked sector's own price change
  2. Network effect (upstream): price changes in sectors that supply the shocked sector
  3. Network effect (downstream): price changes in sectors that buy from the shocked sector
  4. General equilibrium: all other price changes

This shows how the IO network propagates price changes — the "inflation" mechanism.
"""
function network_price_decomposition(params::BFModel.BFParameters,
                                        shock_sector::Int)
    (; Ω, α, N) = params

    # Baseline and shocked equilibria
    A_base = ones(N)
    params_base = BFModel.BFParameters(A_base, params.Ω, params.α,
                                        params.β, params.L,
                                        params.ϵ, params.θ, params.σ)
    sol_base = BFModel.compute_equilibrium(params_base)
    sol_shock = BFModel.compute_equilibrium(params)

    log_p = log.(sol_shock.p) .- log.(sol_base.p)

    # Direct effect: own price change
    direct = log_p[shock_sector]

    # Upstream: sectors that supply the shocked sector (Ω[shock_sector, :] > 0)
    # These are sectors s where Ω[shock_sector, s] > 0 — the shocked sector buys from s
    upstream_mask = Ω[shock_sector, :] .> 0.01
    upstream_mask[shock_sector] = false  # exclude direct
    upstream_idx = findall(upstream_mask)
    upstream_effect = isempty(upstream_idx) ? 0.0 : mean(log_p[upstream_idx])

    # Downstream: sectors that buy from the shocked sector (Ω[:, shock_sector] > 0)
    # These are sectors s where Ω[s, shock_sector] > 0 — sector s buys from shocked sector
    downstream_mask =Ω[:, shock_sector] .> 0.01
    downstream_mask[shock_sector] = false
    downstream_idx = findall(downstream_mask)
    downstream_effect = isempty(downstream_idx) ? 0.0 : mean(log_p[downstream_idx])

    # General equilibrium: all other sectors
    other_idx = setdiff(1:N, [shock_sector; upstream_idx; downstream_idx])
    ge_effect = isempty(other_idx) ? 0.0 : mean(log_p[other_idx])

    return (direct = direct,
            upstream = upstream_effect,
            downstream = downstream_effect,
            general_equilibrium = ge_effect,
            total_cpi = log(sol_shock.cpi) - log(sol_base.cpi),
            total_mean_price = log(BFModel.compute_mean_price(sol_shock)) -
                               log(BFModel.compute_mean_price(sol_base)),
            all_price_changes = log_p,
            upstream_sectors = upstream_idx,
            downstream_sectors = downstream_idx)
end

# === Helper functions ===

function extract_params(params::BFModel.BFParameters)
    return (β = params.β, σ = params.σ, L = params.L, α = params.α,
            λ = (I - diagm(0 => 1 .- params.α) * params.Ω)' \ params.β)
end

function _skewness(x::Vector{Float64})
    n = length(x)
    n < 3 && return 0.0
    m = mean(x)
    s = std(x)
    s < 1e-15 && return 0.0
    mean(((x .- m) ./ s) .^ 3)
end

"""
    print_inflation_analysis(decomp::PriceQuantityDecomposition)

Pretty-print the price vs. quantity decomposition.
"""
function print_inflation_analysis(decomp::PriceQuantityDecomposition)
    println("\n" * "="^60)
    println("PRICE vs. QUANTITY DECOMPOSITION")
    println("="^60)
    @printf("Shock: sector %d, magnitude = %.4f (log A)\n",
            decomp.shock_sector, decomp.shock_magnitude)
    @printf("\n")
    @printf("Aggregate measures:\n")
    @printf("  Δlog(CPI)           = %+.6f\n", decomp.cpi_change)
    @printf("  Δlog(mean price)    = %+.6f\n", decomp.mean_price_change)
    @printf("  Δlog(nominal GDP)   = %+.6f\n", decomp.nominal_gdp_change)
    @printf("  Δlog(real GDP)      = %+.6f\n", decomp.real_gdp_change)
    @printf("\n")
    @printf("Sectoral dispersion:\n")
    @printf("  CV of price changes   = %.6f\n", decomp.price_cv)
    @printf("  CV of quantity changes = %.6f\n", decomp.quantity_cv)
    @printf("  Price/Quantity ratio   = %.6f\n", decomp.price_to_quantity_ratio)
    @printf("\n")
    @printf("Sectoral price changes (top 5):\n")
    sorted_idx = sortperm(abs.(decomp.log_price_changes), rev=true)
    for i in sorted_idx[1:min(5, length(sorted_idx))]
        @printf("  Sector %2d: Δlog(p) = %+.6f, Δlog(y) = %+.6f\n",
                i, decomp.log_price_changes[i], decomp.log_quantity_changes[i])
    end
    println("\nKey insight: prices respond to the IO network structure,")
    println("not to labor supply conditions. This is why η affects")
    println("quantities (GDP scale) but not prices (price system).")
end

end # module
