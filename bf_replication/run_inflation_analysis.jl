"""
    run_inflation_analysis.jl — Test script for inflation dynamics analysis

This script:
  1. Loads the B&F data from simulationData.mat
  2. Solves the baseline equilibrium
  3. Applies a TFP shock to a key sector (e.g., Oil = sector 7)
  4. Analyzes price vs. quantity responses
  5. Computes all implicit "inflation" measures
  6. Tests η-insensitivity of prices
  7. Decomposes price changes into network effects
"""

using LinearAlgebra
using Statistics
using Printf

# Include modules
include("src/BFReplication.jl")
using .BFReplication
using .BFReplication.DataLoader
using .BFReplication.BFModel
using .BFReplication.InflationAnalysis

# === Configuration ===
data_dir = joinpath(@__DIR__, "..", "Replication Files", "GDP Simulatin -- 88 Sector")
csv_path = joinpath(data_dir, "BFdata.csv")
tfp_path = joinpath(data_dir, "stfp.csv")

# Elasticities from B&F (default calibration)
ϵ = 0.5    # elasticity: intermediates vs labor
θ = 0.001  # elasticity: between intermediates (near-Leontief)
σ = 0.9    # elasticity: consumption substitution

# Shock sector: Oil (sector 7 in B&F's 76-sector data)
SHOCK_SECTOR = 7
SHOCK_MAGNITUDE = 0.7  # 30% productivity reduction

println("\n" * "="^60)
println("B&F (2019) REPLICATION — INFLATION DYNAMICS ANALYSIS")
println("="^60)

# === Step 1: Load data ===
println("\n[1] Loading B&F data from BFdata.csv...")
data = load_bf_data(csv_path, year=1980)
describe_data(data)

# === Step 2: Baseline equilibrium ===
println("\n[2] Solving baseline equilibrium (A = 1, no shocks)...")
A_base = ones(data.N)
params_base = BFParameters(A_base, data.Ω, data.α, data.β, data.L, ϵ, θ, σ)
sol_base = compute_equilibrium(params_base)

println("  Converged: $(sol_base.converged)")
println("  Residual:  $(sol_base.residual_norm)")
println("  CPI:       $(sol_base.cpi)")
println("  Nom. GDP:  $(sol_base.nominal_gdp)")
println("  Real GDP:  $(sol_base.real_gdp)")
println("  Mean price: $(compute_mean_price(sol_base))")

# Check baseline: all prices should be ≈ 1
max_price_dev = maximum(abs.(sol_base.p .- 1.0))
println("  Max |p - 1|: $max_price_dev")

# === Step 3: Apply TFP shock ===
println("\n[3] Applying TFP shock to sector $SHOCK_SECTOR (magnitude=$SHOCK_MAGNITUDE)...")
A_shock = ones(data.N)
A_shock[SHOCK_SECTOR] = SHOCK_MAGNITUDE
params_shock = BFParameters(A_shock, data.Ω, data.α, data.β, data.L, ϵ, θ, σ)
sol_shock = compute_equilibrium(params_shock)

println("  Converged: $(sol_shock.converged)")
println("  Residual:  $(sol_shock.residual_norm)")
println("  CPI:       $(sol_shock.cpi)")
println("  Nom. GDP:  $(sol_shock.nominal_gdp)")
println("  Real GDP:  $(sol_shock.real_gdp)")

# === Step 4: Price vs. Quantity decomposition ===
println("\n[4] Price vs. Quantity decomposition...")
decomp = analyze_price_vs_quantity(params_base, params_shock, SHOCK_SECTOR)
print_inflation_analysis(decomp)

# === Step 5: Inflation measures ===
println("\n[5] Computing all implicit inflation measures...")
inflation = compute_inflation_measures(sol_base, sol_shock, params_shock)

println("\n  Inflation Measures (B&F implicit):")
@printf("  CPI level:               %.6f\n", inflation.cpi)
@printf("  CPI inflation (Δlog):   %+.6f\n", inflation.cpi_inflation)
@printf("  Mean price level:       %.6f\n", inflation.mean_price)
@printf("  Mean price change:      %+.6f\n", inflation.mean_price_change)
@printf("  Wage (λ-weighted):      %.6f\n", inflation.weighted_wage_lambda)
@printf("  Wage (β-weighted):      %.6f\n", inflation.weighted_wage_beta)
@printf("  Wage inflation (Δlog):  %+.6f\n", inflation.wage_inflation)
@printf("  Price dispersion (σ):    %.6f\n", inflation.sectoral_price_dispersion)
@printf("  Price skewness:           %.6f\n", inflation.sectoral_price_skew)
@printf("  Max sectoral price chg:  %+.6f\n", inflation.max_price_change)
@printf("  Min sectoral price chg:  %+.6f\n", inflation.min_price_change)

# === Step 6: η price insensitivity ===
println("\n[6] Testing η-insensitivity of prices...")
println("  (In B&F's original model, η is not a parameter — L is fixed)")
println("  (Our extension: L = L̄·w^η → η affects labor supply)")

eta_result = eta_price_insensitivity(data, SHOCK_SECTOR,
    eta_values=[0.0, 0.5, 1.0, 2.0, 5.0],
    elasticities=(ϵ, θ, σ),
    shock_magnitude=SHOCK_MAGNITUDE)

println("\n  η-insensitivity results:")
@printf("  Price CV across η:    %.8f (should be ≈ 0)\n", mean(eta_result.price_cvs))
@printf("  Quantity CV across η:  %.8f\n", mean(eta_result.quantity_cvs))
@printf("  GDP CV across η:      %.8f\n", eta_result.gdp_cv)
@printf("  CPI CV across η:      %.8f (should be ≈ 0)\n", eta_result.cpi_cv)

println("\n  Interpretation:")
println("  In the B&F fixed-labor model, η simply doesn't enter the equations,")
println("  so prices are trivially invariant. In our mobile-labor extension,")
println("  the zero-profit condition under CRTS pins prices as a function of")
println("  (A, Ω, α) alone, independent of the labor supply elasticity.")

# === Step 7: Network price decomposition ===
println("\n[7] Network price decomposition...")
net_decomp = network_price_decomposition(params_shock, SHOCK_SECTOR)

println("\n  Network effects on price system:")
@printf("  Direct (own sector):       %+.6f\n", net_decomp.direct)
@printf("  Upstream (suppliers):       %+.6f (avg of %d sectors)\n",
    net_decomp.upstream, length(net_decomp.upstream_sectors))
@printf("  Downstream (customers):     %+.6f (avg of %d sectors)\n",
    net_decomp.downstream, length(net_decomp.downstream_sectors))
@printf("  General equilibrium:        %+.6f\n", net_decomp.general_equilibrium)
@printf("  Total CPI change:           %+.6f\n", net_decomp.total_cpi)
@printf("  Total mean price change:    %+.6f\n", net_decomp.total_mean_price)

# === Summary ===
println("\n" * "="^60)
println("SUMMARY: INFLATION DYNAMICS IN B&F (2019)")
println("="^60)
println("""
The B&F model is a static comparative-statics equilibrium:
  • No time dynamics, no Phillips curve, no monetary policy
  • Prices are RELATIVE prices determined by (A, Ω, α)
  • "Inflation" = aggregate relative price change, not monetary phenomenon

Key findings:
  1. The MATLAB code does NOT compute CPI inside the solver
  2. The only "inflation" measure in eg.m is a Domar-weighted mean price
  3. Prices respond to the IO network structure, not to labor supply
  4. Under CRTS, zero-profit conditions pin prices independent of η
  5. η only affects the SCALE of activity (GDP, employment)

Implication for dynamic extension:
  • Price dynamics → driven by TFP + IO network
  • Quantity dynamics → driven by labor supply elasticity
  • These are separable: the "supply-side inflation story"
    is independent of the "demand-side employment story"
""")
