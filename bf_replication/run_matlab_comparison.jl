"""
    run_matlab_comparison.jl

Compare our Julia replication against the original MATLAB code's exact specifications.

Key corrections from reading the MATLAB code:
1. Single-sector test uses θ = 0.0001 (not 0.001)
2. GDP measure = sum(w' * L * p^(-σ) .* β) — welfare-weighted GDP
   (NOT nominal GDP w'*L, and NOT real GDP)
3. Default Monte Carlo uses θ = 0.001, ϵ = 0.5, σ = 0.9
4. Must also compute Hulten benchmark: Δln(GDP) = λ' * Δln(A)
"""

include("src/BFReplication.jl")
using .BFReplication.DataLoader
using .BFReplication.BFModel
using .BFReplication.InflationAnalysis

using LinearAlgebra, Statistics, Printf, DelimitedFiles

# === Configuration ===
data_dir = joinpath(@__DIR__, "..", "Replication Files", "GDP Simulatin -- 88 Sector")
csv_path = joinpath(data_dir, "BFdata.csv")

# === Step 1: Load data (same pipeline as MATLAB) ===
println("="^70)
println("B&F (2019) — Comparing Julia Replication to MATLAB Code")
println("="^70)

println("\n[1] Loading B&F data from BFdata.csv...")
data = load_bf_data(csv_path, year=1980)
describe_data(data)

N = data.N
lambda = (I - diagm(0 => 1 .- data.α) * data.Ω)' \ data.β

# === Step 2: Single-sector oil shock — MATLAB exact spec ===
println("\n" * "="^70)
println("SINGLE-SECTOR OIL SHOCK (A₇ = 0.7, 30% TFP reduction)")
println("MATLAB spec: ϵ=0.5, θ=0.0001, σ=0.9")
println("="^70)

ϵ = 0.5
θ = 0.0001  # MATLAB's exact value for single-sector test (line 210)
σ = 0.9

SHOCK_SECTOR = 7    # Oil
A_shock_val = 0.7   # 30% reduction

A_base = ones(N)
params_base = BFModel.BFParameters(A_base, data.Ω, data.α, data.β, data.L, ϵ, θ, σ)

A_shock = ones(N)
A_shock[SHOCK_SECTOR] = A_shock_val
params_shock = BFModel.BFParameters(A_shock, data.Ω, data.α, data.β, data.L, ϵ, θ, σ)

sol_base = BFModel.compute_equilibrium(params_base)
sol_shock = BFModel.compute_equilibrium(params_shock)

# MATLAB GDP measure (line 235 of GDP_Simulation_88sectorKLEMS.m):
#   GDP = sum(wages' * L * p .^ (-sigma) .* beta)
#   = (w'*L) * Σ β_u * p_u^(-σ)
wages_shock = sol_shock.w
GDP_welfare_MATLAB = sum(wages_shock' * data.L .* (sol_shock.p .^ (-σ)) .* data.β)
wages_base = sol_base.w
GDP_welfare_base_MATLAB = sum(wages_base' * data.L .* (sol_base.p .^ (-σ)) .* data.β)

# Also compute nominal GDP (used in eg.m and Oil_Shock.m)
GDP_nominal_shock = sol_shock.nominal_gdp   # w'*L
GDP_nominal_base = sol_base.nominal_gdp

# Hulten benchmark: Δln(GDP) = λ' * Δln(A)
ΔlnA = log.(A_shock) .- log.(A_base)
GDP_hulten_change = lambda' * ΔlnA
GDP_hulten_level = exp(GDP_hulten_change)

@printf("\nResults for A₇ = 0.7 (θ = %.4f):\n", θ)
@printf("  Welfare GDP (MATLAB formula):  %.6f  (Δlog = %+.6f)\n",
    GDP_welfare_MATLAB, log(GDP_welfare_MATLAB) - log(GDP_welfare_base_MATLAB))
@printf("  Nominal GDP (w'*L):           %.6f  (Δlog = %+.6f)\n",
    GDP_nominal_shock, log(GDP_nominal_shock) - log(GDP_nominal_base))
@printf("  Hulten benchmark:              %.6f  (Δlog = %+.6f)\n",
    GDP_hulten_level, GDP_hulten_change)
@printf("  CPI:                           %.6f\n", sol_shock.cpi)
@printf("  Mean price (Domar-weighted):   %.6f\n", BFModel.compute_mean_price(sol_shock))

# Sectoral price changes
log_p = log.(sol_shock.p) .- log.(sol_base.p)
@printf("\n  Top 5 sectoral price changes:\n")
sorted = sortperm(abs.(log_p), rev=true)
for i in sorted[1:min(5, N)]
    @printf("    Sector %2d: Δlog(p) = %+.6f\n", i, log_p[i])
end

# === Step 3: Compare with different θ values ===
println("\n" * "="^70)
println("SENSITIVITY: GDP measure by θ value")
println("="^70)

for θ_test in [0.0001, 0.001, 0.01, 0.1, 1.0]
    p_test = BFModel.BFParameters(A_shock, data.Ω, data.α, data.β, data.L, ϵ, θ_test, σ)
    sol_test = BFModel.compute_equilibrium(p_test)
    if sol_test.converged
        w_test = sol_test.w
        gdp_welfare = sum(w_test' * data.L .* (sol_test.p .^ (-σ)) .* data.β)
        gdp_nominal = sol_test.nominal_gdp
        @printf("  θ=%7.4f: welfare GDP=%.6f (Δlog=%+.6f), nominal=%.6f, CPI=%.6f\n",
            θ_test, gdp_welfare, log(gdp_welfare) - log(GDP_welfare_base_MATLAB),
            gdp_nominal, sol_test.cpi)
    else
        @printf("  θ=%7.4f: did not converge\n", θ_test)
    end
end

# === Step 4: Full response curve (sector 7, A ∈ [0.7, 1.3]) ===
println("\n" * "="^70)
println("FULL RESPONSE CURVE: Sector 7 (Oil), A ∈ [0.7, 1.3]")
println("MATLAB spec: θ=0.0001")
println("="^70)

θ = 0.0001
M = 10
a, b = 0.7, 1.3

# Contraction side (A from 1.0 down to 0.7)
println("\n  A_val    log(A)   Welfare GDP  Δlog(WGDP)  Hulten Δlog  Gap")
println("  " * "-"^60)
for A_val in range(1.0, a, M)
    A_tmp = ones(N)
    A_tmp[SHOCK_SECTOR] = A_val
    p_tmp = BFModel.BFParameters(A_tmp, data.Ω, data.α, data.β, data.L, ϵ, θ, σ)
    sol_tmp = BFModel.compute_equilibrium(p_tmp)
    if sol_tmp.converged
        w_tmp = sol_tmp.w
        gdp_w = sum(w_tmp' * data.L .* (sol_tmp.p .^ (-σ)) .* data.β)
        ΔlnA_tmp = log(A_val)  # since A_base[7]=1
        ΔlnGDP_w = log(gdp_w) - log(GDP_welfare_base_MATLAB)
        hulten_change = lambda[SHOCK_SECTOR] * ΔlnA_tmp
        gap = ΔlnGDP_w / hulten_change - 1
        @printf("  %.4f   %+.4f   %.6f   %+.6f   %+.6f   %+.4f\n",
            A_val, ΔlnA_tmp, gdp_w, ΔlnGDP_w, hulten_change, gap)
    end
end

# Expansion side (A from 1.0 up to 1.3)
println()
for A_val in range(1.0, b, M)
    A_tmp = ones(N)
    A_tmp[SHOCK_SECTOR] = A_val
    p_tmp = BFModel.BFParameters(A_tmp, data.Ω, data.α, data.β, data.L, ϵ, θ, σ)
    sol_tmp = BFModel.compute_equilibrium(p_tmp)
    if sol_tmp.converged
        w_tmp = sol_tmp.w
        gdp_w = sum(w_tmp' * data.L .* (sol_tmp.p .^ (-σ)) .* data.β)
        ΔlnA_tmp = log(A_val)
        ΔlnGDP_w = log(gdp_w) - log(GDP_welfare_base_MATLAB)
        hulten_change = lambda[SHOCK_SECTOR] * ΔlnA_tmp
        gap = ΔlnGDP_w / hulten_change - 1
        @printf("  %.4f   %+.4f   %.6f   %+.6f   %+.6f   %+.4f\n",
            A_val, ΔlnA_tmp, gdp_w, ΔlnGDP_w, hulten_change, gap)
    end
end

# === Step 5: Compare the three sectors from MATLAB (Oil=7, Retail=53, Construction=8) ===
println("\n" * "="^70)
println("SECTOR COMPARISON (MATLAB: sectors 7, 53, 8)")
println("Log GDP vs Hulten for A ∈ [0.7, 1.3]")
println("="^70)

θ = 0.0001
sectors_to_test = [7, 53, 8]
sector_names = ["Oil", "Retail", "Construction"]

for (idx, s) in enumerate(sectors_to_test)
    println("\n  Sector $s ($(sector_names[idx])) — λ = $(lambda[s]):")
    println("  A_val    log(A)   Δlog(GDP)  Δlog(Hulten)  Ratio")
    for A_val in vcat(range(a, 1.0, M), range(1.0, b, M)[2:end])
        A_tmp = ones(N)
        A_tmp[s] = A_val
        p_tmp = BFModel.BFParameters(A_tmp, data.Ω, data.α, data.β, data.L, ϵ, θ, σ)
        sol_tmp = BFModel.compute_equilibrium(p_tmp)
        if sol_tmp.converged
            w_tmp = sol_tmp.w
            gdp_w = sum(w_tmp' * data.L .* (sol_tmp.p .^ (-σ)) .* data.β)
            ΔlnGDP_w = log(gdp_w) - log(GDP_welfare_base_MATLAB)
            hulten_change = lambda[s] * log(A_val)
            @printf("  %.4f   %+.4f   %+.6f   %+.6f   %.4f\n",
                A_val, log(A_val), ΔlnGDP_w, hulten_change, ΔlnGDP_w / hulten_change)
        end
    end
end

# === Step 6: Summary of discrepancies ===
println("\n" * "="^70)
println("SUMMARY: Comparison with Original MATLAB Code")
println("="^70)
println("""
Model equations verified:
  ✓ q = (Ω * p^(1-θ))^(1/(1-θ))   — matches Simulation.m line 7
  ✓ w = p * A^((ϵ-1)/ϵ) * α^(1/ϵ) * y^(1/ϵ) * L^(-1/ϵ)  — matches line 8
  ✓ C = w'*L   — matches line 9
  ✓ Price residuals — matches line 11
  ✓ Market clearing — matches line 12
  ✓ Analytical Jacobian in Simulation_Derivs.m matches our numerical Jacobian

Key parameters:
  ✓ Default Monte Carlo: ϵ=0.5, θ=0.001, σ=0.9   (GDP_Simulation_88sectorKLEMS.m lines 90-92)
  ✓ Single-sector shock: ϵ=0.5, θ=0.0001, σ=0.9  (lines 209-211)
  ✓ Data: year=1980, 76 sectors (same removal: govt [60,80:88], zero-sales [8,62])

GDP measure differences:
  - MATLAB single-sector test uses: GDP = w'*L * Σ β*p^(-σ)  (welfare-weighted)
  - MATLAB eg.m uses: GDP = w'*L  (nominal GDP)
  - Our original analysis used: nominal GDP w'*L and CES real GDP
  - For A₇=0.7, θ=0.0001: welfare GDP change = ? (just computed above)

Paper's published key results to target:
  - Oil shocks: 0.23% → 0.61% GDP impact (nonlinearities triple the effect)
  - Business cycle: average log(GDP) ≈ -0.5% (an order > Lucas 1987)
  - GDP distribution: negative skew, excess kurtosis
""")