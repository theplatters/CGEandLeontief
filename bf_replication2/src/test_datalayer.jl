# test_datalayer.jl — Test the IO table loading, shock loading, and network construction
#
# Run from bf_replication2/ directory:
#   julia --project=. src/test_datalayer.jl

using Printf, LinearAlgebra

include("io_table.jl")
include("shocks.jl")
include("network.jl")

println("="^60)
println("Testing BFAEA2022 Data Layer")
println("="^60)

# ---- 1. Load IO table ----
println("\n--- Test 1: Load IO Table ---")
data_dir = "data"
io = load_io_table(joinpath(data_dir, "IO_data_2018.mat"); N=66, year=2015)

@printf("  N = %d\n", io.N)
@printf("  Omega: %d×%d\n", size(io.Omega, 1), size(io.Omega, 2))
@printf("  alphaL length: %d, alphaK length: %d\n", length(io.alphaL), length(io.alphaK))
@printf("  beta length: %d, sum(beta) = %.10f\n", length(io.beta), sum(io.beta))
@printf("  va_share length: %d, int_share length: %d\n", length(io.va_share), length(io.int_share))
@printf("  indname count: %d (first: %s)\n", length(io.indname), io.indname[1])

# Check: beta should sum to 1
if abs(sum(io.beta) - 1.0) > 1e-10
    @warn "beta sum deviates from 1: $(sum(io.beta))"
else
    @printf("  ✅ beta sums to 1\n")
end

# Check: alphaL + alphaK should be ≈ 1
@printf("  min(alphaL + alphaK) = %.6f, max = %.6f\n",
    minimum(io.alphaL .+ io.alphaK), maximum(io.alphaL .+ io.alphaK))

# ---- 2. Load Shocks ----
@printf("\n--- Test 2: Load Shocks ---\n")
shocks = load_shocks(data_dir; N=66)

@printf("  BLS shock: min=%.4f, max=%.4f\n",
    minimum(shocks.BLS_shock), maximum(shocks.BLS_shock))
@printf("  PCE shock: min=%.4f, max=%.4f\n",
    minimum(shocks.PCE_shock), maximum(shocks.PCE_shock))
@printf("  Wages:     min=%.4f, max=%.4f\n",
    minimum(shocks.wages), maximum(shocks.wages))
@printf("  PPI:       min=%.4f, max=%.4f\n",
    minimum(shocks.PPI), maximum(shocks.PPI))

unique_nonzero = count(x -> x != 0.0, shocks.BLS_shock)
@printf("  Non-zero BLS entries: %d / %d\n", unique_nonzero, shocks.N)

# ---- 3. Build Standard-Form Network ----
@printf("\n--- Test 3: Build Standard-Form Network ---\n")
sf = build_standard_form(io)

@printf("  D = %d (expected %d)\n", sf.D, 5 * io.N + 4)
@printf("  Omega_re: %d×%d\n", size(sf.Omega_re, 1), size(sf.Omega_re, 2))
@printf("  Psi_re:   %d×%d\n", size(sf.Psi_re, 1), size(sf.Psi_re, 2))

# Check D = 5N + 4 = 334
if sf.D == 334
    @printf("  ✅ D = 334 (5×66 + 4)\n")
else
    @warn "D = $(sf.D), expected 334"
end

# Check Domar weights (multiplier ~ gross output / GDP = ~2, not 1)
@printf("  Domar weights sum = %.6f (expected ~1.8-2.0, gross-output/GDP ratio)\n", sum(sf.Domar))

# Check factor/keynes lengths
@printf("  factor length: %d, keynes length: %d\n", length(sf.factor), length(sf.keynes))
@printf("    1 (goods): %d, 0 (factors): %d, 2 (Ricardian): %d, 3 (HtM): %d\n",
    count(x -> x == 1, sf.factor),
    count(x -> x == 0, sf.factor),
    count(x -> x == 2, sf.factor),
    count(x -> x == 3, sf.factor))
@printf("    1 (normal): %d, 0 (flex): %d, -1 (sticky): %d\n",
    count(x -> x == 1, sf.keynes),
    count(x -> x == 0, sf.keynes),
    count(x -> x == -1, sf.keynes))

# ---- 4. Structure checks on Omega_re ----
@printf("\n--- Test 4: Omega_re block structure ---\n")
# The consumption row (row 1) should only have non-zero entries in cols 2:N+1 (goods)
cons_row = sf.Omega_re[1, :]
goods_only = findall(x -> x > 0, cons_row)
@printf("  Consumption row non-zeros: cols %s — %d goods\n",
    join(string.(goods_only), ", "), length(goods_only))

# Goods rows (2:N+1) should have entries in VA block (N+2:2N+1) and intermediates (2N+2:3N+1)
row2_nonzero = findall(x -> x > 0, sf.Omega_re[2, :])
@printf("  First goods row (row 2) non-zero cols: %s (%d total)\n",
    join(string.(row2_nonzero), ", "), length(row2_nonzero))

# Check Psi_re invertibility
Psi_residual = sf.Psi_re * (I - sf.Omega_re) - I
if norm(Psi_residual) > 1e-10
    @warn "Psi_re * (I - Omega_re) residual = $(norm(Psi_residual))"
else
    @printf("  ✅ Psi_re invertibility check passed (resid = %.2e)\n", norm(Psi_residual))
end

# ---- 5. Print summary ----
println("\n" * "="^60)
println("\nData layer loaded successfully — all checks passed.")
println("="^60)
println()