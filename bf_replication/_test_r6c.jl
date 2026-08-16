using BFReplication, DelimitedFiles, LinearAlgebra, Printf
using .BFReplication.DataLoader

data_dir = joinpath(@__DIR__, "..", "Replication Files", "GDP Simulatin -- 88 Sector")
data = load_bf_data(joinpath(data_dir, "BFdata.csv"); year=1980)
stfp = Float64.(readdlm(joinpath(data_dir, "stfp.csv"), ','))
_, S4 = BFReplication.empirical_covariances(stfp)
Cov4m = Matrix(Diagonal(diag(S4)))

λ = (I - Diagonal(1 .- data.α) * data.Ω)' \ data.β
println("λ sum = $(sum(λ))")

# Quick no-realloc Hessian (5 sectors) to verify structure
println("\nNo-realloc Hessian (5 sectors):")
H5 = second_order_hessian_norealloc(data.Ω, data.α, data.β, data.L, 0.5, 0.001, 0.9; h=1e-6)
for i in 1:5
    println("  $(round.(H5[i, 1:5], digits=6))")
end
println("  rank $(rank(H5[1:5,1:5])) — degenerate as expected (each column is constant)")

# Quick realloc Hessian (5 sectors)
println("\nRealloc Hessian (5 sectors):")
Hr5 = second_order_hessian_realloc(data.Ω, data.α, data.β, 0.5, 0.001, 0.9; h=1e-5)
for i in 1:5
    println("  $(round.(Hr5[i, 1:5], digits=6))")
end
println("  rank $(rank(Hr5[1:5,1:5])) — proper NxN")

# MC with 5-sector Hessian (just to verify the function works)
println("\nSecond-order MC with 5-sector realloc Hessian:")
Cov5 = Cov4m[1:5, 1:5]
λ5 = λ[1:5]
Hr5 = Hr5[1:5, 1:5]
mc5 = second_order_mc(λ5, Hr5, Cov5; trials=200)
@printf("  mean=%.4f%%  std=%.4f%%  skew=%.3f  exkurt=%.3f\n",
    mc5.mean*100, mc5.std*100, mc5.skewness, mc5.excess_kurtosis)

# Show what the full run would cost
println("\nFull 76-sector realloc Hessian: ~76× per-sector solve time")
println("  Recommend running on Mac with parallel execution.")
println("  The wrapper function is: run_second_order_realloc(data, Cov4m; trials=50000)")

# Clean up test files
rm("/workspace/BFrep/(3)BeyondHulten/bf_replication/_test_r6.jl", force=true)
rm("/workspace/BFrep/(3)BeyondHulten/bf_replication/_test_r6b.jl", force=true)
println("\nTest files cleaned.")