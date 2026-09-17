# ═══════════════════════════════════════════════════════════════════════════════
# verify_retained_pipeline.jl — validation of the rebuilt sector pipeline
# (cbase2/src/calibration.jl, 2026-09-17 repair of review findings 2 + 3).
#
# Checks, for every DATASET_VARIANTS configuration:
#   1. Ω_raw rows sum to 1 (the CES/CD price-index contract);
#   2. Σ labor_share = 1 (baseline factor income at w = 1 is the income unit);
#   3. the intermediate price index at p = 1 equals 1 for every θ, including
#      the 0.999 / 1.0 / 1.0001 neighbourhood where the old pipeline showed the
#      artificial CES-vs-CD jump (2.2e-13 vs 1.0);
#   4. p = w = 1 is a zero-profit equilibrium (max|p - cost| ≈ machine zero)
#      for θ below, at and above the CD branch;
#   5. the open-economy calibration is internally consistent
#      (E_h0 = 1 - ΣgG, saving rate and clamp mass reported).
# Run from anywhere:  julia --project=. cbase2/scripts/verify_retained_pipeline.jl
# ═══════════════════════════════════════════════════════════════════════════════

using CSV, DataFrames, LinearAlgebra, Printf
const CB = normpath(joinpath(@__DIR__, ".."))
for f in ["interface.jl", "solution.jl", "ces.jl", "mobile_labor.jl", "leontief.jl", "util.jl"]
	include(joinpath(CB, "src", "core", f))
end
include(joinpath(CB, "src", "calibration.jl"))

const THETAS = (0.5, 0.9, 0.999, 1.0, 1.0001, 1.1, 2.0)

function main()
	data_full = read_data(joinpath(CB, "data_raw", "I-O_DE2019_formatiert.csv"))
	println("full table: N = ", length(data_full.λ),
		", Σ labor_share = ", sum(data_full.labor_share))

	failures = String[]
	for (name, drops) in [("full", DATASET_VARIANTS["full"]),
						  ("70s", DATASET_VARIANTS["70s"]),
						  ("reduced", DATASET_VARIANTS["reduced"])]
		d = retained_dataset(data_full, drops)
		de = recalibrate_open(d; exo_scale = 1.0)
		n = length(d.λ)
		cover = dataset_coverage(data_full, drops)

		# (1) Ω_raw row normalization; old slice kept for contrast
		rs = vec(sum(d.Ω_raw; dims = 2))
		keep = setdiff(1:length(data_full.λ), drops)
		old_rs = vec(sum(data_full.Ω_raw[keep, keep]; dims = 2))

		# (3) price index at p = 1, and (4) zero-profit at p = w = A = 1
		ip_dev = maximum(maximum(abs.(_intermediate_price(d.Ω_raw, ones(n), th) .- 1)) for th in THETAS)
		zp = Dict{eltype(THETAS),Float64}()
		for th in THETAS
			ip = _intermediate_price(d.Ω_raw, ones(n), th)
			cost = _ces_unit_cost(ones(n), min.(d.factor_share, 1 - 1e-12), ones(n), ip, 0.5)
			zp[th] = maximum(abs.(1 .- cost))
		end

		# (5) open-economy consistency
		E_h0 = 1.0 - sum(de.gov_demand)
		Mλ = d.Ω_raw' * ((1 .- d.factor_share) .* d.λ)
		c0_dom = d.λ .- Mλ .- (1 .- de.import_margin) .* (de.gov_demand .+ de.exo_demand) .- de.exports_demand
		clamp_mass = sum(abs.(min.(c0_dom, 0)))

		# (6) the rebuilt table preserves the full-table domestic/import split
		# of each final-demand category (the FD import/tax rows are scaled by the
		# retained category share in retained_io_table).
		fd_full = final_demand_split(data_full.io, length(data_full.λ))
		fd_red = final_demand_split(d.io, n)
		den_full = max.(vec(sum(fd_full.tot; dims = 1)), 1e-12)
		den_red = max.(vec(sum(fd_red.tot; dims = 1)), 1e-12)
		domfrac_full = vec(sum(fd_full.dom; dims = 1)) ./ den_full
		domfrac_red = vec(sum(fd_red.dom; dims = 1)) ./ den_red
		domfrac_dev = maximum(abs.(domfrac_full .- domfrac_red))

		println("\n══ ", name, " (N=", n, ", drops=", drops, ") ══")
		println("  coverage: GO ", cover.gross_share_kept, "%, VA ", cover.va_share_kept,
			"%, FD ", cover.fd_share_kept, "%")
		@printf("  Ω_raw row sums: rebuilt [%.12f, %.12f]   old slice [%.6f, %.6f]\n",
			minimum(rs), maximum(rs), minimum(old_rs), maximum(old_rs))
		@printf("  Σ labor_share = %.12f   Σ λ·fs = %.12f\n",
			sum(d.labor_share), sum(d.λ .* d.factor_share))
		@printf("  max|ip(p=1) - 1| over θ ∈ {%s} : %.3e\n",
			join(THETAS, ", "), ip_dev)
		@printf("  max|p - cost| at p=w=1: θ=0.5 %.1e, θ=1.0 %.1e, θ=1.1 %.1e, θ=2.0 %.1e\n",
			zp[0.5], zp[1.0], zp[1.1], zp[2.0])
		@printf("  E_h0 = 1 - ΣgG = %.6f   saving rate s = %.6f   clamp mass = %.6f\n",
			E_h0, de.saving_rate, clamp_mass)
		@printf("  FD domestic-share drift vs full table (max over categories) = %.3e\n",
			domfrac_dev)

		minimum(rs) >= 1 - 1e-10 && maximum(rs) <= 1 + 1e-10 ||
			push!(failures, "$name: Ω_raw rows deviate from 1")
		isapprox(sum(d.labor_share), 1.0; atol = 1e-10) ||
			push!(failures, "$name: Σ labor_share ≠ 1")
		ip_dev < 1e-10 || push!(failures, "$name: intermediate price index ≠ 1 at p=1")
		maximum(values(zp)) < 1e-10 ||
			push!(failures, "$name: p=w=1 is not a zero-profit equilibrium")
		isapprox(sum(d.λ .* d.factor_share), 1.0; atol = 1e-10) ||
			push!(failures, "$name: baseline factor income Σλ·fs ≠ 1")
		domfrac_dev < 1e-10 ||
			push!(failures, "$name: FD domestic/import split drifted from the full table")
	end

	isempty(failures) || error("retained-pipeline validation FAILED:\n  " * join(failures, "\n  "))
	println("\nALL RETAINED-PIPELINE VALIDATION CHECKS PASSED")
end

main()
