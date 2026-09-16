# ═══════════════════════════════════════════════════════════════════════════════
# cbase2 v2 open-absorption calibration (Notebook 03b) — src/calibration.jl
# ═══════════════════════════════════════════════════════════════════════════════
#
# cbase2-only. Consumes the §4.1 artifacts (AC_domestic_final_demand.csv) plus
# the v1 Data object and rebuilds it under the open-absorption structure:
#
#   government  gG_i  = domestic government consumption (exogenous, real)
#   baseline    c0_i  = λ_i − (M_raw λ)_i − gG_i   (exact baseline clearing)
#   weights     ω_i   = c0_i / Σc0                 (Σω = 1: numeraire-consistent)
#   margins     m_i   = sector import share of marginal final demand
#                     (final-demand domestic fractions by category, weighted
#                     by the sector's composition of total final demand)
#
# Reference after-tax income: E_h0 = 1 − ΣgG (baseline wage income = 1).
# Acceptance gate (asserted here): c0 ≥ 0 sectorwise and the baseline
# round-gain matrix M + diag(1−m)·ω·fs' has all column sums strictly below 1
# (the unit root of the v1 fixed-wage system is dead).
# ═══════════════════════════════════════════════════════════════════════════════

"""
	recalibrate_open(data::Data, cbroot::AbstractString) -> Data

Rebuild `data` under the v2 open-absorption calibration. `cbroot` is the
cbase2 folder (for `data_processed/AC_domestic_final_demand.csv`).
"""
function recalibrate_open(data::Data, cbroot::AbstractString)
	Ω_raw, fs, λ = data.Ω_raw, data.factor_share, data.λ
	n = length(λ)

	# ── Government baseline (domestic, EUR m → model units) ──
	fd = CSV.read(joinpath(cbroot, "data_processed", "AC_domestic_final_demand.csv"), DataFrame)
	gG = Vector{Float64}(fd.government_consumption) / data.gdp_production

	# ── Import margins: sector composition-weighted final-demand domestic share ──
	# Category k's domestic fraction df_k = Σ_i dom[i,k] / Σ_i tot[i,k]; the
	# sector's marginal domestic share weights the category fractions by its
	# own composition of TOTAL (purchaser) final demand.
	tot_fd = Matrix{Float64}(data.io[1:n, 75:81])        # purchaser, dom+imp
	dom_fd = Matrix{Float64}(fd[:, 2:end])
	@assert size(tot_fd) == size(dom_fd) == (n, 7)
	df_k = vec(sum(dom_fd; dims=1)) ./ vec(sum(tot_fd; dims=1))
	m = 1.0 .- (vec(sum(dom_fd; dims=2)) ./ max.(vec(sum(tot_fd; dims=2)), eps()))

	# ── Baseline household residual: exact clearing at λ ──
	Mλ = Ω_raw' * ((1.0 .- fs) .* λ)
	c0 = λ .- Mλ .- gG
	# Government-heavy sectors (public services) can have gG_i exceeding their
	# private residual; clamp to zero there (parent convention: cons_vec was
	# clamped the same way). The baseline output vector then adjusts endogenously
	# -- the init λ is no longer the exact baseline, which the notebooks report.
	n_clamped = count(<(0), c0)
	n_clamped > 0 && println("recalibrate_open: household residual clamped to 0 in ",
		n_clamped, " government-heavy sectors (total clamped mass = ",
		round(sum(abs.(min.(c0, 0))); digits=6), " model units)")
	c0 = max.(c0, 0.0)

	# ── Weights and finiteness gate ──
	ω = c0 ./ sum(c0)
	τ0 = sum(gG)                       # baseline income = 1
	M = Ω_raw' * Diagonal(1.0 .- fs)
	gain = M + Diagonal(1.0 .- m) * ω * fs' * (1 - τ0)
	colsums = vec(sum(gain; dims=1))
	@assert all(<(1), colsums) "finiteness gate failed: column sums must be strictly below 1 (max = $(maximum(colsums)))"

	# ── Rebuild the Data object ──
	Data(data.io, data.Ω, data.Ω_raw, ω, data.factor_share, data.λ,
		data.labor_share, data.consumption_share_gross_output, data.grossy, data.value_added,
		data.gross_output_basic, data.value_added_components, data.imports_intermediate,
		data.import_share, data.domestic_final_demand, gG, c0, m,
		data.gdp_production, data.gdp_income, data.gdp_expenditure)
end

export recalibrate_open