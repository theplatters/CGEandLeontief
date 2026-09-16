# ═══════════════════════════════════════════════════════════════════════════════
# cbase2 v3 open-economy calibration (Notebook 03b) — src/calibration.jl
# ═══════════════════════════════════════════════════════════════════════════════
#
# National-accounts-consistent demand structure (fixes the v2 inconsistency):
#
#   Y = C + I + G + X - M      (expenditure side, domestic output valuation)
#   Household:Disposable income splits into consumption (1-s)E and saving sE.
#   Exogenous injections: government G (tax-financed, F2 balanced),
#   investment I and exports X (backed by saving: S = I + X - M at baseline).
#   Imports M: fixed sector margins m_i on consumption, government,
#   investment and the programme bundle. Exports carry NO margin (domestic
#   sales abroad).
#
# The saving rate is NOT assumed -- it is the data-implied residual
#   s = 1 - sum(c0_gross)/E_h0,  E_h0 = 1 - sum(gG),
# with c0 the baseline household block calibrated so the model clears
# exactly at y = lambda. The identity S = I + X - M then holds at baseline
# up to the clamping/accounting residual (validated in the notebooks).
#
# POST-KEYNESIAN ALTERNATIVE (major option, see process_comments.md): the
# single s lumps household and corporate saving (Kaldor-style split is the
# natural refinement: s_h on wage income, s_c on profits, two households).
# ═══════════════════════════════════════════════════════════════════════════════

"""
	recalibrate_open(data, cbroot) -> Data

v3 open-economy recalibration of a v1 `Data` object. Consumes the §4.1
artifacts (AC_domestic_final_demand.csv) plus the raw io table for the
gross category columns. Returns a new Data with:
  gov_demand      gross government consumption (margin applies)
  exo_demand      gross investment: equipment + construction + inventories
  exports_demand  exports (no margin)
  household_baseline  c0_gross: residual baseline household demand
  consumption_share   omega = c0_gross / sum(c0_gross)  (CPI weights, sum 1)
  import_margin       sector margin m_i of total final demand
  saving_rate         s = 1 - sum(c0_gross)/(1 - sum(gG))
Asserts the finiteness gate: round-gain column sums strictly below 1.
"""
function recalibrate_open(data::Data, cbroot::String; exo_scale::Real = 1.0)
	(; Ω_raw, factor_share, λ, gross_output_basic) = data
	n = length(factor_share)
	fs = factor_share
	0 <= exo_scale <= 1 || throw(ArgumentError("exo_scale must be in [0, 1]"))

	# ── §4.1 artifact: domestic final demand by sector and category ──
	fd = CSV.read(joinpath(cbroot, "data_processed", "AC_domestic_final_demand.csv"), DataFrame)
	@assert size(fd, 1) == n "AC_domestic_final_demand rows must match sectors"
	scale = 1.0 / data.gdp_production          # EUR m -> model units (GDP_P = 1)

	# ── Gross category columns from the raw io table (FD = 75:81) ──
	# 1 private_consumption, 2 private_orgs, 3 government, 4 equipment,
	# 5 construction, 6 inventories, 7 exports
	FD = 75:81
	tot = Matrix{Float64}(data.io[1:n, FD])
	dom = Matrix{Float64}(fd[1:n, 2:8])

	# Sector import margin of total final demand (proportional §4.1 split)
	m = 1.0 .- vec(sum(dom; dims=2)) ./ max.(vec(sum(tot; dims=2)), 1e-12)
	m = clamp.(m, 0.0, 1.0)

	gG = tot[:, 3] .* scale                        # government (gross, margin applies)
	inv = (tot[:, 4] .+ tot[:, 5] .+ tot[:, 6]) .* scale .* exo_scale  # investment (gross, margin applies)
	expo = tot[:, 7] .* scale .* exo_scale         # exports (NO margin -- domestic sales abroad)

	# ── Baseline household residual (gross), exact clearing at λ ──
	Mλ = Ω_raw' * ((1.0 .- fs) .* λ)
	c0_dom = λ .- Mλ .- (1.0 .- m) .* (gG .+ inv) .- expo
	n_clamped = count(<(0), c0_dom)
	n_clamped > 0 && println("recalibrate_open: household residual clamped to 0 in ",
		n_clamped, " sectors (total clamped mass = ",
		round(sum(abs.(min.(c0_dom, 0))); digits=6), " model units)")
	c0_dom = max.(c0_dom, 0.0)
	c0_gross = c0_dom ./ max.(1.0 .- m, 1e-6)

	# ── Saving rate (data-implied) and CPI weights ──
	E_h0 = 1.0 - sum(gG)                # baseline after-tax income (baseline income = 1)
	s = 1.0 - sum(c0_gross) / E_h0
	@assert -1.0 < s < 1.0 "calibrated saving rate out of range: $s"
	ω = c0_gross ./ sum(c0_gross)

	# ── Finiteness gate: worst-case round-gain column sums (F2, marginal tau = 0) ──
	# Enforced only for s >= 0 (negative-s calibrations occur only inside the
	# exo_scale bisection and are never solved).
	colsums = (1.0 .- fs) .+ (1.0 .- m) .* (1.0 - s) .* fs
	s >= 0 && @assert all(<(1), colsums) "finiteness gate failed: max column sum = $(maximum(colsums))"

	println("recalibrate_open (v3): saving rate s = ", round(s; digits=4),
		", tau0 = ", round(sum(gG); digits=4),
		", export share = ", round(sum(expo); digits=4),
		", investment share = ", round(sum(inv); digits=4))

	return Data(data.io, data.Ω, data.Ω_raw, ω, data.factor_share, data.λ,
		data.labor_share, data.consumption_share_gross_output, data.grossy,
		data.value_added, data.gross_output_basic, data.value_added_components,
		data.imports_intermediate, data.import_share, data.domestic_final_demand,
		gG, c0_gross, m, inv, expo, s,
		data.gdp_production, data.gdp_income, data.gdp_expenditure)
end
