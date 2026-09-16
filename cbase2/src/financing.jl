# ═══════════════════════════════════════════════════════════════════════════════
# cbase2 financing closures (Foundation II, Stage 1.1) — src/financing.jl
# ═══════════════════════════════════════════════════════════════════════════════
#
# This file is cbase2-only (no parent counterpart). It defines the three
# financing closures selected in `docs/DOCS_ASSESSMENT.md` (Foundation II)
# and the demand-side hooks consumed by the equilibrium systems in
# `src/core/mobile_labor.jl`. The supertype `AbstractFinancing` is defined
# in `src/core/interface.jl` (surgical core edit, recorded DIFF — see
# process_comments.md, chapter "Notebook 03").
#
# Design (agreed 2026-09-15, recorded in process_comments.md):
#
#   Programme: real bundle g_i = m · G0 · ψ_i  (numeraire units), with ψ the
#   sectoral incidence (year-1 impulse shares) and G0 = 40,300 EUR m scaled
#   into model units by the notebook (1 model income unit = GDP at basic
#   prices ≈ 2,864,724 EUR m). Same ψ and scale for F1/F2/F3.
#
#   F1  PreferenceReallocation(d): compositional household demand shift.
#       Preference weights become demand_shock .* d; the household CES
#       normalizer keeps Σ p_i c_i = E exactly (budget closed by
#       construction). NO additive demand; at baseline prices the shift
#       redirects m·G0·ψ_i toward programme sectors and shrinks all other
#       categories pro rata within fixed E.
#
#   F2  TaxFinanced(g): additive real public/investment demand g_i;
#       endogenous public budget T(p) = Σ p_i g_i levied as a lump-sum tax:
#       household expenditure E = w·ΣL − T. Budget closes by construction
#       (Walras: Σ p(y − int − c − g) = wΣL − (wΣL − T) − T = 0).
#
#   F3  ExternalDebt(g): identical additive demand; household untaxed; the
#       external balance F = Σ p_i g_i is recorded post-solve (accounting
#       open-economy closure, import leakage mechanical via Ω_raw).
#
#   The retired unfinanced autonomous/investment shocks are NOT part of this
#   interface; Shocks keeps its legacy fields (parent compatibility) but
#   cbase2 experiments always pass zero vectors.
# ═══════════════════════════════════════════════════════════════════════════════

"""Baseline reference: no programme, household untaxed, no additive demand."""
struct NoFinancing <: AbstractFinancing end

"""
	PreferenceReallocation(shift)

F1 — compositional household demand shift. `shift` is a vector of strictly
positive preference shifters d_i (length N); the equilibrium budget closes
by construction through the household CES normalizer.
"""
struct PreferenceReallocation <: AbstractFinancing
	shift::Vector{Float64}
	function PreferenceReallocation(shift::Vector{<:Real})
		n = length(shift)
		n > 0 || throw(ArgumentError("preference shift vector must be non-empty"))
		all(>(0), shift) ||
			throw(ArgumentError("preference shifters must be strictly positive (a zero weight cannot be redirected)"))
		new(Float64.(shift))
	end
end

"""
	TaxFinanced(g)

F2 — additive real public/investment demand `g_i` (numeraire units),
financed by an endogenous lump-sum tax T(p) = Σ p_i g_i on household income.
"""
struct TaxFinanced <: AbstractFinancing
	g::Vector{Float64}
	function TaxFinanced(g::Vector{<:Real})
		n = length(g)
		n > 0 || throw(ArgumentError("public demand vector must be non-empty"))
		all(≥(0), g) ||
			throw(ArgumentError("public demand must be non-negative (use ExternalDebt/PreferenceReallocation for other designs)"))
		new(Float64.(g))
	end
end

"""
	ExternalDebt(g)

F3 — additive real public/investment demand `g_i` financed externally: the
household is untaxed and the external balance F = Σ p_i g_i records the
resource inflow post-solve.
"""
struct ExternalDebt <: AbstractFinancing
	g::Vector{Float64}
	function ExternalDebt(g::Vector{<:Real})
		n = length(g)
		n > 0 || throw(ArgumentError("public demand vector must be non-empty"))
		all(≥(0), g) ||
			throw(ArgumentError("public demand must be non-negative"))
		new(Float64.(g))
	end
end

# ── Demand-side hooks (called by problem() / problem_fixed() / solve()) ──

"""
	preference_weights(fin, demand_shock)

Effective household preference weights. F1 multiplies its shift into the
shock weights; all other closures pass the weights through unchanged.
"""
preference_weights(::AbstractFinancing, demand_shock::Vector{Float64}) = demand_shock
preference_weights(f::PreferenceReallocation, demand_shock::Vector{Float64}) =
	demand_shock .* f.shift

"""
	household_expenditure(fin, wage_income, p)

Household expenditure base E: wage income minus the lump-sum tax
T(p) = Σ p_i g_i under F2; unchanged otherwise (F1 composition shifts
operate within E; F3 is externally financed).
"""
household_expenditure(::AbstractFinancing, wage_income::Real, p) = wage_income
household_expenditure(f::TaxFinanced, wage_income::Real, p) =
	wage_income - dot(f.g, p)

"""
	additive_demand(fin, N)

Additive real public/investment demand entering market clearing: the bundle
g under F2/F3, zeros otherwise (F1 shifts composition only).
"""
additive_demand(::AbstractFinancing, N::Int) = zeros(N)
additive_demand(f::Union{TaxFinanced, ExternalDebt}, N::Int) = f.g

"""
	public_budget(fin, p)

Domestic public budget T(p) = Σ p_i g_i: positive under F2 (levied as a
lump-sum tax) and F3 (funded externally); zero for F1 and the baseline.
"""
public_budget(::AbstractFinancing, p) = 0.0
public_budget(f::Union{TaxFinanced, ExternalDebt}, p) = dot(f.g, p)

"""
	external_balance(fin, p)

External balance F = Σ p_i g_i: positive only under F3 (externally financed
resources entering via net imports); zero otherwise.
"""
external_balance(::AbstractFinancing, p) = 0.0
external_balance(f::ExternalDebt, p) = dot(f.g, p)

"""
	has_additive_anchor(fin)

Whether the financing carries an additive demand anchor that breaks the
homogeneity of the fixed-wage η ≈ 1 system (only F2/F3 do; F1 is purely
compositional and does NOT anchor scale — see the scale-indeterminacy
guard in `_solve_fixed`).
"""
has_additive_anchor(::AbstractFinancing) = false
has_additive_anchor(::Union{TaxFinanced, ExternalDebt}) = true

export AbstractFinancing, NoFinancing, PreferenceReallocation, TaxFinanced,
	ExternalDebt, preference_weights, household_expenditure, additive_demand,
	public_budget, external_balance, has_additive_anchor
