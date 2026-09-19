# probe19_gamma_dual_labour.jl -- GAMMA with a DUAL labour market: the shocked
# (capacity-pressured) sectors are inflexible in QUANTITY, the rest adjust.
#
# READ-ONLY probe: no run, no src/ change, nothing committed.
#   julia --project=. experiments/probes/probe19_gamma_dual_labour.jl
#
# The variant the fixed-wage row is missing: the wage is pinned (w = 1) AND a
# group of sectors cannot expand its labour -- the insiders. Employment in the
# insider sectors is the datum Lbar_i, employment elsewhere is the
# cost-minimizing demand, so the dualism is in QUANTITIES, not in wages.
#
# The mechanism follows from zero profit with the quantity constraint. In a
# sector that must produce y_i with L_i = Lbar_i, the firm can only substitute
# towards intermediates, and the CES production function (elasticity eps in
# labour and the intermediate aggregate, factor share fs) gives the constrained
# unit cost
#
#   y = A [ fs^(1/eps) L^((eps-1)/eps) + (1-fs)^(1/eps) Z^((eps-1)/eps) ]^(eps/(eps-1))
#   => Z(y, L) solved from it,  c_con = (w L + ip Z) / y  >=  c_uncon.
#
# So the insider price RISES to cover the cost of the constraint (and the sector
# has a hard output CEILING y_max = A fs^(eps/(eps-1)) Lbar whenever eps < 1,
# i.e. when labour and intermediates are complements -- our calibration is
# eps = 0.5), while the demand it cannot serve spills over to the outsiders,
# whose employment expands at the unchanged wage. That is the dual-labour-market
# story in one system: insider prices up, outsider employment up, and the wage
# identical in both segments.
#
# Reported: aggregate max|p-1| and the deflator; the price response SPLIT into
# the insider and outsider groups (the dualism signature); total employment and
# the outsider-only employment change; consumption; and the insider output
# relative to its ceiling.

using BeyondHulten
using CSV, DataFrames, LinearAlgebra, Printf
import NonlinearSolve

const BH = BeyondHulten
const ROOT = joinpath(@__DIR__, "..", "..") |> normpath

data = recalibrate_open(retained_dataset(read_data("I-O_DE2019_formatiert.csv"; datadir = ROOT),
    Vector{Int}([])); exo_scale = 1.0)
N = length(data.factor_share)
fs = data.factor_share
ϵ = 0.5; θ = 0.5; σ = 0.9

imp = CSV.read(joinpath(ROOT, "cbase2/data_raw/impulses.csv"), DataFrame)
rows = imp[imp.year .== 2024, :]
ψ = Matrix{Float64}(rows[1:1, 3:73])[:][collect(1:N)]
ψ = ψ ./ sum(ψ)
g = (40300.0 / data.gdp_production) .* ψ

shocks = Shocks(ones(N), ones(N), zeros(N))
ref_sol = solve(mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 1.0);
    init = [ones(N); data.λ; 1.0; 0.0])
ref_cons = real_consumption(ref_sol)

function f1_tilt(baseline, ψ, g)
    pos = baseline .> 0
    ψ1 = ψ .* pos; ψ1 = ψ1 ./ sum(ψ1)
    return 1.0 .+ sum(g) .* ψ1 ./ max.(baseline, 1e-12)
end
fin_of(id) = id == "F1" ? PreferenceReallocation(f1_tilt(data.household_baseline, ψ, g)) :
             id == "F2" ? TaxFinanced(g) : ExternalDebt(g)

"""Constrained unit cost: labour fixed at `L`, output `y` required."""
function constrained_cost(A, fs, w, ip, ϵ, L, y)
    e = (ϵ - 1) / ϵ
    rhs = (y ./ A) .^ e
    bracket = rhs .- (fs .^ (1 / ϵ)) .* (L .^ e)
    Z = BH._positive_floor(bracket ./ ((1 .- fs) .^ (1 / ϵ))) .^ (1 / e)
    return (w .* L .+ ip .* Z) ./ y
end

"""Output ceiling with L fixed: y_max = A alpha^(1/(eps-1)) L (finite when eps < 1).

Derivation: as Z -> inf with eps < 1, Z^((eps-1)/eps) -> 0, so
y -> A [fs^(1/eps) L^((eps-1)/eps)]^(eps/(eps-1)) = A fs^(1/(eps-1)) L.
The baseline sits at the fraction fs of this ceiling (L = fs y at the
calibration), so the diagnostic reads how far a constrained sector has been
pushed toward its capacity limit.
"""
ceiling(A, fs, ϵ, L) = A .* (fs .^ (1 / (ϵ - 1))) .* L

"""The dual-labour-market residual: 2N unknowns [p; y], w = 1, F = 0."""
function dual_residual!(out, X, pars)
    model, mask, wbar = pars
    (; data, options, shocks) = model
    (; θ, ϵ, σ) = options.elasticities
    p = BH._positive_floor(X[1:N]); y = BH._positive_floor(X[N+1:2N])
    ip = BH._intermediate_price(data.Ω_raw, p, θ)
    Lcm = BH._cost_minimizing_labor(p, y, wbar, model)
    L = ifelse.(mask, data.labor_share, Lcm)          # insiders frozen, outsiders free
    fin = model.financing
    ds_eff = BH.preference_weights(fin, shocks.demand_shock)
    E = household_expenditure(fin, model, BH._wage_bill(wbar, L), p, sum(L);
        external_transfer = 0.0)
    agg = sum(data.consumption_share .* ds_eff .* p .^ (1 - σ))
    c_dom = (1 .- data.saving_rate) .* (1 .- data.import_margin) .*
            (data.consumption_share .* ds_eff) .* E .* p .^ (-σ) ./ agg
    tfd = c_dom .+ (1 .- data.import_margin) .* BH.additive_demand(fin, N) .+
          (1 .- data.import_margin) .* (data.gov_demand .+ data.exo_demand) .+
          data.exports_demand
    intd = p .^ (-θ) .* (data.Ω_raw' * (p .^ ϵ .* shocks.supply_shock .^ (ϵ - 1) .*
        ip .^ (θ - ϵ) .* (data.A_bill ./ data.λ) .* y))
    c_uncon = BH._ces_unit_cost(shocks.supply_shock, fs, wbar, ip, ϵ)
    c_con = constrained_cost(shocks.supply_shock, fs, wbar, ip, ϵ, data.labor_share, y)
    out[1:N] .= p .- ifelse.(mask, c_con, c_uncon)
    out[N+1:2N] .= y .- intd .- tfd
    nothing
end

rmax_of(f, x, pars) = maximum(abs, (o = similar(x); f(o, x, pars); o))

function solve_ladder(f, init, pars; tol = 1e-12)
    res = NonlinearSolve.solve(NonlinearSolve.NonlinearProblem(f, copy(init), pars),
        reltol = 1e-9, abstol = 1e-9, maxiters = 60000)
    x = res.u; rmax = rmax_of(f, x, pars)
    for _ in 1:10
        rmax <= tol && break
        res = NonlinearSolve.solve(NonlinearSolve.NonlinearProblem(f, x, pars),
            NonlinearSolve.LevenbergMarquardt(); reltol = 1e-13, abstol = 1e-13, maxiters = 60000)
        r2 = rmax_of(f, res.u, pars); r2 < rmax || break
        x, rmax = res.u, r2
    end
    x, rmax
end

println("N = ", N, "  programme sectors = ", count(>(0.0), ψ), "  eps = ", ϵ)

function run_case(tag, mask, fin_id)
    fin = fin_of(fin_id)
    mo = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 1.0; closure = :fixed, financing = fin)
    # warm start from the unconstrained fixed-wage solution
    gsol = solve(mo; init = [ones(N); data.λ])
    x, rmax = solve_ladder(dual_residual!, [gsol.prices_raw; gsol.quantities], (mo, mask, 1.0))
    p = BH._positive_floor(x[1:N]); y = BH._positive_floor(x[N+1:2N])
    ip = BH._intermediate_price(data.Ω_raw, p, θ)
    Lcm = BH._cost_minimizing_labor(p, y, 1.0, mo)
    L = ifelse.(mask, data.labor_share, Lcm)
    cei = ceiling.(1.0, fs, ϵ, data.labor_share)
    # the kernel's welfare index (src/core/equilibrium.jl, fixed-wage branch)
    fin2 = mo.financing
    ds_eff2 = BH.preference_weights(fin2, shocks.demand_shock)
    E2 = household_expenditure(fin2, mo, BH._wage_bill(1.0, L), p, sum(L);
        external_transfer = 0.0)
    agg2 = sum(data.consumption_share .* ds_eff2 .* p .^ (1 - σ))
    cons_vec = (1 .- data.saving_rate) .* (data.consumption_share .* ds_eff2 .* E2 .* p .^ (-σ)) ./ agg2
    cons_idx = tornqvist_quantity_index(p, cons_vec, ones(N), data.household_baseline)
    cons = cons_idx / real_consumption(ref_sol)
    yin = isempty(findall(mask)) ? 0.0 : maximum((y ./ cei)[mask])
    @printf("%-22s %-3s |r|=%.1e  max|p-1|=%.6f  defl=%.6f  L=%.6f  dL_out=%+.5f  cons=%+.5f  y_in/y_ceil=%.4f\n",
        tag, fin_id, rmax, maximum(abs, p .- 1),
        sum(data.consumption_share .* p .^ (1 - σ))^(1 / (1 - σ)),
        sum(L), sum(L[.!mask]) - sum(data.labor_share[.!mask]),
        cons - 1, yin)
    pin = isempty(findall(mask)) ? 0.0 : maximum(abs, (p .- 1)[mask])
    pout = isempty(findall(.!mask)) ? 0.0 : maximum(abs, (p .- 1)[.!mask])
    @printf("%-22s %-3s    insider max|p-1|=%.6f   outsider max|p-1|=%.6f\n",
        "", "", pin, pout)
    return nothing
end

mask_prog = ψ .> 0
rank_emp = sortperm(data.labor_share; rev = true)
mask_half = falses(N); mask_half[rank_emp[1:cld(N, 2)]] .= true

println("\n── dual labour market: insiders = the shocked (capacity-pressured) sectors ──")
for fin_id in ("F1", "F2", "F3")
    run_case("dual (shocked rigid)", mask_prog, fin_id)
end
println("\n── dual labour market: insiders = the largest half by employment ──")
for fin_id in ("F2",)
    run_case("dual (largest half)", mask_half, fin_id)
end
println("\n── no dualism (all outsiders = the baseline GAMMA row), for contrast ──")
for fin_id in ("F2",)
    run_case("no dualism", falses(N), fin_id)
end
println("\ndone.")