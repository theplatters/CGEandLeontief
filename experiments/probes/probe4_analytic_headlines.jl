# probe4_analytic_headlines.jl — how much of the matrix is arithmetic?
#
# Run from the repository root:
#   julia --project=. experiments/probes/probe4_analytic_headlines.jl
#
# Everything the paper reports as a "result" is checked here against a
# CLOSED FORM that needs no solve:
#
#   1. The reported real-GDP metric is the household-consumption Tornqvist
#      index (src/core/equilibrium.jl: `real_gdp_index`). Under a demand-only
#      shock the price block is p = 1 and w = 1, so
#         F2:  real GDP = (1 - tau0 - G0) / (1 - tau0)   (a pure tax effect),
#         F3:  real GDP = 1                              (households untouched),
#         F1:  real GDP = Tornqvist(c0 -> c1) with c1 from the preference tilt
#                         alone (a second-order composition term).
#      Predicted vs the stored manifests of the mobile rows.
#   2. GAMMA / DELTA vs the analytic Leontief multiplier: with p = 1 the CES
#      production block collapses to fixed coefficients, so the sticky-wage
#      rows ARE the Leontief multiplier of `leontief_multiplier` -- the
#      GAMMA = DELTA equivalence is not a coincidence of two solves, it is one
#      solve and its analytic limit.

using BeyondHulten
using CSV, DataFrames, LinearAlgebra, TOML, Printf

const ROOT = joinpath(@__DIR__, "..", "..") |> normpath
design = TOML.parsefile(joinpath(ROOT, "experiments", "designs", "matrix_5x3_v3.toml"))
dat, prog = design["data"], design["programme"]

drops = Vector{Int}(dat["drops"])
data_full = read_data(String(dat["source_table"]) |> f -> last(splitpath(f)); datadir = ROOT)
data = recalibrate_open(retained_dataset(data_full, drops); exo_scale = 1.0)
N = length(data.factor_share)

imp = CSV.read(joinpath(ROOT, prog["source"]), DataFrame)
rows = imp[imp.year .== prog["year"], :]
c1r, c2r = Int.(prog["column_slice"])
v = Float64.(Matrix(rows[1:1, c1r:c2r])[:])[1:N]
ψ = v ./ sum(v)
g = (Float64(prog["total_eur_m"]) / data.gdp_production) .* ψ
G0 = sum(g)
τ0 = sum(data.gov_demand)
s = data.saving_rate
cs = data.consumption_share
c0 = data.household_baseline
σ = 0.9

println("G0 = ", round(G0; digits = 8), "  tau0 = ", round(τ0; digits = 8),
    "  1 - tau0 = ", round(1 - τ0; digits = 8), "  s = ", round(s; digits = 8))

# ── 1. closed-form mobile headlines ───────────────────────────────────────
E0 = (1 - τ0)                       # L = w = 1
pred_F2 = (1 - τ0 - G0) / (1 - τ0)
pred_F3 = 1.0

# F1 tilt, verbatim the run.jl rule; consumption = (1-s) E (cs .* d) / agg at p = 1
function tilt()
    pos = c0 .> 0
    ψ1 = ψ .* pos
    ψ1 ./= sum(ψ1)
    return 1.0 .+ G0 .* ψ1 ./ max.(c0, 1e-12)
end
d = tilt()
agg = sum(cs .* d)
c1_F1 = (1 - s) .* E0 .* (cs .* d) ./ agg
pred_F1 = tornqvist_quantity_index(ones(N), c1_F1, ones(N), c0)
@printf("\npredicted (closed form, no solve): F1 = %.12f  F2 = %.12f  F3 = %.12f\n",
    pred_F1, pred_F2, pred_F3)

println("\nstored manifests (matrix_5x3-v3):")
for L in ("BF", "ALPHA", "BETA"), F in ("F1", "F2", "F3")
    id = "matrix_5x3-v3-" * L * "-" * F
    man = TOML.parsefile(joinpath(ROOT, "runs", id, "manifest.toml"))
    m = man["metrics"]
    pred = F == "F1" ? pred_F1 : F == "F2" ? pred_F2 : pred_F3
    @printf("  %-20s realGDP = %.12f  predicted = %.12f  diff = %+.3e   L = %.12f\n",
        L * "-" * F, Float64(m["real_gdp"]), pred, Float64(m["real_gdp"]) - pred,
        Float64(m["employment"]))
end

# Closed form of the F1 index: with c1/c0 = d/agg and shares s1 = s0 * d/agg,
#   ln Q = 0.5 [ sum_i cs_i ln(d_i/agg) + sum_i cs_i (d_i/agg) ln(d_i/agg) ]
#        = 0.5 [ E(ln d) + E_d(ln d) ] - ln(agg),   E over the consumption
#   shares cs (the baseline is proportional to cs at p = 1).
# The FIRST-ORDER term vanishes identically because the tilt is budget
# neutral:  sum_i cs_i (d_i/agg - 1) = 0.
a = d ./ agg
active = cs .> 0
lnQ_closed = 0.5 * (sum(cs[active] .* log.(a[active])) + sum(cs[active] .* a[active] .* log.(a[active])))
@printf("\nF1 closed form: ln Q = %+.10e (from the tilt alone) vs Tornqvist %+.10e (diff %.2e)\n",
    lnQ_closed, log(pred_F1), lnQ_closed - log(pred_F1))
@printf("F1 budget neutrality: sum_i cs_i (d_i/agg - 1) = %+.3e   (the first-order term)\n",
    sum(cs .* (a .- 1)))
@printf("F1 as a fraction of the impulse: index change = %.4f %% of G0   (F2: %.1f %%, F3: 0 %%)\n",
    100 * (pred_F1 - 1) / G0, 100 * (pred_F2 - 1) / G0)
@printf("F1 tilt dispersion: max d_i = %.3f, sectors with d_i > 1.01: %d of %d, their consumption mass = %.4f\n",
    maximum(d), count(>(1.01), d), length(d), sum(cs[d .> 1.01]))

# ── 2. GAMMA / DELTA vs the analytic Leontief multiplier ──────────────────
println("\n── GAMMA / DELTA vs the analytic Leontief multiplier (no solve) ──")
for (F, mode) in (("F3", :F3), ("F2", :F2))
    ana = leontief_multiplier(data, g; mode = mode)
    for L in ("GAMMA", "DELTA")
        id = "matrix_5x3-v3-" * L * "-" * F
        sol = CSV.read(joinpath(ROOT, "runs", id, "solution.csv"), DataFrame)
        q = Float64.(sol.quantity)
        @printf("  %-12s max|y - y_leontief| = %.3e   rel = %.3e   L_num = %.10f  L_ana = %.10f\n",
            L * "-" * F, maximum(abs.(q .- ana.y)), maximum(abs.(q .- ana.y)) / maximum(abs.(ana.y)),
            dot(data.factor_share, q), ana.L)
    end
end
