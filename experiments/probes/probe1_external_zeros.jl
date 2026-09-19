# probe1_external_zeros.jl — why is the net external position ~ 0 in the
# GAMMA / DELTA rows?
#
# Run from the repository root:
#   julia --project=. experiments/probes/probe1_external_zeros.jl
#
# For every executed matrix_5x3_v3 cell this rebuilds the calibrated model,
# loads the STORED solution (runs/<cell>/solution.csv + manifest wage) and
# prints, at that point:
#   * max|residual| of the cell's own system,
#   * |p . market_clearing_residuals|  (the value of aggregate excess demand),
#   * the canary  S - (I+X-M) + T      (external_balance_canary(...).diff),
#   * their difference  (the identity gap),
#   * the F3 decomposition: the canary with the M_prog = -F term switched off,
#   * the gross programme values F = p.g and its import content.
#
# Read-only: no model is solved, no run is created, no src/ file is touched.

using BeyondHulten
using CSV, DataFrames, LinearAlgebra, TOML, Printf

const ROOT = joinpath(@__DIR__, "..", "..") |> normpath
const DESIGN = "matrix_5x3_v3"
const CELLS = ["BF", "ALPHA", "BETA", "GAMMA", "DELTA"]
const FINS = ["F1", "F2", "F3"]

design = TOML.parsefile(joinpath(ROOT, "experiments", "designs", DESIGN * ".toml"))
dat = design["data"]
prog = design["programme"]

# ── calibration (mirrors experiments/run.jl build_reference) ──────────────
drops = Vector{Int}(dat["drops"])
data_full = read_data(String(dat["source_table"]) |> f -> last(splitpath(f)); datadir = ROOT)
data_v1 = retained_dataset(data_full, drops)
N = length(data_v1.factor_share)
# the reference continuation ends at exo_scale = 1
data = recalibrate_open(data_v1; exo_scale = 1.0)
shocks = Shocks(ones(N), ones(N), zeros(N))

# ── programme vectors (mirrors programme_vectors) ─────────────────────────
imp = CSV.read(joinpath(ROOT, prog["source"]), DataFrame)
rows = imp[imp.year .== prog["year"], :]
c1, c2 = Int.(prog["column_slice"])
v = Float64.(Matrix(rows[1:1, c1:c2])[:])
v = v[1:N]
ψ = v ./ sum(v)
g = (Float64(prog["total_eur_m"]) / data.gdp_production) .* ψ
G0 = sum(g)

"F1 tilt, verbatim copy of run.jl's `tilt_g0_over_c0` rule."
function f1_tilt(baseline, ψ, g)
    pos = baseline .> 0
    ψ1 = ψ .* pos
    ψ1 = ψ1 ./ sum(ψ1)
    return 1.0 .+ sum(g) .* ψ1 ./ max.(baseline, 1e-12)
end

financing(id) = id == "F1" ? PreferenceReallocation(f1_tilt(data.household_baseline, ψ, g)) :
                id == "F2" ? TaxFinanced(g) :
                id == "F3" ? ExternalDebt(g) :
                error("unknown financing $id")

function cell_model(cell)
    fin = financing(cell["financing"])
    θ, ϵ, σ, η = Float64(cell["theta"]), Float64(cell["epsilon"]),
                 Float64(cell["sigma"]), Float64(cell["eta"])
    lab = cell["labor"]
    if lab in ("BF", "ALPHA")
        return mobile_labor_model(data, shocks, θ, ϵ, σ, η; financing = fin)
    elseif lab == "BETA"
        return mobile_labor_model(data, shocks, θ, ϵ, σ, η; financing = fin,
            eta_s = Float64(cell["eta_s"]))
    elseif lab == "GAMMA"
        return mobile_labor_model(data, shocks, θ, ϵ, σ, η; closure = :fixed, financing = fin)
    elseif lab == "DELTA"
        return delta_model(data, shocks; ε = Float64(cell["delta_epsilon"]), financing = fin)
    end
    error("unknown labor $lab")
end

# canary with the F3 programme-inflow term M_prog = -dot(p, additive) removed
function canary_no_mprog(model, X)
    c = external_balance_canary(model, X)
    blocks = BeyondHulten._mobile_market_demand(model, X[1:N], X[N+1:2N], X[2N+1])
    F = dot(X[1:N], blocks.additive)
    return model.financing isa ExternalDebt ? c.diff + F : c.diff
end

println("calibration: N = ", N, ", s = ", round(data.saving_rate; digits = 6),
    ", tau0 = ", round(sum(data.gov_demand); digits = 6),
    ", G0 = ", round(G0; digits = 6),
    ", c0 gross = ", round(sum(data.household_baseline); digits = 6))
println()

fmt(x) = @sprintf("%+.6e", x)
hdr = rpad("cell", 16) * rpad("|resid|", 12) * rpad("p.mktr", 12) *
      rpad("canary", 12) * rpad("gap", 11) * rpad("canary+F", 12) *
      rpad("F=p.g", 11) * rpad("imp.prog", 11) * "max|dp|"
println(hdr); println("-"^length(hdr))

rowsout = DataFrame(cell = String[], resid = Float64[], pmkt = Float64[],
    canary = Float64[], gap = Float64[], canary_no_prog = Float64[],
    F = Float64[], import_prog = Float64[], maxdp = Float64[])

for L in CELLS, F in FINS
    id = "matrix_5x3-v3-" * L * "-" * F
    rundir = joinpath(ROOT, "runs", id)
    isdir(rundir) || continue
    cell = design["cells"][id]
    sol_df = CSV.read(joinpath(rundir, "solution.csv"), DataFrame)
    man = TOML.parsefile(joinpath(rundir, "manifest.toml"))
    p = Float64.(sol_df.price)
    q = Float64.(sol_df.quantity)
    w = Float64(man["metrics"]["wage"])
    model = cell_model(cell)
    fixed = labor_closure(model.options) isa FixedWageClosure
    Xsys = fixed ? [p; q] : [p; q; w]
    Xcan = fixed ? [p; q; 1.0] : [p; q; w]
    resid = maximum(abs, equilibrium_residuals(model, Xsys))
    mktr = market_clearing_residuals(model, Xcan)
    pm = dot(p, mktr)
    c = external_balance_canary(model, Xcan)
    cn = canary_no_mprog(model, Xcan)
    Fv = dot(p, g)
    imp = external_balance(model.financing, model, p)
    maxdp = maximum(abs.(p .- 1))
    println(rpad(id, 16) * rpad(@sprintf("%.3e", resid), 12) * rpad(fmt(pm), 12) *
            rpad(fmt(c.diff), 12) * rpad(fmt(pm - c.diff), 11) * rpad(fmt(cn), 12) *
            rpad(@sprintf("%.6f", Fv), 11) * rpad(@sprintf("%.6f", imp), 11) *
            @sprintf("%.3e", maxdp))
    push!(rowsout, (id, resid, pm, c.diff, pm - c.diff, cn, Fv, imp, maxdp))
end

CSV.write(joinpath(ROOT, "experiments", "probes", "probe1_external_zeros.csv"), rowsout)
println()
println("wrote experiments/probes/probe1_external_zeros.csv")

# ── the tracking test: does the fixed-row canary follow the solver's stop? ─
println()
println("tracking test (fixed rows): canary / max|residual|")
for L in ("GAMMA", "DELTA"), F in FINS
    r = rowsout[findfirst(==( "matrix_5x3-v3-" * L * "-" * F), rowsout.cell), :]
    println("  ", rpad(r.cell, 18), "canary = ", fmt(r.canary),
        " |resid| = ", @sprintf("%.3e", r.resid),
        " ratio = ", @sprintf("%+.3f", r.canary / max(r.resid, 1e-300)))
end

# ── the same test across generations: v2 stops looser than v3 (ADR-0015) ──
println()
println("generation comparison (v2 vs v3): canary and residual from the manifests")
for L in ("GAMMA", "DELTA", "ALPHA"), F in FINS
    for gen in ("v2", "v3")
        id = "matrix_5x3-" * gen * "-" * L * "-" * F
        f = joinpath(ROOT, "runs", id, "manifest.toml")
        isfile(f) || continue
        m = TOML.parsefile(f)
        println("  ", rpad(id, 20), "canary = ", fmt(Float64(m["diagnostics"]["canary_diff"])),
            " |resid| = ", @sprintf("%.3e", Float64(m["gates"]["residual"]["value"])),
            " max|dp| = ", @sprintf("%.3e", Float64(m["metrics"]["max_abs_price_dev"])))
    end
end
