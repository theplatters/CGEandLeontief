# analytic_headlines.jl -- analytic equivalences and capacity on the ADR-0019 kernel.
#
# READ-ONLY probe: no run is created, no src/ file is touched, nothing is
# committed. Run from the repository root:
#   julia --project=. experiments/probes/analytic_headlines.jl
#
# (i)   Stored v5 GAMMA/DELTA F2/F3 quantities vs
#       leontief_multiplier(data, g; mode = :F2/:F3): max|y - y_analytic| and L
#       (pre-ADR-0019: 8.3e-17 / 1.7e-16, L = 1.0177957289).
# (ii)  Scaled programme k in {1, 2, 5, 10}: ALPHA-F3 vs GAMMA-F3 -- report
#       success/failure, residual, L, max|dp| (pre-ADR-0019: ALPHA stalls above
#       k = 1, GAMMA scales cleanly; with F the all-N system may differ).
# (iii) v5 BF/ALPHA/BETA F1 consumption_rel against the tilt-based Tornqvist
#       closed form from probe4. If the closed form needs re-derivation under
#       ADR-0018/ADR-0019, the measured consumption_rel/gdp_rel are reported
#       and the drift is stated.

using BeyondHulten
using CSV, DataFrames, LinearAlgebra, TOML, Printf

const ROOT = joinpath(@__DIR__, "..", "..") |> normpath
const DESIGN = "matrix_5x3_v5"

design = TOML.parsefile(joinpath(ROOT, "experiments", "designs", DESIGN * ".toml"))
dat, prog = design["data"], design["programme"]

# -- calibration (mirrors experiments/run.jl build_reference final state) --
drops = Vector{Int}(dat["drops"])
data_full = read_data("I-O_DE2019_formatiert.csv"; datadir = ROOT)
data_v1 = retained_dataset(data_full, drops)
N = length(data_v1.factor_share)
data = recalibrate_open(data_v1; exo_scale = 1.0)
shocks0 = Shocks(ones(N), ones(N), zeros(N))
tau0 = sum(data.gov_demand)
s = data.saving_rate
cs = data.consumption_share
c0 = data.household_baseline

# -- programme incidence (mirrors run.jl programme_vectors + kept) --
n_full = length(data.factor_share) + length(drops)
kept = sort(setdiff(1:n_full, drops))
imp = CSV.read(joinpath(ROOT, prog["source"]), DataFrame)
rows = imp[imp.year .== Int(prog["year"]), :]
c1, c2 = Int(prog["column_slice"][1]), Int(prog["column_slice"][2])
raw = Matrix{Float64}(rows[1:1, c1:c2])[:]
v = raw[kept]
psi = v ./ sum(v)
g = (Float64(prog["total_eur_m"]) / data.gdp_production) .* psi
G0 = sum(g)
@printf("N=%d G0=%.8f tau0=%.8f s=%.8f\n", N, G0, tau0, s)

function stored_X(id)
    sol = CSV.read(joinpath(ROOT, "runs", id, "solution.csv"), DataFrame)
    man = TOML.parsefile(joinpath(ROOT, "runs", id, "manifest.toml"))
    me = man["metrics"]
    return [Float64.(sol.price); Float64.(sol.quantity);
            Float64(me["wage"]); Float64(get(me, "external_transfer", 0.0))], me
end

println()
println("== (i) GAMMA/DELTA vs the analytic Leontief multiplier (no solve) ==")
println("   pre-ADR-0019: 8.3e-17 / 1.7e-16, L = 1.0177957289")
for (F, mode) in (("F3", :F3), ("F2", :F2))
    ana = leontief_multiplier(data, g; mode = mode)
    for L in ("GAMMA", "DELTA")
        local lid = "matrix_5x3-v5-" * L * "-" * F
        sol = CSV.read(joinpath(ROOT, "runs", lid, "solution.csv"), DataFrame)
        q = Float64.(sol.quantity)
        @printf("  %-12s max|y - y_leontief| = %.3e   rel = %.3e   L_num = %.10f  L_ana = %.10f\n",
            L * "-" * F, maximum(abs.(q .- ana.y)),
            maximum(abs.(q .- ana.y)) / maximum(abs.(ana.y)),
            dot(data.factor_share, q), ana.L)
    end
end

println()
println("== (ii) scaled programme: ALPHA-F3 (mobile) vs GAMMA-F3 (fixed) ==")
println("   pre-ADR-0019: ALPHA stalls above k = 1, GAMMA scales cleanly")
function stored_X_nofail(id)
    try
        return stored_X(id)
    catch e
        return nothing
    end
end
for (tag, closure) in (("ALPHA", :mobile), ("GAMMA", :fixed))
    r = stored_X_nofail("matrix_5x3-v5-" * tag * "-F3")
    r === nothing && continue
    X0, _ = r
    warm = X0
    for k in (1.0, 2.0, 5.0, 10.0)
        m = mobile_labor_model(data, shocks0, 0.5, 0.5, 0.9, 1.0;
            closure = closure, financing = ExternalDebt(k .* g))
        Xinit = closure == :fixed ? warm[1:2N] : warm
        try
            sol = solve(m; init = Xinit)
            p, q, w = sol.prices_raw, sol.quantities, sol.wages_raw[1]
            Xc = closure == :fixed ? [p; q] : [p; q; w; sol.external_transfer]
            Ld = sum(sectoral_labor_demand(p, q, w, m))
            @printf("%s F3 xk%-4.0f: resid = %.2e  L = %.8f  max|dp| = %.2e  welfare = %.8f  F = %+.6e\n",
                tag, k, maximum(abs, equilibrium_residuals(m, Xc)),
                Ld, maximum(abs.(p .- 1)), real_consumption(sol), sol.external_transfer)
            warm = closure == :fixed ? [p; q; 1.0; 0.0] : [p; q; w; sol.external_transfer]
        catch e
            @printf("%s F3 xk%-4.0f: FAILED (%s)\n", tag, k, first(split(sprint(showerror, e), "\n")))
        end
    end
end

println()
println("== (iii) tilt-based Tornqvist closed form vs v5 manifests ==")
E0 = (1 - tau0)
pred_F2 = (1 - tau0 - G0) / (1 - tau0)
pred_F3 = 1.0
function tilt()
    pos = c0 .> 0
    psi1 = psi .* pos
    psi1 ./= sum(psi1)
    return 1.0 .+ G0 .* psi1 ./ max.(c0, 1e-12)
end
d = tilt()
agg = sum(cs .* d)
c1_F1 = (1 - s) .* E0 .* (cs .* d) ./ agg
pred_F1 = tornqvist_quantity_index(ones(N), c1_F1, ones(N), c0)
@printf("closed form (no solve): F1 = %.12f  F2 = %.12f  F3 = %.12f\n",
    pred_F1, pred_F2, pred_F3)
println("stored v5 manifests (consumption_rel+1 = welfare vs ref; gdp_rel+1 = income GDP vs ref):")
for L in ("BF", "ALPHA", "BETA"), F in ("F1", "F2", "F3")
    local lid = "matrix_5x3-v5-" * L * "-" * F
    man = TOML.parsefile(joinpath(ROOT, "runs", lid, "manifest.toml"))
    me = man["metrics"]
    pred = F == "F1" ? pred_F1 : F == "F2" ? pred_F2 : pred_F3
    cons = Float64(me["consumption"])
    @printf("  %-10s cons = %.12f  pred = %.12f  diff = %+.3e   gdp = %.12f   L = %.12f\n",
        L * "-" * F, cons, pred, cons - pred,
        Float64(me["gdp"]), Float64(me["employment"]))
end
a = d ./ agg
active = cs .> 0
lnQ_closed = 0.5 * (sum(cs[active] .* log.(a[active])) + sum(cs[active] .* a[active] .* log.(a[active])))
@printf("F1 second-order check: ln Q (tilt) = %+.10e vs Tornqvist %+.10e (diff %.2e)\n",
    lnQ_closed, log(pred_F1), lnQ_closed - log(pred_F1))
println()
println("VERDICT (iii): diff ~ 0 for a cell => the tilt closed form still describes")
println("the measured welfare there. A nonzero diff under ADR-0018/ADR-0019 means the")
println("closed form needs re-derivation (E now carries F; GDP and welfare are split).")
