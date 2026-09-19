# supply_identification.jl -- supply-shock BETA identification on the ADR-0019 kernel.
#
# READ-ONLY probe: no run is created, no src/ file is touched, nothing is
# committed. Run from the repository root:
#   julia --project=. experiments/probes/supply_identification.jl
#
# A +20% productivity shock in sector 1 (A[1] = 1.2) moves the real wage, so
# the BETA labour-supply elasticity eta_s becomes identified. For F in
# F1/F2/F3: solve ALPHA (eta_s = 0) warm-started from the stored v5 ALPHA-Fx
# solution, then BETA at eta_s in {0.5, 1.0, 2.0, 5.0} via solve_beta
# warm-started from that shocked ALPHA solution (steps = 8).
#
# Per cell: max|resid|, w, CPI, real wage w/CPI, L, implied elasticity
# log(L)/log(w/CPI) (plus the raw-wage variant), gdp_rel, consumption_rel,
# max|dp|, external_transfer, gate pass (mobile gate 1e-5).
#
# Pre-ADR-0019 pilot for comparison: w = 1.0049029745, L = 1.0024484897 /
# 1.0098299882 / 1.0247564458 for eta_s = 0.5 / 2 / 5, elasticities recovered
# exactly.

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
Lbar = sum(data.labor_share)

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

function f1_tilt(baseline, psi, g)
    pos = baseline .> 0
    psi1 = psi .* pos
    psi1 = psi1 ./ sum(psi1)
    return 1.0 .+ sum(g) .* psi1 ./ max.(baseline, 1e-12)
end
fin_of(id) = id == "F1" ? PreferenceReallocation(f1_tilt(data.household_baseline, psi, g)) :
             id == "F2" ? TaxFinanced(g) :
             id == "F3" ? ExternalDebt(g) : error("unknown financing $id")

function stored_X(id)
    sol = CSV.read(joinpath(ROOT, "runs", id, "solution.csv"), DataFrame)
    man = TOML.parsefile(joinpath(ROOT, "runs", id, "manifest.toml"))
    me = man["metrics"]
    return [Float64.(sol.price); Float64.(sol.quantity);
            Float64(me["wage"]); Float64(get(me, "external_transfer", 0.0))]
end

# cheap NoFinancing ALPHA reference for gdp/consumption relatives
shocks0 = Shocks(ones(N), ones(N), zeros(N))
ref_model = mobile_labor_model(data, shocks0, 0.5, 0.5, 0.9, 1.0)
ref_sol = solve(ref_model; init = [ones(N); data.λ; 1.0; 0.0])
ref_cons = real_consumption(ref_sol)

cpi_of(p, sig) = sum(data.consumption_share .* p .^ (1 - sig))^(1 / (1 - sig))

println("N = ", N, "  Lbar = ", Lbar, "  ref_cons = ", ref_cons)
println()
println("supply shock: A = ones(N), A[1] = 1.2; theta = eps = 0.5, sigma = 0.9, eta = 1")
println("pilot (pre-ADR-0019): w = 1.0049029745, L = 1.0024484897 / 1.0098299882 / 1.0247564458")
println("for eta_s = 0.5 / 2 / 5, elasticities recovered exactly")
println()

A = ones(N); A[1] = 1.2
shock_supply = Shocks(A, ones(N), zeros(N))

for F in ("F1", "F2", "F3")
    println("---- financing ", F, " ----")
    # shocked ALPHA warm-started from the stored v5 ALPHA-Fx solution
    mA = mobile_labor_model(data, shock_supply, 0.5, 0.5, 0.9, 1.0; financing = fin_of(F))
    x0 = stored_X("matrix_5x3-v5-ALPHA-" * F)
    sA = try
        solve(mA; init = x0)
    catch e
        println("ALPHA-", F, " shocked solve FAILED: ", first(split(sprint(showerror, e), "\n")))
        nothing
    end
    if sA === nothing
        continue
    end
    xA = [sA.prices_raw; sA.quantities; sA.wages_raw[1]; sA.external_transfer]
    for (tag, eta_s) in (("ALPHA", 0.0), ("BETA-0.5", 0.5), ("BETA-1.0", 1.0),
                         ("BETA-2.0", 2.0), ("BETA-5.0", 5.0))
        sol = try
            if eta_s == 0.0
                sA
            else
                solve_beta(data, shock_supply, 0.5, 0.5, 0.9, 1.0; eta_s = eta_s,
                    financing = fin_of(F), init = xA, steps = 8)
            end
        catch e
            println(tag, "-", F, " FAILED: ", first(split(sprint(showerror, e), "\n")))
            nothing
        end
        if sol === nothing
            continue
        end
        m = eta_s == 0.0 ? mA :
            mobile_labor_model(data, shock_supply, 0.5, 0.5, 0.9, 1.0;
                financing = fin_of(F), eta_s = eta_s)
        p, q, w, Fv = sol.prices_raw, sol.quantities, sol.wages_raw[1], sol.external_transfer
        resid = maximum(abs, equilibrium_residuals(m, [p; q; w; Fv]))
        cpi = cpi_of(p, 0.9)
        rw = w / cpi
        L = sum(sectoral_labor_demand(p, q, w, m))
        elas_rw = log(L / Lbar) / log(rw)
        elas_raw = log(L / Lbar) / log(w)
        gd = gdp_income(sol, ref_sol) - 1
        co = real_consumption(sol) / ref_cons - 1
        gate = resid < 1e-5 ? "pass" : "FAIL"
        @printf("%-9s resid=%.2e w=%.10f CPI=%.10f rw=%.10f L=%.10f elas_rw=%+.6f elas_raw=%+.6f gdp=%+.6e cons=%+.6e maxdp=%.2e F=%+.6e %s\n",
            tag, resid, w, cpi, rw, L, elas_rw, elas_raw, gd, co,
            maximum(abs.(p .- 1)), Fv, gate)
    end
    println()
end

println("INTERPRETATION (read from the elas_rw column above):")
println("  elas_rw == eta_s to solver precision => the supply-shock B-grades survive")
println("  on the ADR-0019 kernel. elas_raw is reported to show the CPI-deflation")
println("  wedge (ADR-0014: supply is on the real wage, not the raw wage).")
