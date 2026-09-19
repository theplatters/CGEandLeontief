# numeraire_invariance.jl -- demand-only gauge on the ADR-0019 kernel (DE-0011 evidence).
#
# READ-ONLY probe: no run is created, no src/ file is touched, nothing is
# committed. Run from the repository root:
#   julia --project=. experiments/probes/numeraire_invariance.jl
#
# Demand-only programme (unit shocks), ALPHA (eta = 1, eta_s = 0), F2 and F3,
# scales k in {1, 2, 5, 10}, solved warm-started (chained). Per cell: w,
# max|dp|, L, gdp_rel, consumption_rel, external_transfer.
#
# Questions: (i) is w = 1 with no price response (numeraire invariance)?
# (ii) are the F2 and F3 real allocations identical (ADR-0019 financing
# neutrality)? (iii) does the old closed form consumption = 1 - k*G0/(1-tau0)
# still describe the measured consumption_rel? Measure, do not assume: under
# ADR-0019 E = (1-tau)*w*sum(L) + F changed the system.

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
@printf("N=%d tau0=%.8f G0=%.8f 1-tau0=%.8f\n", N, tau0, G0, 1 - tau0)

# cheap NoFinancing ALPHA reference for gdp/consumption relatives
ref_model = mobile_labor_model(data, shocks0, 0.5, 0.5, 0.9, 1.0)
ref_sol = solve(ref_model; init = [ones(N); data.λ; 1.0; 0.0])
ref_cons = real_consumption(ref_sol)

function stored_X(id)
    sol = CSV.read(joinpath(ROOT, "runs", id, "solution.csv"), DataFrame)
    man = TOML.parsefile(joinpath(ROOT, "runs", id, "manifest.toml"))
    me = man["metrics"]
    return [Float64.(sol.price); Float64.(sol.quantity);
            Float64(me["wage"]); Float64(get(me, "external_transfer", 0.0))]
end

println()
println("== demand-only gauge: ALPHA eta_s = 0, F2/F3, k = 1/2/5/10 ==")
println("old closed form: cons(F2) = 1 - k*G0/(1-tau0); cons(F3) = 1 (households untouched)")
println()
warm = stored_X("matrix_5x3-v5-ALPHA-F3")
neutral = Dict{Float64,Any}()
for k in (1.0, 2.0, 5.0, 10.0)
    gk = k .* g
    sols = Dict{String,Any}()
    for F in ("F2", "F3")
        fin = F == "F2" ? TaxFinanced(gk) : ExternalDebt(gk)
        m = mobile_labor_model(data, shocks0, 0.5, 0.5, 0.9, 1.0; financing = fin)
        sol = try
            solve(m; init = warm)
        catch e
            println("ALPHA-", F, " xk", k, " FAILED: ", first(split(sprint(showerror, e), "\n")))
            nothing
        end
        if sol === nothing
            continue
        end
        p, q, w, Fv = sol.prices_raw, sol.quantities, sol.wages_raw[1], sol.external_transfer
        resid = maximum(abs, equilibrium_residuals(m, [p; q; w; Fv]))
        L = sum(sectoral_labor_demand(p, q, w, m))
        gd = gdp_income(sol, ref_sol) - 1
        co = real_consumption(sol) / ref_cons - 1
        pred = F == "F2" ? 1 - k * G0 / (1 - tau0) : 1.0
        gate = resid < 1e-5 ? "pass" : "FAIL"
        @printf("ALPHA-%s xk%-4.0f resid=%.2e w=%.10f maxdp=%.2e L=%.10f gdp=%+.6e cons=%+.6e (pred %+.6e, diff %+.2e) F=%+.6e Bgov=%.6f %s\n",
            F, k, resid, w, maximum(abs.(p .- 1)), L, gd, co, pred - 1, co - (pred - 1),
            Fv, (F == "F3" ? dot(p, gk) : 0.0), gate)
        sols[F] = (p = p, q = q, w = w, F = Fv, sol = sol, model = m)
        global warm = [p; q; w; Fv]
    end
    if haskey(sols, "F2") && haskey(sols, "F3")
        a, b = sols["F2"], sols["F3"]
        @printf("   F2-vs-F3 real allocation: max|dp2-dp3|=%.2e max|dq|=%.2e |dw|=%.2e |dL|=%.2e |dcons|=%.2e |dnet|=%.2e\n",
            maximum(abs.(a.p .- b.p)), maximum(abs.(a.q .- b.q)), abs(a.w - b.w),
            abs(sum(sectoral_labor_demand(a.p, a.q, a.w, a.model)) - sum(sectoral_labor_demand(b.p, b.q, b.w, b.model))),
            abs(real_consumption(a.sol) - real_consumption(b.sol)),
            abs(a.F - (b.F + dot(b.p, gk))))
        neutral[k] = (F2 = a.F, F3 = b.F, Bgov = dot(b.p, gk))
    end
end

println()
println("== financing-neutrality check: F_F3 =?= F_F2 - B_gov, net position F + B_gov ==")
for k in (1.0, 2.0, 5.0, 10.0)
    if haskey(neutral, k)
        t = neutral[k]
        @printf("  k=%-4.0f F2=%+.6e F3=%+.6e Bgov=%.6f F3-(F2-Bgov)=%+.2e net2=%+.6e net3=%+.6e\n",
            k, t.F2, t.F3, t.Bgov, t.F3 - (t.F2 - t.Bgov), t.F2, t.F3 + t.Bgov)
    end
end
println()
println("VERDICTS (read from the columns above):")
println("  w = 1 and maxdp ~ 0 at every k => numeraire invariance holds.")
println("  max|dq|(F2,F3) ~ 0 => ADR-0019 financing neutrality holds.")
println("  cons - pred ~ 0 for F2 (and cons ~ 0 for F3) => the old closed form survives;")
println("  otherwise report the measured drift (E now carries F).")
