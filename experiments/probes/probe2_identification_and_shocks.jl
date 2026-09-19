# probe2_identification_and_shocks.jl — is the BF row identified, and what
# would move the real wage?
#
# Run from the repository root:
#   julia --project=. experiments/probes/probe2_identification_and_shocks.jl
#
# Four questions, all answered by re-solving (never by editing src/):
#
#   A. BF (eta = 0) identification. Under eta = 0 the labour-market equation
#      is `sum(labor_share) - labor_bar`, which is identically zero, so the
#      mobile system has 2N effective equations for 2N+1 unknowns. Solve the
#      same cell from two different warm starts and see whether the sectoral
#      quantities move while every reported aggregate stays put.
#   B. The BF factor-market gap: the wage bill implied by production
#      (cost-minimising labour at the solved y) minus the wage income the
#      closure pays out (sum(labor_share)). This gap is why the canary
#      identity fails in the BF rows (probe1: gap up to 7.8e-3).
#   C. Demand-independence of the price block: re-solve ALPHA with a 10x
#      programme. If p = 1 and w = 1 survive, the price block is atechnology
#      block and no demand-side modification of the matrix can move it.
#   D. The channel that DOES move the real wage: a +20 % productivity shock in
#      sector 1 under BETA with eta_s in {0.5, 2, 5} against ALPHA.

using BeyondHulten
using CSV, DataFrames, LinearAlgebra, TOML, Printf

const ROOT = joinpath(@__DIR__, "..", "..") |> normpath
const DESIGN = "matrix_5x3_v3"

design = TOML.parsefile(joinpath(ROOT, "experiments", "designs", DESIGN * ".toml"))
dat, prog = design["data"], design["programme"]

drops = Vector{Int}(dat["drops"])
data_full = read_data(String(dat["source_table"]) |> f -> last(splitpath(f)); datadir = ROOT)
data_v1 = retained_dataset(data_full, drops)
N = length(data_v1.factor_share)
data = recalibrate_open(data_v1; exo_scale = 1.0)

imp = CSV.read(joinpath(ROOT, prog["source"]), DataFrame)
rows = imp[imp.year .== prog["year"], :]
c1, c2 = Int.(prog["column_slice"])
v = Float64.(Matrix(rows[1:1, c1:c2])[:])[1:N]
ψ = v ./ sum(v)
g = (Float64(prog["total_eur_m"]) / data.gdp_production) .* ψ
G0 = sum(g)

function f1_tilt(baseline, ψ, g)
    pos = baseline .> 0
    ψ1 = ψ .* pos
    ψ1 = ψ1 ./ sum(ψ1)
    return 1.0 .+ sum(g) .* ψ1 ./ max.(baseline, 1e-12)
end
financing(id) = id == "F1" ? PreferenceReallocation(f1_tilt(data.household_baseline, ψ, g)) :
                id == "F2" ? TaxFinanced(g) :
                id == "F3" ? ExternalDebt(g) : error("unknown financing $id")

# warm start: the stored ALPHA baseline-near solution of the same cell
function stored_X(id)
    sol = CSV.read(joinpath(ROOT, "runs", id, "solution.csv"), DataFrame)
    man = TOML.parsefile(joinpath(ROOT, "runs", id, "manifest.toml"))
    return [Float64.(sol.price); Float64.(sol.quantity); Float64(man["metrics"]["wage"])]
end

function report(tag, model, sol)
    fixed = labor_closure(model.options) isa FixedWageClosure
    p, q, w = sol.prices_raw, sol.quantities, sol.wages_raw[1]
    X = fixed ? [p; q] : [p; q; w]
    resid = maximum(abs, equilibrium_residuals(model, X))
    Ld = sum(sectoral_labor_demand(p, q, w, model))
    Lc = sum(BeyondHulten._positive_floor(BeyondHulten._cost_minimizing_labor(p, q, w, model)))
    cpi_v = sum(data.consumption_share .* p .^ (1 - 0.9))^(1 / (1 - 0.9))
    @printf("%-26s |resid|=%.2e w=%.10f realw=%.10f L_d=%.10f L_cost=%.10f maxdp=%.2e realGDP=%.10f\n",
        tag, resid, w, w / cpi_v, Ld, Lc, maximum(abs.(p .- 1)), real_gdp(sol))
    return (p = p, q = q, w = w, Ld = Ld, Lc = Lc, resid = resid, rgdp = real_gdp(sol))
end

println("N = ", N, "  G0 = ", round(G0; digits = 6), "  s = ", round(data.saving_rate; digits = 6))

# ── A. BF identification: two warm starts, same cell ──────────────────────
println("\n── A. BF (eta = 0) solved from two different warm starts ──")
cell = design["cells"]["matrix_5x3-v3-BF-F3"]
θ, ϵ, σ, η = Float64(cell["theta"]), Float64(cell["epsilon"]), Float64(cell["sigma"]), Float64(cell["eta"])
mdl = mobile_labor_model(data, Shocks(ones(N), ones(N), zeros(N)), θ, ϵ, σ, η; financing = financing("F3"))

x0 = stored_X("matrix_5x3-v3-ALPHA-F3")
x1 = copy(x0)
x1[N+1:2N] .*= (1 .+ 0.01 .* cos.((1:N) .* 0.7))   # 1 % jitter on quantities
s0 = solve(mdl; init = x0)
s1 = solve(mdl; init = x1)
r0, r1 = report("BF-F3 warm start A", mdl, s0), report("BF-F3 warm start B", mdl, s1)
@printf("  quantity difference: max|dq| = %.3e   sum|dq| = %.3e\n",
    maximum(abs.(r0.q .- r1.q)), sum(abs.(r0.q .- r1.q)))
@printf("  aggregate difference: |dw| = %.3e   |dL| = %.3e   |d realGDP| = %.3e\n",
    abs(r0.w - r1.w), abs(r0.Ld - r1.Ld), abs(r0.rgdp - r1.rgdp))
@printf("  labour-market equation at eta = 0: sum(labor_share) - labor_bar = %.3e (identically zero)\n",
    sum(data.labor_share) - mdl.options.labor_bar)

# ── B. the BF factor-market gap ───────────────────────────────────────────
println("\n── B. factor-market gap: wage bill of production vs wage income paid out ──")
for L in ("BF", "ALPHA", "BETA"), F in ("F1", "F2", "F3")
    id = "matrix_5x3-v3-" * L * "-" * F
    cl = design["cells"][id]
    m = if L == "BETA"
        mobile_labor_model(data, Shocks(ones(N), ones(N), zeros(N)), Float64(cl["theta"]),
            Float64(cl["epsilon"]), Float64(cl["sigma"]), Float64(cl["eta"]);
            financing = financing(F), eta_s = Float64(cl["eta_s"]))
    else
        mobile_labor_model(data, Shocks(ones(N), ones(N), zeros(N)), Float64(cl["theta"]),
            Float64(cl["epsilon"]), Float64(cl["sigma"]), Float64(cl["eta"]); financing = financing(F))
    end
    X = stored_X(id)
    p, q, w = X[1:N], X[N+1:2N], X[2N+1]
    Lc = sum(BeyondHulten._positive_floor(BeyondHulten._cost_minimizing_labor(p, q, w, m)))
    Ld = sum(sectoral_labor_demand(p, q, w, m))
    canary = external_balance_canary(m, X)
    pm = dot(p, market_clearing_residuals(m, X))
    @printf("%-20s L_demand = %.10f  L_costmin(y) = %.10f  gap = %+.3e  (p.mktr - canary = %+.3e)\n",
        L * "-" * F, Ld, Lc, Lc - Ld, pm - canary.diff)
end

# ── C. demand-independence of the price block (10x programme) ─────────────
println("\n── C. ALPHA with a 10x programme: does the price block move? ──")
# Analytic prediction (see the session log): p = 1 and w = 1 are pinned by the
# zero-profit block and the CPI numeraire alone, so under a demand-only shock
# E = (1 - tau0 - k*G0) under F2 and the consumption Tornqvist index is
# E/E0 = 1 - k*G0/(1 - tau0); under F3 households are untouched and it is 1.
function run_C(x_start)
    τ0 = sum(data.gov_demand)
    warm = x_start
    for k in (1.0, 2.0, 5.0, 10.0)
        gk = k .* g
        for F in ("F2", "F3")
            fin = F == "F2" ? TaxFinanced(gk) : ExternalDebt(gk)
            m = mobile_labor_model(data, Shocks(ones(N), ones(N), zeros(N)), 0.5, 0.5, 0.9, 1.0; financing = fin)
            r = try
                report("ALPHA-$(F) x$(k)", m, solve(m; init = warm))
            catch e
                println("ALPHA-$(F) x$(k)  solve failed: ", sprint(showerror, e))
                nothing
            end
            if r !== nothing
                warm = [r.p; r.q; r.w]
                pred = F == "F2" ? 1 - k * G0 / (1 - τ0) : 1.0
                @printf("      analytic prediction (p = w = 1): realGDP = %.10f   difference = %+.2e\n",
                    pred, r.rgdp - pred)
            end
        end
    end
end
run_C(x0)

# ── D. the channel that moves the real wage: a supply shock ───────────────
println("\n── D. +20 % productivity in sector 1: BETA (eta_s) vs ALPHA ──")
A = ones(N); A[1] = 1.2
shock_supply = Shocks(A, ones(N), zeros(N))
for F in ("F1", "F2", "F3")
    for (tag, eta_s) in (("ALPHA", 0.0), ("BETA eta_s=0.5", 0.5), ("BETA eta_s=2.0", 2.0), ("BETA eta_s=5.0", 5.0))
        m = mobile_labor_model(data, shock_supply, 0.5, 0.5, 0.9, 1.0;
            financing = financing(F), eta_s = eta_s)
        s = begin
            if eta_s == 0.0
                solve(m; init = x0)
            else
                solve_beta(data, shock_supply, 0.5, 0.5, 0.9, 1.0; eta_s = eta_s,
                    financing = financing(F), init = x0, steps = 8)
            end
        end
        r = report("$(tag) $(F) A1=1.2", m, s)
        @printf("      implied elasticity log(L)/log(w) = %+.6f\n", log(r.Ld) / log(r.w))
    end
end
