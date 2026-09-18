# probe3_sector_specific_wage_foresight.jl — what would sector-specific wages
# at eta = 0 buy?  (PROPOSAL SUPPORT — no src/ change, no new closure.)
#
# Run from the repository root:
#   julia --project=. experiments/probes/probe3_sector_specific_wage_foresight.jl
#
# E. Feasibility of a scaled-up F3 programme: ALPHA pins L = Lbar, so the
#    external-financing column can only be solved while the programme fits
#    inside full-employment capacity. GAMMA (employment endogenous) is the
#    control.
# F. First-order estimate of the true BF benchmark (sector-specific wages at
#    eta = 0, the recorded open gate):
#      1. employment the programme asks for, by sector (Leontief round, as in
#         `leontief_multiplier`, mode :F3);
#      2. the sectoral wage increase that clears each sectoral market when
#         labour is immobile, for a sectoral supply elasticity eps_L;
#      3. the cost push through the domestic intermediate bill,
#         dln p = (I - B')^-1 (fs .* dln w);
#      4. the CPI effect and the implied real-consumption (real GDP) effect.
#    This is a FIRST-ORDER estimate -- it is not a solve and it is not a run.

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
c1, c2 = Int.(prog["column_slice"])
v = Float64.(Matrix(rows[1:1, c1:c2])[:])[1:N]
ψ = v ./ sum(v)
g = (Float64(prog["total_eur_m"]) / data.gdp_production) .* ψ
G0 = sum(g)
τ0 = sum(data.gov_demand)

sectors = String.(data_full.io.Sektoren)[1:N]

# ── E. scaled F3: capacity limit under ALPHA vs GAMMA ─────────────────────
println("── E. scaled-up F3 programme: capacity limit ──")
function stored_X(id)
    sol = CSV.read(joinpath(ROOT, "runs", id, "solution.csv"), DataFrame)
    man = TOML.parsefile(joinpath(ROOT, "runs", id, "manifest.toml"))
    return [Float64.(sol.price); Float64.(sol.quantity); Float64(man["metrics"]["wage"])]
end
for (tag, closure, start) in (("ALPHA", :mobile, "matrix_5x3-v3-ALPHA-F3"),
                              ("GAMMA", :fixed, "matrix_5x3-v3-GAMMA-F3"))
    warm = stored_X(start)
    for k in (1.0, 2.0, 5.0, 10.0)
        m = mobile_labor_model(data, Shocks(ones(N), ones(N), zeros(N)), 0.5, 0.5, 0.9, 1.0;
            closure = closure, financing = ExternalDebt(k .* g))
        X0 = closure == :fixed ? warm[1:2N] : warm
        try
            s = solve(m; init = X0)
            p, q, w = s.prices_raw, s.quantities, s.wages_raw[1]
            Ld = sum(sectoral_labor_demand(p, q, w, m))
            @printf("%s F3 x%.0f: |resid| = %.2e  L = %.8f  max|dp| = %.2e  realGDP = %.8f\n",
                tag, k, maximum(abs, equilibrium_residuals(m, closure == :fixed ? [p; q] : [p; q; w])),
                Ld, maximum(abs.(p .- 1)), real_gdp(s))
            warm = closure == :fixed ? [p; q; 1.0] : [p; q; w]
        catch e
            @printf("%s F3 x%.0f: FAILED (%s)\n", tag, k, first(split(sprint(showerror, e), "\n")))
        end
    end
end

# ── F. what sector-specific wages at eta = 0 would do (first order) ───────
println("\n── F. first-order estimate of the true BF benchmark (eta = 0, sector wages) ──")
(; Ω_raw, factor_share, λ, A_bill, consumption_share, import_margin, labor_share) = data
fs = factor_share
B = Ω_raw' .* (A_bill ./ λ)'          # B[i,u] = domestic input of i per unit of u
B = Matrix(B)
# check: column sums are the domestic intermediate share (round-gain criterion)
@printf("domestic intermediate share: max column sum of B = %.6f (round-gain input)\n",
    maximum(vec(sum(B; dims = 1))))

# 1. employment the programme asks for (Leontief round with consumption feedback)
M = B                                            # intermediate round
gain = M + Diagonal(1.0 .- import_margin) * (consumption_share * fs') * (1 - data.saving_rate)
injg = (1.0 .- import_margin) .* g
b = injg .+ (1.0 .- import_margin) .* (data.gov_demand .+ data.exo_demand) .+ data.exports_demand .-
    (1.0 .- import_margin) .* consumption_share .* (1 - data.saving_rate) .* τ0
y0 = (I - gain) \ b                              # counterfactual gross output (F3)
y0_base = (I - gain) \ (b - injg)                # without the programme
Δy = y0 - y0_base
Ld0 = fs .* y0_base
ΔL = fs .* Δy
@printf("total employment the programme asks for: dL = %.6f  (multiplier dL/G0 = %.3f)\n",
    sum(ΔL), sum(ΔL) / G0)

order = sortperm(ΔL ./ Ld0; rev = true)
println("\n  sector                                        dL/L      L share   dln w (eps=1)   cons share")
for i in order[1:8]
    @printf("  %-42s %+8.4f%%  %7.4f%%   %+8.4f%%      %7.4f%%\n",
        first(sectors[i], 42), 100 * ΔL[i] / Ld0[i], 100 * Ld0[i] / sum(Ld0),
        100 * ΔL[i] / Ld0[i], 100 * consumption_share[i])
end

# 2./3./4. wage response, cost push, CPI and real consumption
for eps_L in (0.5, 1.0, 2.0)
    dlnw = (ΔL ./ Ld0) ./ eps_L
    dlnp = (I - B') \ (fs .* dlnw)
    dlnCPI = dot(consumption_share, dlnp)
    dlnIncome = dot(Ld0 ./ sum(Ld0), dlnw)
    @printf("\n  sectoral supply elasticity eps_L = %.1f:\n", eps_L)
    @printf("    max sectoral wage increase        %+8.4f%%  (sector %s)\n",
        100 * maximum(dlnw), first(sectors[argmax(dlnw)], 30))
    @printf("    employment-weighted wage bill     %+8.4f%%\n", 100 * dlnIncome)
    @printf("    max price increase                %+8.4f%%\n", 100 * maximum(dlnp))
    @printf("    CPI (consumption-share weights)   %+8.4f%%\n", 100 * dlnCPI)
    @printf("    first-order real consumption      %+8.4f%%\n", 100 * (dlnIncome - dlnCPI))
end
