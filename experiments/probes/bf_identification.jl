# bf_identification.jl -- eta = 0 (BF) identification on the ADR-0019 kernel.
#
# READ-ONLY probe: no run is created, no src/ file is touched, nothing is
# committed. Run from the repository root:
#   julia --project=. experiments/probes/bf_identification.jl
#
# Under ADR-0019 the mobile system is [p; y; w; F] = 2N+2 unknowns and the
# degenerate eta = 0 labour row is replaced by the explicit pin F = 0, so the
# old one-dimensional indeterminacy (probe2 section A on the pre-ADR-0019
# kernel) may now be closed. This script re-measures it:
#
#   A1. Solve each BF-F1/F2/F3 cell from (i) the stored v5 BF solution and
#       (ii) that same vector with 1% cosine jitter on quantities, and report
#       max|dq|, sum|dq|, |dw|, |dL|, |d gdp_rel|, |d consumption_rel|, |dF|.
#       Control: ALPHA-F3 from the same two starts.
#   A2. The factor-market gap L_costmin(y) - sum(labor_share) at the stored v5
#       solutions (BF/ALPHA/BETA x F1/F2/F3), which the canary carries at eta=0.

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

function f1_tilt(baseline, psi, g)
    pos = baseline .> 0
    psi1 = psi .* pos
    psi1 = psi1 ./ sum(psi1)
    return 1.0 .+ sum(g) .* psi1 ./ max.(baseline, 1e-12)
end
fin_of(id) = id == "F1" ? PreferenceReallocation(f1_tilt(data.household_baseline, psi, g)) :
             id == "F2" ? TaxFinanced(g) :
             id == "F3" ? ExternalDebt(g) : error("unknown financing $id")

function cell_model(lab, fin_id; eta_s = 0.0)
    fin = fin_of(fin_id)
    if lab == "BF"
        return mobile_labor_model(data, shocks0, 0.5, 0.5, 0.9, 0.0; financing = fin)
    elseif lab == "ALPHA"
        return mobile_labor_model(data, shocks0, 0.5, 0.5, 0.9, 1.0; financing = fin)
    elseif lab == "BETA"
        return mobile_labor_model(data, shocks0, 0.5, 0.5, 0.9, 1.0; financing = fin, eta_s = eta_s)
    end
    error("unknown labor $lab")
end

# stored v5 solution as a canonical 2N+2 vector [p; q; w; F]
function stored_X(id)
    sol = CSV.read(joinpath(ROOT, "runs", id, "solution.csv"), DataFrame)
    man = TOML.parsefile(joinpath(ROOT, "runs", id, "manifest.toml"))
    me = man["metrics"]
    return [Float64.(sol.price); Float64.(sol.quantity);
            Float64(me["wage"]); Float64(get(me, "external_transfer", 0.0))], me
end

# cheap NoFinancing ALPHA reference (single solve from the baseline init;
# verified: resid 4.4e-16, p = w = 1, F = 0) for gdp/consumption relatives
ref_model = mobile_labor_model(data, shocks0, 0.5, 0.5, 0.9, 1.0)
ref_sol = solve(ref_model; init = [ones(N); data.λ; 1.0; 0.0])
ref_cons = real_consumption(ref_sol)
@printf("reference: resid=%.2e w=%.10f F=%.2e cons0=%.12f Lbar=%.12f\n",
    maximum(abs, equilibrium_residuals(ref_model,
        [ref_sol.prices_raw; ref_sol.quantities; ref_sol.wages_raw[1]; ref_sol.external_transfer])),
    ref_sol.wages_raw[1], ref_sol.external_transfer, ref_cons, sum(data.labor_share))

function unpack(sol)
    return sol.prices_raw, sol.quantities, sol.wages_raw[1], sol.external_transfer
end
function L_of(model, p, q, w)
    return sum(sectoral_labor_demand(p, q, w, model))
end

println()
println("== A1. BF two-start identification (stored v5 BF solution vs 1% jittered) ==")
for F in ("F1", "F2", "F3")
    id = "matrix_5x3-v5-BF-" * F
    mdl = cell_model("BF", F)
    x0, me0 = stored_X(id)
    x1 = copy(x0)
    x1[N+1:2N] .*= (1 .+ 0.01 .* cos.((1:N) .* 0.7))
    s0 = try solve(mdl; init = x0) catch e; println(id, " start A FAILED: ", sprint(showerror, e)); nothing end
    s1 = try solve(mdl; init = x1) catch e; println(id, " start B FAILED: ", sprint(showerror, e)); nothing end
    if s0 === nothing || s1 === nothing
        continue
    end
    p0, q0, w0, F0 = unpack(s0)
    p1, q1, w1, F1v = unpack(s1)
    r0 = maximum(abs, equilibrium_residuals(mdl, [p0; q0; w0; F0]))
    r1 = maximum(abs, equilibrium_residuals(mdl, [p1; q1; w1; F1v]))
    L0, L1 = L_of(mdl, p0, q0, w0), L_of(mdl, p1, q1, w1)
    g0, g1 = gdp_income(s0, ref_sol) - 1, gdp_income(s1, ref_sol) - 1
    c0, c1v = real_consumption(s0) / ref_cons - 1, real_consumption(s1) / ref_cons - 1
    @printf("%-18s residA=%.2e residB=%.2e max|dq|=%.3e sum|dq|=%.3e |dw|=%.3e |dL|=%.3e |dgdp|=%.3e |dcons|=%.3e |dF|=%.3e\n",
        "BF-" * F, r0, r1, maximum(abs.(q0 .- q1)), sum(abs.(q0 .- q1)),
        abs(w0 - w1), abs(L0 - L1), abs(g0 - g1), abs(c0 - c1v), abs(F0 - F1v))
    @printf("%-18s   startA: w=%.10f L=%.10f F=%+.6e gdp_rel=%+.6e cons_rel=%+.6e maxdp=%.2e\n",
        "", w0, L0, F0, g0, c0, maximum(abs.(p0 .- 1)))
    @printf("%-18s   startB: w=%.10f L=%.10f F=%+.6e gdp_rel=%+.6e cons_rel=%+.6e maxdp=%.2e\n",
        "", w1, L1, F1v, g1, c1v, maximum(abs.(p1 .- 1)))
end

println()
println("== A1 control. ALPHA-F3 from stored solution vs 1% jittered quantities ==")
begin
    id = "matrix_5x3-v5-ALPHA-F3"
    mdl = cell_model("ALPHA", "F3")
    x0, _ = stored_X(id)
    x1 = copy(x0)
    x1[N+1:2N] .*= (1 .+ 0.01 .* cos.((1:N) .* 0.7))
    s0 = solve(mdl; init = x0)
    s1 = solve(mdl; init = x1)
    p0, q0, w0, F0 = unpack(s0)
    p1, q1, w1, F1v = unpack(s1)
    @printf("%-18s max|dq|=%.3e sum|dq|=%.3e |dw|=%.3e |dL|=%.3e |dF|=%.3e residA=%.2e residB=%.2e\n",
        "ALPHA-F3", maximum(abs.(q0 .- q1)), sum(abs.(q0 .- q1)),
        abs(w0 - w1), abs(L_of(mdl, p0, q0, w0) - L_of(mdl, p1, q1, w1)),
        abs(F0 - F1v),
        maximum(abs, equilibrium_residuals(mdl, [p0; q0; w0; F0])),
        maximum(abs, equilibrium_residuals(mdl, [p1; q1; w1; F1v])))
end

println()
println("== A2. factor-market gap at the stored v5 solutions ==")
println("   L_costmin(y) = cost-minimizing labour at solved (p,q,w); L_paid = sum(labor_share)")
for L in ("BF", "ALPHA", "BETA"), F in ("F1", "F2", "F3")
    local lid = "matrix_5x3-v5-" * L * "-" * F
    cl = design["cells"][lid]
    m = L == "BETA" ? cell_model("BETA", F; eta_s = Float64(cl["eta_s"])) : cell_model(L == "BF" ? "BF" : "ALPHA", F)
    X, _ = stored_X(lid)
    p, q, w = X[1:N], X[N+1:2N], X[2N+1]
    Lc = sum(BeyondHulten._positive_floor(BeyondHulten._cost_minimizing_labor(p, q, w, m)))
    Ld = sum(sectoral_labor_demand(p, q, w, m))
    canary = external_balance_canary(m, X)
    pm = dot(p, market_clearing_residuals(m, X))
    @printf("%-12s L_paid = %.10f  L_costmin(y) = %.10f  gap = %+.3e  (p.mktr - canary.diff = %+.3e, canary.diff = %+.3e)\n",
        L * "-" * F, Ld, Lc, Lc - Ld, pm - canary.diff, canary.diff)
end

println()
println("IDENTIFICATION VERDICT (read from the A1 max|dq| column above):")
println("  max|dq| ~ 1e-12 or below => the eta = 0 quantity vector is identified to")
println("  ~machine precision (the ADR-0019 F = 0 pin closed the old null direction).")
println("  max|dq| >> solver tolerance => the quantity vector still drifts between starts.")
