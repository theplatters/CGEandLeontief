# v2_verify.jl — acceptance test for the open-absorption recalibration.
# Run from the PARENT repo root:  julia /workspace/agents/v2_verify.jl
using CSV, DataFrames, LinearAlgebra, Statistics, NonlinearSolve

CB = "/workspace/git/BFRep/(3)BeyondHulten/cbase2"
for f in ["interface.jl", "solution.jl", "ces.jl", "mobile_labor.jl", "leontief.jl", "util.jl"]
    include(joinpath(CB, "src", "core", f))
end
include(joinpath(CB, "src", "closures.jl"))
include(joinpath(CB, "src", "financing.jl"))
include(joinpath(CB, "src", "calibration.jl"))

# Documented decision (Notebook 03b): drop sector 71 ("Other personal service
# activities") -- residual catch-all, 1.8 percent of gross output, 37.5 percent
# self-loop whose price spirals under theta < 1 complementarity. The model is a
# 70-sector system; sector 71's VA and demand are excluded from the accounting.
N_FULL = 71
data_full = read_data("I-O_DE2019_formatiert.csv")
data_v1 = drop_sectors(data_full, [71])
N = 70
shocks = Shocks(ones(N), ones(N), zeros(N))

# ── 1. Injection continuation: exo_scale from s~0 to 1, warm-starting each step ──
# Each intermediate economy is internally exact (c0 and s recalibrated); the
# injections grow smoothly, tracking the good root and avoiding the spurious
# price-explosion root the direct solve found (max|p-1| ~ 991).
# The start scale is where the data-implied saving rate crosses zero
# (bisection: s(exo_scale) is monotone increasing); below it s < 0 and the
# round-gain matrix is non-contractive.
s_of(esc) = (d = recalibrate_open(data_v1, CB; exo_scale=esc, drops=[71]); d.saving_rate)
lo, hi = 0.0, 1.0
@assert s_of(hi) > 0
for _ in 1:40
    global lo, hi
    mid = (lo + hi) / 2
    s_of(mid) < 0 ? (lo = mid) : (hi = mid)
end
esc0 = hi   # smallest exo_scale with s >= 0 (within bisection tolerance)
println("continuation start: exo_scale* = ", round(esc0; digits=4),
        " (s there = ", round(s_of(esc0); digits=4), ")")

K = 6
THETAS = [2.0, 1.5, 1.2, 1.0, 0.8, 0.65, 0.5]   # gross substitutes -> the target 0.5
t_cont = time()
init_warm = nothing
data = nothing
ref = nothing
ref_sol = nothing
for k in 0:K
    global data, ref, ref_sol, init_warm
    exo_scale = esc0 + (1.0 - esc0) * k / K
    data = recalibrate_open(data_v1, CB; exo_scale=exo_scale, drops=[71])
    for θ in THETAS
        global ref, ref_sol, init_warm
        ref = mobile_labor_model(data, shocks, θ, 0.5, 0.9, 0.5)
        if init_warm === nothing
            # Very first solve: warm start from the linear fixed point
            (; Ω_raw, factor_share, consumption_share, import_margin,
             gov_demand, exo_demand, exports_demand) = data
            M0 = Ω_raw' * Diagonal(1.0 .- factor_share)
            Gk = M0 + Diagonal(1.0 .- import_margin) * consumption_share *
                 factor_share' * (1 - data.saving_rate) * (1 - sum(gov_demand))
            y0 = (I - Gk) \ ((1.0 .- import_margin) .* (gov_demand .+ exo_demand) .+ exports_demand)
            init_warm = [ones(N); y0; 1.0]
        end
        ref_sol = solve(ref; init=init_warm)
        mxp = maximum(abs.(ref_sol.prices_raw .- 1))
        init_warm = [ref_sol.prices_raw; ref_sol.quantities; ref_sol.wages_raw[1]]
        println("cont k=$k exo=", round(exo_scale; digits=4), " θ=", θ,
                ": resid = ", round(maximum(abs, equilibrium_residuals(ref, init_warm)); digits=10),
                ", w* = ", round(ref_sol.wages_raw[1]; digits=4),
                ", max|p-1| = ", round(mxp; digits=4))
        mxp > 10 && (println("  EXPLODED BRANCH at θ=$θ -- stopping the θ ladder for k=$k");
                     init_warm = init_warm; break)
    end
end
println("continuation done in ", round(time() - t_cont; digits=1), " s")
τ0 = sum(data.gov_demand)
println("v3 calibration OK: tau0 = ", round(τ0; digits=4),
        ", saving rate s = ", round(data.saving_rate; digits=4),
        ", export share = ", round(sum(data.exports_demand); digits=4),
        ", investment share = ", round(sum(data.exo_demand); digits=4))
init_warm = [ref_sol.prices_raw; ref_sol.quantities; ref_sol.wages_raw[1]]

# ── 1b. S = I + X - M canary at the v3 baseline (clamping -> tolerance) ──
L_ref = sum(sectoral_labor_demand(ref_sol.prices_raw, ref_sol.quantities, ref_sol.wages_raw[1], ref))
w_ref = ref_sol.wages_raw[1]
E_ref = household_expenditure(NoFinancing(), ref, w_ref * L_ref, ref_sol.prices_raw, L_ref)
S_lhs = data.saving_rate * E_ref
p_ref = ref_sol.prices_raw
cg_vec = (1 - data.saving_rate) * E_ref .* (data.consumption_share .* p_ref .^ (1-0.9)) ./ sum(data.consumption_share .* p_ref .^ (1-0.9))
M_rhs = dot(data.import_margin, cg_vec ./ max.(p_ref, 1e-12) .+ data.gov_demand .+ data.exo_demand)
I_X = dot(data.exo_demand .+ data.exports_demand, p_ref)
println("S = I + X - M canary: S = ", round(S_lhs; digits=6), " vs I+X-M = ",
        round(I_X - M_rhs; digits=6), " (diff ", round(S_lhs - (I_X - M_rhs); digits=6), ")")
println("ref diagnostics: w* = ", round(w_ref; digits=4), ", L = ", round(L_ref; digits=6),
        ", nominal GDP wL = ", round(w_ref * L_ref; digits=6),
        ", max|p-1| = ", round(maximum(abs.(p_ref .- 1)); digits=6))

# ── Programme calibration (as notebook 03) ──
imp = CSV.read(joinpath(CB, "data_raw", "impulses.csv"), DataFrame)
r24 = imp[imp.year .== 2024, :]
imp24_full = Matrix{Float64}(r24[1:1, 3:73])[:]
imp24 = imp24_full[1:N]   # sector 71 dropped with the sector; shares renormalized below
G0_MODEL = 40_300.0 / data.gdp_production
ψ = imp24 ./ sum(imp24)
g = G0_MODEL .* ψ
fin1 = PreferenceReallocation(fill(1.0, N))          # composition-neutral stand-in (m=1 test)
fin2 = TaxFinanced(g)
fin3 = ExternalDebt(g)

function headline(tag, model, sol)
    p, q, w = sol.prices_raw, sol.quantities, sol.wages_raw[1]
    fixed = labor_closure(model) isa FixedWageClosure
    X = fixed ? [p; q] : [p; q; w]
    rmax = maximum(abs, equilibrium_residuals(model, X))
    L = sum(sectoral_labor_demand(p, q, w, model))
    E = household_expenditure(model.financing, model, w * L, p, L)
    # v3 budget identity: gross household consumption = (1-s)E; saving sE leaks
    bud = dot(p, sol.consumption) - (1 - data.saving_rate) * E
    S_num = data.saving_rate * E
    F = external_balance(model.financing, model, p)
    println(tag, ": |resid| = ", round(rmax; digits=10), " | real_gdp = ",
            round(real_gdp(sol) / real_gdp(ref_sol) - 1, digits=6), " (rel v3 ref) | L = ",
            round(L; digits=6), " | E = ", round(E; digits=6), " | F = ", round(F; digits=6),
            " | S = ", round(S_num; digits=6),
            " | Σpc - (1-s)E = ", round(bud; digits=12))
    @assert rmax < 1e-6 "$tag residuals"
    @assert abs(bud) < 1e-9 "$tag budget identity (gross consumption must equal (1-s)E)"
    (tag = tag, resid = rmax, gdp_rel = real_gdp(sol) / real_gdp(ref_sol) - 1,
     employment = L, F = F)
end

res = DataFrame()
for (tag, fin) in [("F1_preference_reallocation", fin1), ("F2_tax_financed", fin2), ("F3_external_debt", fin3)]
    mdl = mobile_labor_model(data, shocks, 1.0, 0.5, 0.9, 0.5; financing=fin)
    push!(res, headline(tag * " mobile", mdl, solve(mdl; init=init_warm)))
end

# ── 3. DELTA corner: analytic equivalence restored (finiteness!) ──
for mode in (:F2, :F3)
    fin = mode === :F2 ? fin2 : fin3
    ana = leontief_multiplier(data, g; mode=mode)
    mdl = delta_model(data, shocks; ε=1e-4, financing=fin)
    sol = solve(mdl)
    p, q, w = sol.prices_raw, sol.quantities, sol.wages_raw[1]
    rmax = maximum(abs, equilibrium_residuals(mdl, [p; q]))
    L_num = sum(sectoral_labor_demand(p, q, w, mdl))
    rel = maximum(abs.(q .- ana.y)) / maximum(abs.(ana.y))
    println("DELTA ", mode, " ε=1e-4: |resid| = ", round(rmax; digits=10),
            " | L num/ana = ", round(L_num; digits=6), "/", round(ana.L; digits=6),
            " | rel y error = ", round(rel; digits=6), " | F ana = ", round(ana.F; digits=6))
    @assert rmax < 1e-6
    @assert rel < 5e-3 "DELTA solve must match the analytic Leontief system (rel error $rel)"
    # ε-convergence
    mdl3 = delta_model(data, shocks; ε=1e-3, financing=fin)
    sol3 = solve(mdl3)
    e3 = maximum(abs.(sol3.quantities .- ana.y)); e4 = maximum(abs.(q .- ana.y))
    println("   convergence err(1e-3)/err(1e-4) = ", round(e3 / e4; digits=2), " (expect ≈ 10)")
    @assert e4 <= e3 "DELTA must converge toward the analytic system as ε → 0"
end

# ── 4. BETA elasticity (single-point identification) ──
for η_s in (0.5, 1.0)
    mdl = mobile_labor_model(data, shocks, 1.0, 0.5, 0.9, 0.5; financing=fin3, eta_s=η_s)
    sol = solve_beta(data, shocks, 1.0, 0.5, 0.9, 0.5; financing=fin3, eta_s=η_s, init=init_warm)
    w = sol.wages_raw[1]
    L = sum(sectoral_labor_demand(sol.prices_raw, sol.quantities, w, mdl))
    est = log(L / 1.0) / log(w / 1.0)
    println("BETA η_s = $η_s: implied elasticity = ", round(est; digits=5))
    @assert abs(est - η_s) < 5e-2
end

# ── 5. Cobb-Douglas guard: ϵ = 1 exact, continuous ──
cd = mobile_labor_model(data, shocks, 1.0, 1.0, 0.9, 0.5; financing=fin3)
sol_cd = solve(cd; init=init_warm)
cd_near = mobile_labor_model(data, shocks, 1.0, 1 - 1e-7, 0.9, 0.5; financing=fin3)
sol_near = solve(cd_near; init=init_warm)
println("CD guard: real_gdp(ϵ=1) = ", round(real_gdp(sol_cd); digits=8),
        " vs ϵ=1-1e-7: ", round(real_gdp(sol_near); digits=8))
@assert all(isfinite, sol_cd.prices_raw) && abs(real_gdp(sol_cd) - real_gdp(sol_near)) < 1e-5

println("\nALL V2 ACCEPTANCE TESTS PASSED")
