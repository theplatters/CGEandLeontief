# probe20_pricing_door.jl -- ROADMAP Phase 7 ("the pricing door") measured
# before any ADR: Route 1 (capacity / Kaldor-Verdoorn) and Route 2 Stage A
# (the retained markup), on all three financing columns.
#
# READ-ONLY probe: no run, no kernel change, nothing committed.
#   julia --project=. experiments/probes/probe20_pricing_door.jl [--csv DIR]
# (--csv writes the ladder and the F2 sectoral cross-section as CSV.)
#
# Route 1: A_eff_i = (y_i/lambda_i)^(-delta) in the unit-cost hook, solved as
#   the damped fixed point of probe17 (A <- A .* (y/lambda)^(-delta), re-solve).
# Route 2 Stage A: p_i = mu_i c_i, mu_i = 1 + kappa (y_i/lambda_i - 1), rents
#   retained (household budget unchanged). Emulated WITHOUT a kernel change:
#   every FOC term of the kernel demand hook carries p_u^eps * y_u (labour
#   demand and the intermediate bill alike), so evaluating the hook at
#   y' = y ./ mu.^eps puts the marginal cost c = p/mu into the FOCs while the
#   purchaser prices p stay in the price indices. The same y' is passed to the
#   canary and to gdp_components, whose bill factors carry p^eps too.
#
# Parts:
#   A  fixed-wage ladders  delta in -0.75:0.125:0.75, kappa in 0:0.125:0.75
#   B  Cobb-Douglas pass-through canary (theta = eps = 1): the price block is
#      log p = (I - diag(1-fs) Omega)^(-1) r exactly, r = log mu (markup) or
#      delta log(y/lambda) (capacity); with a CONSTANT uniform markup mu_bar
#      the prices are closed-form and y-free. The ROADMAP's "p = mu" reading
#      ignores the intermediate cascade; its gap is printed.
#   C  mobile-gauge robustness row (eta = 1, ALPHA, CPI = 1): delta in
#      {-0.5, 0.5}, kappa in {0.25, 0.5}; here w/P carries the price response.
#
# Canaries: delta = 0 and kappa = 0 reproduce matrix_5x3-v10-GAMMA-* (fixed)
# and the v10 ALPHA cells (mobile); the external-account gap is ~0 on Route 1
# (rent-free) and equals -Pi_ret on Route 2 (the retained-profit leak).

using OpenCGE
using OpenCGE: _mobile_market_demand, Solution, tornqvist_quantity_index, preference_weights,
    household_expenditure, gdp_components, external_balance_canary, gdp_deflator,
    sectoral_labor_demand, equilibrium_residuals, real_consumption, labor_market_residual,
    labor_closure
using CSV, DataFrames, LinearAlgebra, Printf
import NonlinearSolve

const ROOT = joinpath(@__DIR__, "..", "..") |> normpath
const CSVDIR = let i = findfirst(==("--csv"), ARGS); i === nothing ? nothing : ARGS[i+1] end

data = recalibrate_open(retained_dataset(read_data("I-O_DE2019_formatiert.csv"; datadir = ROOT),
    Vector{Int}([])); exo_scale = 1.0)
N = length(data.factor_share)

imp = CSV.read(joinpath(ROOT, "cbase2/data_raw/impulses.csv"), DataFrame)
rows = imp[imp.year .== 2024, :]
ψ = Matrix{Float64}(rows[1:1, 3:73])[:][collect(1:N)]
ψ = ψ ./ sum(ψ)
g = (40300.0 / data.gdp_production) .* ψ
sector_names = String.(names(imp)[3:73])[1:N]

shocks0 = Shocks(ones(N), ones(N), zeros(N))
ref_sol = solve(mobile_labor_model(data, shocks0, 0.5, 0.5, 0.9, 1.0);
    init = [ones(N); data.λ; 1.0; 0.0])
ref_cons = real_consumption(ref_sol)

function f1_tilt(baseline, ψ, g)
    pos = baseline .> 0
    ψ1 = ψ .* pos; ψ1 = ψ1 ./ sum(ψ1)
    return 1.0 .+ sum(g) .* ψ1 ./ max.(baseline, 1e-12)
end
fin_of(id) = id == "F1" ? PreferenceReallocation(f1_tilt(data.household_baseline, ψ, g)) :
             id == "F2" ? TaxFinanced(g) : ExternalDebt(g)
cell_model(fin; A = ones(N), closure = :fixed, θ = 0.5, ϵ = 0.5) =
    mobile_labor_model(data, Shocks(A, ones(N), zeros(N)), θ, ϵ, 0.9, 1.0;
        closure = closure, financing = fin)
ismobile(m) = m.options.closure != :fixed

"Solution built from (p, y_foc, w, F) exactly as the kernel's solvers build it."
function make_solution(m, p, yfoc, w, F)
    σ = m.options.elasticities.σ
    L_i = sectoral_labor_demand(p, yfoc, w, m)
    fin = m.financing
    ds = preference_weights(fin, m.shocks.demand_shock)
    E = household_expenditure(fin, m, w * sum(L_i), p, sum(L_i); external_transfer = F)
    agg = sum(data.consumption_share .* ds .* p .^ (1 - σ))
    cons = (1 .- data.saving_rate) .* (data.consumption_share .* ds .* E .* p .^ (-σ)) ./ agg
    rg = tornqvist_quantity_index(p, cons, ones(N), data.household_baseline)
    num = (data.consumption_share' * p .^ (1 - σ))^(1 / (1 - σ))
    Solution(p, yfoc, fill(w, N), cons, num, rg, w * sum(L_i), m; external_transfer = F)
end

"Cell metrics; `y` is real output, `yfoc` the FOC-site quantity (= y off Route 2)."
function metrics(m, p, y, yfoc, w, F; μ = ones(N))
    sol = make_solution(m, p, yfoc, w, F)
    X = ismobile(m) ? [p; yfoc; w; F] : [p; yfoc]
    ngdp = sum(gdp_components(m, sol).V)
    Π = sum((μ .- 1) .* (p ./ μ) .* y)
    (; maxdp = maximum(abs, p .- 1), meandp = sum(data.λ .* (p .- 1)) / sum(data.λ),
       deflator = gdp_deflator(sol, ref_sol), realwage = w / sol.numeraire,
       L = sum(sectoral_labor_demand(p, yfoc, w, m)),
       cons = real_consumption(sol) / ref_cons - 1,
       canary = external_balance_canary(m, X).diff, profit = Π, profit_share = Π / ngdp)
end

# ── Route 1: capacity / Kaldor-Verdoorn ───────────────────────────────────────
function capacity(fin, δ; closure = :fixed, θ = 0.5, ϵ = 0.5, maxit = 200, tol = 1e-12, damp = 0.6)
    A = ones(N); init = nothing; lastΔ = NaN
    for it in 1:maxit
        m = cell_model(fin; A, closure, θ, ϵ)
        sol = try solve(m; init = init) catch; return nothing, it end
        y = sol.quantities
        Anew = exp.((1 - damp) .* log.(A) .+ damp .* log.((y ./ data.λ) .^ (-δ)))
        lastΔ = maximum(abs, Anew .- A)
        any(!isfinite, Anew) && return nothing, it
        init = closure == :fixed ? [sol.prices_raw; y] : [sol.prices_raw; y; sol.wages_raw[1]; sol.external_transfer]
        A = Anew
        if lastΔ < tol      # final solve at the converged A: one consistent (p, y, A)
            m = cell_model(fin; A, closure, θ, ϵ)
            sol = solve(m; init = init)
            p, y, w, F = sol.prices_raw, sol.quantities, sol.wages_raw[1], sol.external_transfer
            X = closure == :fixed ? [p; y] : [p; y; w; F]
            return (; m, p, y, yfoc = y, w, F, μ = ones(N),
                resid = maximum(abs, equilibrium_residuals(m, X))), it
        end
    end
    return nothing, maxit
end

# ── Route 2, Stage A: the retained markup ─────────────────────────────────────
# μfun(y) gives the markup vector; the FOC sites see y' = y ./ μ.^ϵ.
function markup_resid!(out, X, (m, μfun))
    ϵ = m.options.elasticities.ϵ
    p = max.(X[1:N], eps()); y = max.(X[N+1:2N], eps())
    μ = μfun(y)
    if ismobile(m)
        w = max(X[2N+1], 1e-10); F = X[2N+2]
        b = _mobile_market_demand(m, p, y ./ μ .^ ϵ, w; external_transfer = F)
        σ = m.options.elasticities.σ
        cpi = sum(data.consumption_share .* p .^ (1 - σ))^(1 / (1 - σ))
        out[2N+1] = labor_market_residual(labor_closure(m.options), m, sum(b.L_i), w, cpi)
        out[2N+2] = cpi - 1
    else
        b = _mobile_market_demand(m, p, y ./ μ .^ ϵ, 1.0)
    end
    out[1:N] .= p .- μ .* b.cost
    out[N+1:2N] .= y .- b.intermediary_demand .- b.total_final_demand
    nothing
end
function markup(fin, μfun; closure = :fixed, θ = 0.5, ϵ = 0.5, init = nothing)
    m = cell_model(fin; closure, θ, ϵ)
    if init === nothing
        s = solve(m)
        init = closure == :fixed ? [s.prices_raw; s.quantities] :
            [s.prices_raw; s.quantities; s.wages_raw[1]; s.external_transfer]
    end
    r = zeros(length(init)); x = init
    best() = (markup_resid!(r, x, (m, μfun)); maximum(abs, r))
    for alg in (NonlinearSolve.NewtonRaphson(), NonlinearSolve.LevenbergMarquardt())
        best() < 1e-12 && break
        res = NonlinearSolve.solve(NonlinearSolve.NonlinearProblem(markup_resid!, x, (m, μfun)), alg;
            reltol = 1e-13, abstol = 1e-13, maxiters = 2000)
        markup_resid!(r, res.u, (m, μfun)); rnew = maximum(abs, r)
        rnew < best() && (x = res.u)
    end
    resid = best()
    p, y = x[1:N], x[N+1:2N]
    w, F = closure == :fixed ? (1.0, 0.0) : (x[2N+1], x[2N+2])
    μ = μfun(y)
    (; m, p, y, yfoc = y ./ μ .^ ϵ, w, F, μ, resid, x)
end
μ_kappa(κ) = y -> 1 .+ κ .* (y ./ data.λ .- 1)

# ── Part A: fixed-wage ladders ────────────────────────────────────────────────
ladder = DataFrame()
sectoral = DataFrame(sector = sector_names, psi = ψ, lambda = data.λ)
function record!(route, fin_id, par, r, it)
    if r === nothing
        @printf("  %-8s %s %+6.3f  --- no convergence (it = %d)\n", route, fin_id, par, it)
        return
    end
    mt = metrics(r.m, r.p, r.y, r.yfoc, r.w, r.F; μ = r.μ)
    push!(ladder, (; route, gauge = ismobile(r.m) ? "mobile" : "fixed", fin = fin_id, param = par,
        iters = it, resid = r.resid, mt...))
    @printf("  %-8s %s %+6.3f %4d  %.1e  %.4f  %+.4f  %.4f  %.4f  %.5f  %+.5f  %+.2e  %.4f\n",
        route, fin_id, par, it, r.resid, mt.maxdp, mt.meandp, mt.deflator, mt.realwage, mt.L,
        mt.cons, mt.canary, mt.profit_share)
end
header = "  route    fin  param   it  resid    max|p-1| mean dp  deflator w/P     L        cons      canary     Pi/Y"

println("N = ", N, "  programme sectors = ", count(>(0.0), ψ), "  G0 = ", @sprintf("%.5f", sum(g)))
println("\nPart A -- fixed-wage gauge (GAMMA), delta and kappa ladders")
println(header)
for fin_id in ("F1", "F2", "F3")
    fin = fin_of(fin_id)
    for δ in -0.75:0.125:0.75
        r, it = capacity(fin, δ)
        record!("capacity", fin_id, δ, r, it)
        if fin_id == "F2" && δ in (-0.5, 0.5) && r !== nothing
            sectoral[!, "p_cap$(δ)"] = r.p; sectoral[!, "gap_cap$(δ)"] = r.y ./ data.λ .- 1
        end
    end
    init = nothing
    for κ in 0.0:0.125:0.75
        r = markup(fin, μ_kappa(κ); init); init = r.x
        record!("markup", fin_id, κ, r, 0)
        if fin_id == "F2" && κ == 0.5
            sectoral[!, "p_mk0.5"] = r.p; sectoral[!, "gap_mk0.5"] = r.y ./ data.λ .- 1
        end
    end
end

# ── Part B: Cobb-Douglas pass-through canary ──────────────────────────────────
println("\nPart B -- Cobb-Douglas pass-through canary (theta = eps = 1, F2)")
B = Diagonal(1 .- data.factor_share) * data.Ω_raw
cascade = (I - B) \ ones(N)
@printf("  intermediate cascade (I-B)^-1 1: min %.3f  lambda-mean %.3f  max %.3f\n",
    minimum(cascade), sum(data.λ .* cascade) / sum(data.λ), maximum(cascade))
println("  case                         max|log p - (I-B)^-1 r|   max|log p - log mu| (ROADMAP p = mu reading)")
for (label, solvecase) in (
        ("capacity delta=+0.5", () -> capacity(fin_of("F2"), 0.5; θ = 1.0, ϵ = 1.0)[1]),
        ("capacity delta=-0.5", () -> capacity(fin_of("F2"), -0.5; θ = 1.0, ϵ = 1.0)[1]),
        ("capacity delta=-0.5 damp 0.2", () -> capacity(fin_of("F2"), -0.5; θ = 1.0, ϵ = 1.0, damp = 0.2, maxit = 1000)[1]),
        ("markup kappa=0.5", () -> markup(fin_of("F2"), μ_kappa(0.5); θ = 1.0, ϵ = 1.0)),
        ("markup mu_bar=1.05", () -> markup(fin_of("F2"), y -> fill(1.05, N); θ = 1.0, ϵ = 1.0)))
    r = solvecase()
    if r === nothing
        @printf("  %-28s --- fixed point did not converge\n", label); continue
    end
    rvec = startswith(label, "capacity") ? -log.(r.m.shocks.supply_shock) : log.(r.μ)
    @printf("  %-28s %.2e  (resid %.1e)       %.4f\n", label,
        maximum(abs, log.(r.p) .- (I - B) \ rvec), r.resid, maximum(abs, log.(r.p) .- log.(r.μ)))
    if label == "markup mu_bar=1.05"
        @printf("  constant mu_bar = 1.05 closed form: max|p - 1.05^cascade| = %.2e; lambda-mean p - 1 = %.4f\n",
            maximum(abs, r.p .- 1.05 .^ cascade), sum(data.λ .* (r.p .- 1)) / sum(data.λ))
    end
end

# ── Part C: mobile-gauge robustness row ───────────────────────────────────────
println("\nPart C -- mobile gauge (eta = 1, ALPHA, CPI = 1): w/P carries the response")
println(header)
for fin_id in ("F1", "F2", "F3")
    fin = fin_of(fin_id)
    for δ in (0.0, -0.5, 0.5)
        r, it = capacity(fin, δ; closure = :mobile)
        record!("capacity", fin_id, δ, r, it)
    end
    for κ in (0.25, 0.5)
        r = markup(fin, μ_kappa(κ); closure = :mobile)
        record!("markup", fin_id, κ, r, 0)
    end
end

if CSVDIR !== nothing
    mkpath(CSVDIR)
    CSV.write(joinpath(CSVDIR, "probe20_ladder.csv"), ladder)
    CSV.write(joinpath(CSVDIR, "probe20_sectoral_F2.csv"), sectoral)
    println("\nwrote CSVs to ", CSVDIR)
end
println("\ndone.")
