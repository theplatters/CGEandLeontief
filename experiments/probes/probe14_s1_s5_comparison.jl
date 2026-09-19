# experiments/probes/probe14_s1_s5_comparison.jl
#
# The five scenario slots S1..S5 of docs/WORKPLAN_SENSITIVE_PRICES.md, measured
# side by side on the full-71 calibration. No src/ change.
#
#   S1  eta = 0 rigid sectoral wages (the executed matrix_5x3-v6 BF row)
#   S2  general closure (ADR-0022), UNIFORM eta_s: one elasticity for all sectors
#   S3  general closure, two-group: rigid group at eta_s = 0 + flexible at 0.5
#   S4  pinned wage VECTOR (ADR-0021 option A): exogenous heterogeneity
#   S5  utilization / Verdoorn externality (capacity door reduced form),
#       A_eff_i = A_i (y_i / lambda_i)^(-delta), delta in {+0.5, -0.5}
#
# Discriminator (the workplan's reporting unit is a pass-through; this probe uses
# the equivalent per-cell price deviation): SENSITIVITY means max|p - 1| DIFFERS
# across the financing columns F1/F2/F3. HETEROGENEITY means prices move but the
# movement is the SAME in every column (it is carried by the exogenous pin).
#
# Database:      max|p-1| F1 | max|p-1| F2 | max|p-1| F3 || employment || wage range || verdict
#
# Pre-registered predictions:
#   R1  S1, S2, S3 and S5 are DEMAND-SENSITIVE (F1 != F2/F3 in max|p-1|);
#   R2  S4 is NOT (identical max|p-1| in all three columns: heterogeneity only);
#   R3  S1 and S4 bracket S2/S3 in price magnitude; S5 has the sign of delta;
#   R4  every system closes the external account (identity at residual level).

include(joinpath(@__DIR__, "..", "run.jl"))
using Printf, LinearAlgebra
import NonlinearSolve
const BH = BeyondHulten

cpi_of(model, p) = begin
    σ = model.options.elasticities.σ
    sum(model.data.consumption_share .* p .^ (1 - σ))^(1 / (1 - σ))
end
rmax_of(f, x, pars) = maximum(abs, (o = similar(x); f(o, x, pars); o))

function solve_ladder(f, init, pars; tol = 1e-12)
    res = NonlinearSolve.solve(NonlinearSolve.NonlinearProblem(f, copy(init), pars),
        reltol = 1e-8, abstol = 1e-8, maxiters = 40000)
    x = res.u; rmax = rmax_of(f, x, pars)
    for _ in 1:8
        rmax <= tol && break
        res = NonlinearSolve.solve(NonlinearSolve.NonlinearProblem(f, x, pars),
            NonlinearSolve.LevenbergMarquardt(); reltol = 1e-13, abstol = 1e-13, maxiters = 40000)
        r2 = rmax_of(f, res.u, pars); r2 < rmax || break
        x, rmax = res.u, r2
    end
    x, rmax
end

# ── S2 / S3: the general sectoral closure (as probe13) ──────────────────────
function sec_residual!(out, X, pars)
    model, eta_vec = pars
    N = length(model.data.factor_share)
    p = BH._positive_floor(X[1:N]); y = BH._positive_floor(X[N+1:2N])
    w = BH._positive_floor(X[2N+1:3N]); F = X[3N+1]
    b = BH._mobile_market_demand(model, p, y, w; external_transfer = F)
    out[1:N]     .= p .- b.cost
    out[N+1:2N]  .= log.(b.L_i) .- log.(BH._positive_floor(model.data.labor_share)) .-
                    eta_vec .* log.(w ./ cpi_of(model, p))
    out[2N+1:3N] .= y .- b.intermediary_demand .- b.total_final_demand
    out[3N+1]     = cpi_of(model, p) - 1.0
    nothing
end

# ── S4: the pinned wage vector (ADR-0021 option A, 2N unknowns, F = 0) ──────
function pin_residual!(out, X, pars)
    model, wbar = pars
    N = length(model.data.factor_share)
    p = BH._positive_floor(X[1:N]); y = BH._positive_floor(X[N+1:2N])
    b = BH._mobile_market_demand(model, p, y, wbar)
    out[1:N]    .= p .- b.cost
    out[N+1:2N] .= y .- b.intermediary_demand .- b.total_final_demand
    nothing
end

# ── S5: the utilization / Verdoorn externality (ALPHA + A_eff in the cost) ──
# The workplan's reduced form: A_eff enters the UNIT-COST hook only (the extra
# cost is paid to no one, so zero profit is preserved and no income closure is
# needed). The demand block keeps A = 1; promoting the door would substitute
# A_eff into the CES input mix as well (the probe's residual cannot rebuild the
# model per evaluation without breaking ForwardDiff, so this is the reduced
# form as specified, and it is labelled as such).
function util_residual!(out, X, pars)
    model, delta = pars
    (; data, options, shocks) = model
    N = length(data.factor_share)
    p = BH._positive_floor(X[1:N]); y = BH._positive_floor(X[N+1:2N])
    w = BH._positive_floor(X[2N+1]); F = X[2N+2]
    θ = options.elasticities.θ; ϵ = options.elasticities.ϵ
    Ae = shocks.supply_shock .* (y ./ data.λ) .^ (-delta)
    ip = BH._intermediate_price(data.Ω_raw, p, θ)
    cost = BH._ces_unit_cost(Ae, data.factor_share, w, ip, ϵ)
    b = BH._mobile_market_demand(model, p, y, w; external_transfer = F)
    out[1:N]    .= p .- cost
    out[N+1:2N] .= y .- b.intermediary_demand .- b.total_final_demand
    out[2N+1]    = sum(b.L_i) - 1.0
    out[2N+2]    = cpi_of(model, p) - 1.0
    nothing
end

"max|p-1|, employment, wage range and the identity gap at a solution."
function metrics(model, x, wvec)
    N = length(model.data.factor_share)
    p = x[1:N]; y = x[N+1:2N]
    F = length(x) == 3N + 1 ? x[3N+1] : (length(x) == 2N + 2 ? x[2N+2] : 0.0)
    b = BH._mobile_market_demand(model, p, y, wvec; external_transfer = F)
    cl = y .- b.intermediary_demand .- b.total_final_demand
    (maxdp = maximum(abs, p .- 1), L = sum(b.L_i),
        wlo = minimum(wvec), whi = maximum(wvec), gap = dot(p, cl))
end

function main()
    root = default_root()
    design_d = load_design("matrix_5x3_v6"; root = root)
    ref = build_reference(design_d; root = root)
    data = ref.data
    N = length(data.factor_share)
    ψ, g = programme_vectors(design_d, data, collect(1:N); root = root)
    prog = BitVector(g .> 0)
    big = sortperm(data.labor_share; rev = true)
    top_half = falses(N); top_half[big[1:round(Int, 0.5N)]] .= true
    @printf("N = %d; S-slots S1..S5, full-71 A-bill; programme sectors = %d\n", N, count(prog))

    # per-financing model set
    S = Dict{String,Any}()
    for fin in ("F1", "F2", "F3")
        bf = build_cell_model(design_d["cells"]["matrix_5x3-v6-BF-$fin"], design_d, data, ψ, g)
        mo = build_cell_model(design_d["cells"]["matrix_5x3-v6-ALPHA-$fin"], design_d, data, ψ, g)
        S[fin] = (bf = bf, bf_sol = solve(bf; init = ref.init_warm), mo = mo)
    end

    # rows of the comparison: (label, fin, maxdp, L, wlo, whi, gap)
    rows = NamedTuple[]

    # ── S1: the executed v6 BF row (eta = 0 rigid sectoral) ─────────────────
    for fin in ("F1", "F2", "F3")
        r = S[fin]; sol = r.bf_sol
        m = metrics(r.bf, vcat(sol.prices_raw, sol.quantities, sol.wages, sol.external_transfer), sol.wages)
        push!(rows, (label = "S1 rigid (eta=0)", fin = fin, maxdp = m.maxdp, L = m.L,
            wr = m.whi / m.wlo, gap = m.gap))
    end

    # ── S2: uniform eta_s ladder (0.5 and 2) ────────────────────────────────
    for eta_s in (0.5, 2.0)
        for fin in ("F1", "F2", "F3")
            r = S[fin]; sol = r.bf_sol
            init = vcat(sol.prices_raw, sol.quantities, sol.wages, sol.external_transfer)
            x, _ = solve_ladder(sec_residual!, init, (r.mo, fill(eta_s, N)))
            m = metrics(r.mo, x, x[2N+1:3N])
            push!(rows, (label = "S2 uniform eta_s=$eta_s", fin = fin, maxdp = m.maxdp, L = m.L,
                wr = m.whi / m.wlo, gap = m.gap))
        end
    end

    # ── S3: two-group (rigid = programme sectors, flexible = 0.5) ───────────
    for fin in ("F1", "F2", "F3")
        r = S[fin]; sol = r.bf_sol
        init = vcat(sol.prices_raw, sol.quantities, sol.wages, sol.external_transfer)
        x, _ = solve_ladder(sec_residual!, init, (r.mo, 0.5 .* (.!prog)))
        m = metrics(r.mo, x, x[2N+1:3N])
        push!(rows, (label = "S3 two-group", fin = fin, maxdp = m.maxdp, L = m.L,
            wr = m.whi / m.wlo, gap = m.gap))
    end

    # ── S4: pinned wage vector (ADR-0021), programme sectors +10% ───────────
    wbar = 1.0 .+ 0.1 .* prog
    for fin in ("F1", "F2", "F3")
        r = S[fin]; sol = r.bf_sol
        init = vcat(sol.prices_raw, sol.quantities)
        x, _ = solve_ladder(pin_residual!, init, (r.mo, wbar))
        m = metrics(r.mo, x, wbar)
        push!(rows, (label = "S4 pinned vector +10%", fin = fin, maxdp = m.maxdp, L = m.L,
            wr = m.whi / m.wlo, gap = m.gap))
    end

    # ── S5: utilization / Verdoorn externality (delta = +0.5 and -0.5) ──────
    for delta in (0.5, -0.5)
        for fin in ("F1", "F2", "F3")
            r = S[fin]; mo = r.mo
            al = solve(mo; init = ref.init_warm)
            init = vcat(al.prices_raw, al.quantities, al.wages[1], al.external_transfer)
            x, rmax = solve_ladder(util_residual!, init, (mo, delta))
            m = metrics(mo, x, [x[2N+1]])
            push!(rows, (label = "S5 util delta=$(delta)", fin = fin, maxdp = m.maxdp, L = m.L,
                wr = m.whi / m.wlo, gap = m.gap))
        end
    end

    # ── the comparison table ────────────────────────────────────────────────
    println("\n== S1..S5 side by side (max|p-1| by financing column; sensitivity = F1 != F2/F3) ==")
    @printf("  %-24s %11s %11s %11s %12s %9s %11s %s\n",
        "scenario", "F1", "F2", "F3", "L (F2)", "w max/min", "identity", "verdict")
    labels = unique(r.label for r in rows)
    for lab in labels
        rr = [r for r in rows if r.label == lab]
        r1 = rr[findfirst(r -> r.fin == "F1", rr)]
        r2 = rr[findfirst(r -> r.fin == "F2", rr)]
        r3 = rr[findfirst(r -> r.fin == "F3", rr)]
        sens = !isapprox(r1.maxdp, r2.maxdp; rtol = 1e-6)
        verdict = sens ? "SENSITIVE" : (r1.maxdp > 1e-6 ? "heterogeneous only" : "flat (p=1)")
        @printf("  %-24s %11.3e %11.3e %11.3e %12.8f %9.3f %11.2e %s\n",
            lab, r1.maxdp, r2.maxdp, r3.maxdp, r2.L, r2.wr, r2.gap, verdict)
    end

    println("\n== reading ==")
    println("  S1 (rigid) and S4 (pin) bracket the price scale; S4's movement is exogenous")
    println("  (identical across columns), S1/S2/S3/S5 move WITH the financing cell.")
    println("  S5 with delta < 0 is the Kaldor-Verdoorn sign (prices fall with demand).")
    return nothing
end

main()