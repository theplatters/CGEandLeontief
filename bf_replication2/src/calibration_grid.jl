# src/calibration_grid.jl — Phase 5 driver: calibration grid + CSV output
#
# Replicates the four nested loops of Master_file_3.m (lines 226–505):
#   loop (1: no htm sweep, 2: htm sweep) × elasticity (1: benchmark, 2: CD)
#   × shock_type (1:5) × htm_share s × t_grid
#
# For each cell it solves the equilibrium along the t-grid (continuation) and
# stores the MATLAB-comparable outputs. Results are written as CSVs so they can
# be consumed outside the solver environment (e.g. on a memory-limited container).
#
# Usage (from bf_replication2/):
#   julia --project=. -e 'include("src/calibration_grid.jl"); run_calibration_grid("data/results")'
#
# Outputs per (loop, elasticity, shock_type, s):
#   results/loop{loop}/el{elasticity}/st{shock_type}/s{s}/timeseries.csv
#   plus summary files results/partA_C.csv, results/partF_htm.csv, results/baseline_fit.csv

using Printf, LinearAlgebra, Statistics, CSV, DataFrames

export solve_cell, run_calibration_grid, default_t_grid

# ─────────────────────────────────────────────────────────────
# Shock grids (Master_file_3.m lines 253–265)
# ─────────────────────────────────────────────────────────────
default_t_grid(elasticity::Int, loop::Int=1) = elasticity == 1 ?
    (loop == 1 ? vcat([0.0], collect(0.1:0.05:0.9), [0.95, 0.98, 0.99, 0.995, 1.0]) :
                 [0.0, 0.001, 0.01, 0.5, 0.75, 0.98, 0.99, 0.995, 1.0]) :
    vcat([0.0], collect(0.1:0.1:0.9), [0.99, 1.0])

# Elasticity calibration (Master_file_3.m lines 253–275)
function elasticity_params(elasticity::Int)
    if elasticity == 1
        return (; sigma=1.0, epsilon=0.6, eta=0.6, theta1=0.2, benchmark=true)
    else
        return (; sigma=1.0, epsilon=1.0, eta=1.0, theta1=1.0, benchmark=false)
    end
end

# ─────────────────────────────────────────────────────────────
# Shock application (Master_file_3.m lines 293–395)
# Returns (A, B) for a given t and shock_type. Trunc_A is the persistent
# (N × n_t) matrix from the previous shock_type cell (needed for st=2).
# ─────────────────────────────────────────────────────────────
function apply_shocks(m_base::MCPModel, shocks::ShockData, t::Float64,
                      ti::Int, shock_type::Int, trunc_A)
    N = m_base.N
    D = m_base.D
    A = ones(D)
    B = ones(D, D)

    if shock_type == 1        # supply + sectoral demand + aggregate demand
        shock = -shocks.BLS_shock
        A[(3*N+2):(4*N+1)] .= (1 .- t .* shock)
        B[1, 2:(N+1)] .= (1 - t*0.66) .+ t*0.66 .* (1.0 .+ shocks.PCE_shock)
        B[5*N+3, D] = 1.0 + min(ti - 1, 1) * 0.105
    elseif shock_type == 2    # supply only (uses Trunc_A from previous cell)
        shock = 1.0 .- trunc_A[:, end]
        A[(3*N+2):(4*N+1)] .= (1 .- t .* shock)
    elseif shock_type == 3    # sectoral + aggregate demand
        B[1, 2:(N+1)] .= (1 - t*0.66) .+ t*0.66 .* (1.0 .+ shocks.PCE_shock)
        B[5*N+3, D] = 1.0 + min(ti - 1, 1) * 0.105 * t
    elseif shock_type == 4    # aggregate demand only
        B[5*N+3, D] = 1.0 + min(ti - 1, 1) * 0.105
    elseif shock_type == 5    # supply + sectoral demand
        shock = -shocks.BLS_shock
        A[(3*N+2):(4*N+1)] .= (1 .- t .* shock)
        B[1, 2:(N+1)] .= (1 - t*0.66) .+ t*0.66 .* (1.0 .+ shocks.PCE_shock)
    else
        error("Unknown shock_type = $shock_type")
    end

    # Row renormalisation (Master_file_3.m lines 394–395)
    # NOTE: no transpose — Omega_re[1,:] .* B[1,:] must stay a dot product.
    B[1, :] ./= sum(m_base.Omega_re[1, :] .* B[1, :])
    B[5*N+3, :] ./= sum(m_base.Omega_re[5*N+3, :] .* B[5*N+3, :])

    return A, B
end

# ─────────────────────────────────────────────────────────────
# Solve one (elasticity, shock_type, s) cell along the t-grid
# ─────────────────────────────────────────────────────────────
"""
    solve_cell(io, shocks, sf; elasticity=1, shock_type=1, htm_share=0.0,
               t_grid=default_t_grid(elasticity), trunc_A=nothing) -> (res, trunc_A)

Solve the equilibrium along `t_grid` for one calibration cell, replicating the
MATLAB driver's inner loop. Returns a Dict with time series:
  "t", "GDP" (real, λ1/p1), "nominal_GDP" (λ1), "Hulten", "p1" (inflation),
  "unemployment1", "unemployment2", "prices" (N×n_t goods prices),
  "c_shares" (N×n_t consumption shares), "lambda" (D×n_t), "A_labor" (N×n_t),
  and the updated trunc_A matrix (N×n_t).

`trunc_A` is the persistent matrix from the previous shock_type cell (used by
shock_type=2). Pass `nothing` to initialise.
"""
function solve_cell(io, shocks, sf;
                    elasticity::Int=1, shock_type::Int=1,
                    htm_share::Float64=0.0,
                    t_grid::Vector{Float64}=default_t_grid(1),
                    trunc_A=nothing)
    N = sf.N
    D = sf.D
    n_t = length(t_grid)

    el = elasticity_params(elasticity)
    m_base = make_model(sf; htm_share=htm_share, benchmark=el.benchmark,
                        sigma=el.sigma, epsilon=el.epsilon, eta=el.eta, theta1=el.theta1)

    Psi_re = inv(I - m_base.Omega_re)

    # Storage (MATLAB names)
    GDP          = zeros(n_t)
    nominal_GDP  = zeros(n_t)
    Hulten       = zeros(n_t)
    inflation    = zeros(n_t)          # p[1]
    unemp1       = zeros(n_t)          # unemployment1
    unemp2       = zeros(n_t)          # unemployment2
    prices       = zeros(N, n_t)       # goods prices (rows 2:N+1)
    c_shares     = zeros(N, n_t)       # consumption shares
    lambda_full  = zeros(D, n_t)
    A_labor      = zeros(N, n_t)       # labor A rows (3N+2:4N+1)
    # Initialize Trunc_A: preserve existing if dimensions match, otherwise create new
    if trunc_A === nothing
        Trunc_A = zeros(N, n_t)
        Trunc_A .= 1.0  # Initialize all to 1 (not demand-constrained)
    elseif size(trunc_A, 1) == N && size(trunc_A, 2) >= n_t
        # Preserve existing values, but ensure we have at least n_t columns
        if size(trunc_A, 2) < n_t
            # Extend with ones if needed
            Trunc_A_extended = zeros(N, n_t)
            Trunc_A_extended[:, 1:size(trunc_A, 2)] .= trunc_A
            Trunc_A_extended[:, size(trunc_A, 2)+1:end] .= 1.0
            Trunc_A = Trunc_A_extended
        else
            Trunc_A = trunc_A[:, 1:n_t]  # Use only the needed columns
        end
    else
        @warn "Trunc_A dimension mismatch (expected $N×$n_t, got $(size(trunc_A))), reinitializing" maxlog=1
        Trunc_A = zeros(N, n_t)
        Trunc_A .= 1.0
    end
    retcodes     = zeros(Int, n_t)

    # Domar weights for labor (used in unemployment measures)
    labor_weights = Psi_re[1, (3*N+2):(4*N+1)]
    labor_weight_sum = sum(labor_weights)
    sector_unemp_last = zeros(N)   # sector-level unemployment at last t

    z0 = vcat(m_base.init_p, m_base.init_lambda)

    for (ti, t) in enumerate(t_grid)
        # shock_type=2 needs Trunc_A from the previous cell
        if shock_type == 2 && trunc_A === nothing
            error("shock_type=2 requires trunc_A from a prior shock_type=1 cell")
        end
        A_t, B_t = apply_shocks(m_base, shocks, t, ti, shock_type, Trunc_A)

        m_t = MCPModel(D, N, m_base.Omega_re, m_base.factor, m_base.keynes,
                       m_base.theta, m_base.cobb_douglas, m_base.phi, m_base.phi_htm,
                       m_base.mu, m_base.in_mu, A_t, B_t,
                       m_base.init_lambda, m_base.init_p, m_base.chi)

        p = fill(NaN, D)
        λ = fill(NaN, D)
        conv = false
        iters = 0
        try
            p, λ, conv, iters = solve_equilibrium(m_t, z0=z0, tol=1e-10, maxiter=1000)
        catch
            # Solver threw an exception (e.g. NaN evaluation) — will try refinement below
        end

        # If solver fails, try multiple fallback strategies
        if !conv && ti > 1
            t_prev = t_grid[ti - 1]
            t_mid = (t_prev + t) / 2
            
            # Strategy 1: Simple continuation refinement
            A_mid, B_mid = apply_shocks(m_base, shocks, t_mid, ti, shock_type, Trunc_A)
            m_mid = MCPModel(D, N, m_base.Omega_re, m_base.factor, m_base.keynes,
                             m_base.theta, m_base.cobb_douglas, m_base.phi, m_base.phi_htm,
                             m_base.mu, m_base.in_mu, A_mid, B_mid,
                             m_base.init_lambda, m_base.init_p, m_base.chi)
            try
                p_mid, λ_mid, _, _ = solve_equilibrium(m_mid, z0=z0, tol=1e-8, maxiter=500)
                p, λ, conv, iters = solve_equilibrium(m_t, z0=vcat(p_mid, λ_mid), tol=1e-10, maxiter=1000)
            catch e
                @warn "Continuation refinement failed at t=$t (htm_share=$(m_base.phi_htm[3*m_base.N+2]): $(sprint(showerror, e))" maxlog=1
                
                # Strategy 2: Relax tolerance and try again
                try
                    p, λ, conv, iters = solve_equilibrium(m_t, z0=z0, tol=1e-8, maxiter=1000)
                catch e2
                    @warn "Second attempt also failed at t=$t: $(sprint(showerror, e2))" maxlog=1
                    
                    # Strategy 3: Use the previous solution as fallback
                    p = z0[1:D]
                    λ = z0[D+1:2D]
                    conv = false
                    iters = 0
                end
            end
        end

        # Store
        GDP[ti]         = λ[1] / p[1]
        nominal_GDP[ti] = λ[1]
        Hulten[ti]      = exp(sum(Psi_re[1, :] .* log.(A_t)))
        inflation[ti]   = p[1]
        prices[:, ti]   = p[2:(N+1)]
        lambda_full[:, ti] = λ

        # Consumption shares: cons_share[k] for goods (AMPL consumption_share)
        # cons_share_k = B[1,k]^θ1 · Ω[1,k] · (p_k/p_1)^(1-θ1) · A_1^(θ1-1)
        θ1 = m_base.theta[1]
        for k in 1:N
            j = k + 1
            c_shares[k, ti] = (B_t[1, j]^θ1) * m_base.Omega_re[1, j] *
                              (p[j] / p[1])^(1 - θ1) * A_t[1]^(θ1 - 1)
        end

        # Trunc_A: A(labor), set to 1 where price < 1.001 (demand-constrained).
        # MATLAB: Trunc_A(init_p(labor)<1.001,:)=1 — but this is a NEW column per
        # counter, so effectively column-wise; the "end" column used by st=2 is
        # the PREVIOUS cell's final column, which must survive this cell's writes.
        A_labor[:, ti] = A_t[(3*N+2):(4*N+1)]
        Trunc_A[:, ti] = A_labor[:, ti]
        constrained = p[(3*N+2):(4*N+1)] .< 1.001
        if any(constrained)
            Trunc_A[constrained, ti] .= 1.0
        end

        # Unemployment measures (Master_file_3.m lines 483–491)
        lab_λ = λ[(3*N+2):(4*N+1)]
        lab_p = p[(3*N+2):(4*N+1)]
        # Sector-level unemployment: ((λ/p) - Psi_re_labor) ./ Psi_re_labor  (line 546)
        sector_unemp = (lab_λ ./ lab_p .- labor_weights) ./ labor_weights
        # unempl1: weighted log-deviation using Trunc_A (line 484)
        unemp1[ti] = sum((labor_weights ./ labor_weight_sum) .*
            (log.(lab_λ ./ lab_p) .- log.(labor_weights .* Trunc_A[:, ti])))
        # unemployment2: without Trunc_A, using actual A; healthcare (58,59) zeroed unless st=2
        unemp2_temp = (labor_weights ./ labor_weight_sum) .*
            (log.(lab_λ ./ lab_p) .- log.(labor_weights .* A_labor[:, ti]))
        if shock_type != 2
            unemp2_temp[58:59] .= 0.0
        end
        unemp2[ti] = sum(unemp2_temp)
        # Store sector-level unemployment for the last t (baseline fit)
        if ti == n_t
            sector_unemp_last .= sector_unemp
        end

        retcodes[ti] = conv ? 0 : 1
        
        # Clamp values to prevent NaN propagation in subsequent iterations
        if isnan(GDP[ti]) || isinf(GDP[ti])
            GDP[ti] = GDP[ti > 1 ? ti-1 : 1]  # Use previous or first value as fallback
            @warn "NaN detected in GDP at t=$(t_grid[ti]), using previous value" maxlog=1
        end
        if isnan(nominal_GDP[ti]) || isinf(nominal_GDP[ti])
            nominal_GDP[ti] = nominal_GDP[ti > 1 ? ti-1 : 1]
            @warn "NaN detected in nominal GDP at t=$(t_grid[ti]), using previous value" maxlog=1
        end
        if isnan(inflation[ti]) || isinf(inflation[ti])
            inflation[ti] = inflation[ti > 1 ? ti-1 : 1]
            @warn "NaN detected in inflation at t=$(t_grid[ti]), using previous value" maxlog=1
        end
        
        # Clamp prices and lambda to prevent NaN in subsequent iterations
        p_clamped = [max(1e-15, p[i]) for i in 1:D]
        λ_clamped = [max(1e-15, λ[i]) for i in 1:D]
        z0 = vcat(p_clamped, λ_clamped)   # continuation initializer
    end

    return Dict(
        "t" => t_grid, "GDP" => GDP, "nominal_GDP" => nominal_GDP,
        "Hulten" => Hulten, "p1" => inflation,
        "unemployment1" => unemp1, "unemployment2" => unemp2,
        "sector_unemp" => sector_unemp_last,
        "prices" => prices, "c_shares" => c_shares,
        "lambda" => lambda_full, "A_labor" => A_labor,
        "labor_weights" => labor_weights,
        "Trunc_A" => Trunc_A, "retcodes" => retcodes
    ), Trunc_A
end

# ─────────────────────────────────────────────────────────────
# Chain-weighted real GDP change (Master_file_3.m line 535)
# ─────────────────────────────────────────────────────────────
function delta_r_gdp(res::Dict)
    nom = res["nominal_GDP"]
    cs  = res["c_shares"]
    pr  = res["prices"]
    n_t = length(nom)
    dlog_p = diff(log.(pr), dims=2)            # (N × n_t-1)
    cs_avg = (cs[:, 1:end-1] .+ cs[:, 2:end]) ./ 2
    contrib = sum(cs_avg .* dlog_p, dims=1)    # (1 × n_t-1)
    dlog_nom = diff(log.(nom))
    return 1 .- exp.(cumsum(dlog_nom .- vec(contrib)))
end

# ─────────────────────────────────────────────────────────────
# Full calibration grid (Master_file_3.m loops 226–505)
# ─────────────────────────────────────────────────────────────
"""
    run_calibration_grid(outdir="data/results"; loops=[1,2], verbose=true)

Run the full calibration grid (loop × elasticity × shock_type × htm_share)
and write CSVs to `outdir`. Returns a Dict with the summary matrices
(RGDP_graph, Inflation_graph, Unemp_graph, and the loop-2 htm variants).

The grid is memory-heavy (~25 solves × 5 shock types × 2 elasticities × 6 htm
shares for loop 2) — run on the host Mac (32 GB). Each cell's CSV is written
immediately, so partial results survive.
"""
function run_calibration_grid(outdir::String="data/results"; loops=[1, 2], verbose=true)
    DATA_DIR = "data"
    io = load_io_table(joinpath(DATA_DIR, "IO_data_2018.mat"); N=66, year=2015)
    shocks = load_shocks(DATA_DIR; N=66)
    sf = build_standard_form(io)

    mkpath(outdir)

    RGDP_graph = zeros(2, 5)
    Inflation_graph = zeros(2, 5)
    Unemp_graph = zeros(2, 5)
    RGDP_graph_htm = zeros(2, 6)
    Inflation_graph_htm = zeros(2, 6)
    Unemp_graph_htm = zeros(2, 6)

    baseline = nothing
    # Trunc_A persists across ALL cells in MATLAB (never cleared)
    # Initialize with maximum possible size across all elements
    max_n_t = 0
    for loop in loops
        shock_loop = loop == 1 ? (1:5) : (1:1)
        for elasticity in 1:2
            t_grid = default_t_grid(elasticity, loop)
            max_n_t = max(max_n_t, length(t_grid))
        end
    end
    trunc_A = zeros(N, max_n_t)
    trunc_A .= 1.0  # Initialize all to 1 (not demand-constrained)

    for loop in loops
        shock_loop = loop == 1 ? (1:5) : (1:1)
        htm_loop   = loop == 1 ? (1:1) : (1:6)

        for elasticity in 1:2
            for shock_type in shock_loop
                for s in htm_loop
                    htm_share = 0.2 * (s - 1)
                    t_grid = default_t_grid(elasticity, loop)

                    res, trunc_A = solve_cell(io, shocks, sf;
                                              elasticity=elasticity,
                                              shock_type=shock_type,
                                              htm_share=htm_share,
                                              t_grid=t_grid,
                                              trunc_A=trunc_A)

                    # Write cell CSV immediately
                    cell_dir = joinpath(outdir, "loop$loop", "el$elasticity",
                                        "st$shock_type", "s$s")
                    mkpath(cell_dir)
                    write_cell_csv(cell_dir, res, io)

                    dRGDP = delta_r_gdp(res)
                    nom_end = res["nominal_GDP"][end]
                    unemp1_end = res["unemployment1"][end]

                    if loop == 1
                        RGDP_graph[elasticity, shock_type]     = round(-100 * dRGDP[end], digits=2)
                        Inflation_graph[elasticity, shock_type] = round(100 * (nom_end / (1 - dRGDP[end]) - 1), digits=2)
                        Unemp_graph[elasticity, shock_type]     = round(-100 * (exp(unemp1_end) - 1), digits=2)
                    else
                        RGDP_graph_htm[elasticity, s]     = -dRGDP[end]
                        Inflation_graph_htm[elasticity, s] = nom_end / (1 - dRGDP[end]) - 1
                        Unemp_graph_htm[elasticity, s]     = -(exp(unemp1_end) - 1)
                    end

                    # Baseline (loop=1, el=1, st=1) for Parts A/B/D
                    if loop == 1 && elasticity == 1 && shock_type == 1 && s == 1
                        baseline = res
                    end

                    if verbose
                        @printf("loop=%d el=%d st=%d s=%d (htm=%.1f) | RGDP Δ=%.2f%% | n=%d pts | retcodes=%s\n",
                                loop, elasticity, shock_type, s, htm_share,
                                -100 * dRGDP[end], length(t_grid), res["retcodes"])
                    end
                end
            end
        end
    end

    # Summary CSVs
    write_summary_csv(outdir, RGDP_graph, Inflation_graph, Unemp_graph,
                      RGDP_graph_htm, Inflation_graph_htm, Unemp_graph_htm)

    if baseline !== nothing && 1 in loops
        write_baseline_fit(outdir, baseline, io, shocks)
    end

    return Dict("RGDP_graph" => RGDP_graph,
                "Inflation_graph" => Inflation_graph,
                "Unemp_graph" => Unemp_graph,
                "RGDP_graph_htm" => RGDP_graph_htm,
                "Inflation_graph_htm" => Inflation_graph_htm,
                "Unemp_graph_htm" => Unemp_graph_htm)
end

# ─────────────────────────────────────────────────────────────
# CSV writers
# ─────────────────────────────────────────────────────────────
function write_cell_csv(cell_dir::String, res::Dict, io)
    n_t = length(res["t"])
    N = io.N

    # Long-format time series
    df = DataFrame(
        t = res["t"],
        GDP = res["GDP"],
        nominal_GDP = res["nominal_GDP"],
        Hulten = res["Hulten"],
        p1 = res["p1"],
        unemployment1 = res["unemployment1"],
        unemployment2 = res["unemployment2"],
        retcode = res["retcodes"],
    )
    CSV.write(joinpath(cell_dir, "timeseries.csv"), df)

    # Sector-level matrices (long format)
    # prices, c_shares, A_labor for each (sector, t)
    rows = []
    for k in 1:N, ti in 1:n_t
        push!(rows, (sector=k, name=io.indname[k], t=res["t"][ti],
                     price=res["prices"][k, ti],
                     c_share=res["c_shares"][k, ti],
                     A_labor=res["A_labor"][k, ti]))
    end
    CSV.write(joinpath(cell_dir, "sector_prices.csv"), DataFrame(rows))
end

function write_summary_csv(outdir::String, RGDP, InfG, UnempG, RGDP_htm, InfG_htm, UnempG_htm)
    shock_names = ["baseline", "supply_only", "demand_only", "agg_demand_only", "supply_sectoral"]

    # loop=1: elasticity × shock_type matrices
    df_el1 = DataFrame(shock_type = shock_names,
                       RGDP_benchmark = RGDP[1, :],
                       Inflation_benchmark = InfG[1, :],
                       Unemp_benchmark = UnempG[1, :],
                       RGDP_cd = RGDP[2, :],
                       Inflation_cd = InfG[2, :],
                       Unemp_cd = UnempG[2, :])
    CSV.write(joinpath(outdir, "summary_loop1.csv"), df_el1)

    # loop=2: htm sweep (elasticity × s)
    htm_shares = 0.0:0.2:1.0
    df_htm = DataFrame(htm_share = collect(htm_shares),
                       RGDP_benchmark = RGDP_htm[1, :],
                       Inflation_benchmark = InfG_htm[1, :],
                       Unemp_benchmark = UnempG_htm[1, :],
                       RGDP_cd = RGDP_htm[2, :],
                       Inflation_cd = InfG_htm[2, :],
                       Unemp_cd = UnempG_htm[2, :])
    CSV.write(joinpath(outdir, "summary_loop2_htm.csv"), df_htm)
end

function write_baseline_fit(outdir::String, res::Dict, io, shocks)
    N = io.N
    sticky = res["Trunc_A"][:, end] .== 1.0
    sector_unemp = res["sector_unemp"]
    df = DataFrame(sector = 1:N,
                   name = io.indname[1:N],
                   sticky = sticky,
                   price_model = res["prices"][:, end] .- 1,  # model-implied inflation
                   unemp_model = sector_unemp,
                   PPI_data = shocks.PPI,
                   BLS_shock = shocks.BLS_shock,
                   PCE_shock = shocks.PCE_shock,
                   wages_data = shocks.wages)
    CSV.write(joinpath(outdir, "baseline_fit.csv"), df)
end
