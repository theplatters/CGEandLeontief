# tests/test_pricing_design.jl — the Phase 7 pricing-door harness (ADR-0029,
# ADR-0030) on the paper side: the `[pricing]` design schema, its
# preregistration pinning, and the real-table contract the synthetic kernel
# fixtures cannot see. The kernel's own contract tests are
# OpenCGE@tests/test_pricing_closure.jl.
#
# The real-table testsets need data/I-O_DE2019_formatiert.csv and
# cbase2/data_raw/impulses.csv; they skip cleanly without them.

using OpenCGE, Test, LinearAlgebra, TOML
using CSV, DataFrames
isdefined(Main, :run_design) || include(joinpath(@__DIR__, "..", "experiments", "run.jl"))
isdefined(Main, :ext_test_data) || include(joinpath(@__DIR__, "test_external_closure.jl"))

@testset "pricing design: cell_pricing resolves the [pricing] schema" begin
    # Absent block: NoPricing and no manifest fields (byte-identical manifests).
    pr, f = cell_pricing(Dict{String,Any}(), Dict{String,Any}())
    @test pr isa NoPricing
    @test isempty(f)
    # Design defaults, overridden per cell.
    d = Dict{String,Any}("pricing" => Dict{String,Any}("kind" => "capacity", "delta" => 0.0))
    pr, f = cell_pricing(d, Dict{String,Any}("delta" => -0.5))
    @test pr isa CapacityPricing && pr.δ == -0.5
    @test f["pricing"] == "P1" && f["delta"] == -0.5
    # Markup: s_pi defaults to 1 (Stage A, retained).
    d = Dict{String,Any}("pricing" => Dict{String,Any}("kind" => "markup", "kappa" => 0.25))
    pr, f = cell_pricing(d, Dict{String,Any}())
    @test pr isa MarkupPricing && pr.κ == 0.25 && pr.s_pi == 1.0
    @test f["pricing"] == "P2" && f["s_pi"] == 1.0
    pr, _ = cell_pricing(d, Dict{String,Any}("s_pi" => 0.5))
    @test pr.s_pi == 0.5
    # Per-cell kind override (the mobile robustness design).
    d = Dict{String,Any}("pricing" => Dict{String,Any}("kind" => "none"))
    pr, _ = cell_pricing(d, Dict{String,Any}("kind" => "markup", "kappa" => 0.5))
    @test pr isa MarkupPricing && pr.κ == 0.5
    @test_throws ArgumentError cell_pricing(
        Dict{String,Any}("pricing" => Dict{String,Any}("kind" => "import")), Dict{String,Any}())
    # The continuation path scales the door parameter only.
    @test scaled_pricing(CapacityPricing(-0.5), 0.5).δ == -0.25
    s = scaled_pricing(MarkupPricing(0.5; s_pi = 0.3), 0.25)
    @test s.κ == 0.125 && s.s_pi == 0.3
    @test scaled_pricing(NoPricing(), 0.5) isa NoPricing
    @test scaled_pricing(CapacityPricing(-0.5), 1.0).δ == -0.5
end

@testset "pricing design: the four Phase 7 designs parse and pin their door" begin
    root = normpath(joinpath(@__DIR__, ".."))
    kinds = Dict("pricing_capacity_v1" => CapacityPricing, "pricing_markup_v1" => MarkupPricing,
        "pricing_profit_v1" => MarkupPricing)
    for design in ("pricing_capacity_v1", "pricing_markup_v1", "pricing_profit_v1", "pricing_mobile_v1")
        d = load_design(design; root = root)
        cells = cell_order(design; root = root)
        @test !isempty(cells)
        for id in cells
            @test startswith(id, replace(design, "_v1" => "-v1-"))
            pr, f = cell_pricing(d, design_cell(d, id))
            haskey(kinds, design) && @test pr isa kinds[design]
            @test f["pricing"] in ("P1", "P2")
        end
        pins = pinned_params(d)
        @test haskey(pins, "pricing_kind")
        @test pins["pricing_multistart_tol"] == 1e-6
    end
    # A pricing-free design keeps its pre-pricing preregistration pins.
    @test !any(startswith("pricing_"), keys(pinned_params(load_design("matrix_5x3_v10"; root = root))))
end

@testset "pricing design: the preregistration hash covers the door parameters" begin
    root = normpath(joinpath(@__DIR__, ".."))
    src = read(joinpath(root, "experiments", "designs", "pricing_capacity_v1.toml"), String)
    tmp = mktempdir()
    mkpath(joinpath(tmp, "experiments", "designs"))
    path = joinpath(tmp, "experiments", "designs", "pricing_capacity_v1.toml")
    write(path, src)
    h0 = design_sha256("pricing_capacity_v1"; root = tmp)
    write(path, replace(src, "delta = -0.5" => "delta = -0.45"; count = 1))
    @test design_sha256("pricing_capacity_v1"; root = tmp) != h0
end

"""Real-table GAMMA cell with a pricing door, or `nothing` without the data."""
function pricing_real_cell(fin_id::AbstractString, pricing::AbstractPricing)
    data = ext_test_data()
    data === nothing && return nothing
    prog = ext_matrix_programme(data)
    prog === nothing && return nothing
    g, ψ = prog
    fin = cell_financing(fin_id, ψ, g, "tilt_g0_over_c0", data)
    N = length(data.factor_share)
    mk(pr) = mobile_labor_model(data, Shocks(ones(N), ones(N), zeros(N)), 0.5, 0.5, 0.9, 1.0;
        closure = :fixed, financing = fin, pricing = pr)
    base = solve(mk(NoPricing()))
    x = OpenCGE._solution_vector(mk(NoPricing()), base)
    sol = base
    for k in 1:4  # the run.jl continuation path
        m = mk(scaled_pricing(pricing, k / 4))
        sol = solve(m; init = x)
        x = OpenCGE._solution_vector(m, sol)
    end
    return (; data, base, sol, X = x)
end

@testset "pricing design: real-table external account with the door open" begin
    # The bill-factor FOC sites (ADR-0029 #4, sites 3-4) value M_int and T_int,
    # which every synthetic fixture sets to zero: only the real table exercises
    # them. Each door must close the external account as written (with the
    # retained rent booked) and keep the companion identity gdp_wedge = -diff.
    data = ext_test_data()
    if data === nothing || ext_matrix_programme(data) === nothing
        @test_skip "real IO table or impulse file absent"
    else
        @test sum(data.M_int) > 0 && sum(data.T_int) > 0
        for (fin_id, pr) in (("F2", CapacityPricing(-0.5)), ("F3", MarkupPricing(0.5)),
                ("F1", MarkupPricing(0.5; s_pi = 0.1199)))
            c = pricing_real_cell(fin_id, pr)
            m = c.sol.model
            @test maximum(abs, equilibrium_residuals(m, c.X)) < 1e-9
            @test assert_external_account(m, c.sol) === nothing
            can = external_balance_canary(m, c.X)
            @test abs(can.diff) < 1e-9
            @test isapprox(gdp_wedge(c.sol), -can.diff; atol = 1e-9)
            # The door moves prices under a demand-only shock.
            @test maximum(abs.(c.sol.prices_raw .- 1)) > 1e-2
            if pr isa MarkupPricing
                # The retained rent is booked: dropping it leaves the identity
                # short by exactly Π_ret.
                @test can.retained_profit > 0
                @test isapprox(retained_profit(m, c.sol), can.retained_profit; rtol = 1e-12)
            else
                @test can.retained_profit == 0.0
            end
        end
    end
end

@testset "pricing design: evaluate_gates on a real-table priced cell" begin
    data = ext_test_data()
    if data === nothing || ext_matrix_programme(data) === nothing
        @test_skip "real IO table or impulse file absent"
    else
        # Stage B: the budget gate must include the distributed rent.
        c = pricing_real_cell("F2", MarkupPricing(0.5; s_pi = 0.5))
        design_d = Dict{String,Any}(
            "gates" => Dict{String,Any}("residual_tol" => 1e-6, "budget_tol" => 1e-9,
                "labour_tol" => 1e-6, "wage_tol" => 1e-8),
            "pricing" => Dict{String,Any}("kind" => "markup", "kappa" => 0.5, "s_pi" => 0.5,
                "multistart_tol" => 1e-6))
        cell = Dict{String,Any}("labor" => "GAMMA", "financing" => "F2")
        ev = evaluate_gates(cell, design_d, c.sol.model, c.sol, c.base)
        @test ev.gates["overall"] == "pass"
        @test ev.gates["budget"]["pass"]
        @test ev.gates["uniqueness"]["pass"]
        @test ev.metrics["profit"] > 0
        @test 0 < ev.metrics["profit_share"] < 1
        @test ev.metrics["retained_profit"] ≈ 0.5 * ev.metrics["profit"]
        @test ev.metrics["markup_min"] > 0
        @test ev.metrics["price_dev_min"] <= ev.metrics["price_dev_mean"] <= ev.metrics["price_dev_max"]
        @test ev.diagnostics["jacobian_rank"] == ev.diagnostics["jacobian_n"]
        # A pricing-free cell carries no pricing fields.
        ev0 = evaluate_gates(cell, Dict{String,Any}("gates" => design_d["gates"]),
            c.base.model, c.base, c.base)
        @test !haskey(ev0.gates, "uniqueness")
        @test !haskey(ev0.metrics, "profit")
    end
end
