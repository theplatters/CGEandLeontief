# experiments/probes/probe10_repin_eta0_fixtures.jl
#
# ADR-0020 promotion, step 3 support: measure the new η = 0 goldens on the
# test fixtures so the kernel-regression / closure tests can be re-pinned.
# The fixture constructors are copied verbatim from tests/test_kernel_regression.jl
# and tests/test_promoted_closures.jl (they are file-local there).

using BeyondHulten, Printf, DataFrames

function three_sector_fixture()
    io = DataFrame("Sektoren" => ["a", "b", "c"],
        "Letzte Verwendung von Gütern zusammen" => [1.0, 1.0, 1.0])
    Ω = [0.6 0.3 0.1; 0.2 0.5 0.3; 0.25 0.25 0.5]
    consumption_share = [0.5, 0.3, 0.2]
    factor_share = [0.6, 0.5, 0.4]
    λ = [1.2, 1.0, 0.8]
    labor_share = λ .* factor_share
    consumption_share_gross_output = [0.4, 0.35, 0.3]
    grossy = [10.0, 8.0, 6.0]
    value_added = factor_share .* grossy
    Data(io, Ω, consumption_share, factor_share, λ, labor_share,
        consumption_share_gross_output, grossy, value_added)
end

function v3_fixture()
    io = DataFrame("Sektoren" => ["a", "b", "c"],
        "Letzte Verwendung von Gütern zusammen" => [1.0, 1.0, 1.0])
    Ω = [0.6 0.3 0.1; 0.2 0.5 0.3; 0.25 0.25 0.5]
    factor_share = [0.6, 0.5, 0.4]
    grossy = [2.0, 1.5, 1.0]
    value_added = factor_share .* grossy
    gdp = sum(value_added)
    λ = grossy ./ gdp
    labor_share = λ .* factor_share
    m = 0.2
    gov = [0.03, 0.015, 0.005]
    inv = [0.03, 0.075, 0.045]
    expo = [0.075, 0.045, 0.03]
    c0_dom = λ .- Ω' * ((1 .- factor_share) .* λ) .-
        (1 - m) .* (gov .+ inv) .- expo
    c0_gross = c0_dom ./ (1 - m)
    saving_rate = 1 - sum(c0_gross) / (1 - sum(gov))
    ω = c0_gross ./ sum(c0_gross)
    data = Data(io, Ω, Ω, ω, factor_share, λ, labor_share, ω, grossy,
        value_added, grossy, DataFrame(), zeros(3), zeros(3), zeros(3),
        gov, c0_gross, fill(m, 3), inv, expo, saving_rate,
        (1 .- factor_share) .* λ, zeros(3), zeros(3), gdp, gdp, gdp)
    (; data = data, g = [0.02, 0.0, 0.0], shift = [1.5, 1.0, 0.8],
        saving_rate = saving_rate)
end

const _V3_θ, _V3_ϵ, _V3_σ = 1.0, 0.5, 0.9
_v3_shocks() = Shocks(ones(3), ones(3), zeros(3))

function tiny_fixture()
    io = DataFrame("Sektoren" => ["a", "b"], "Letzte Verwendung von Gütern zusammen" => [1.0, 1.0])
    Ω = [0.5 0.5; 0.5 0.5]
    consumption_share = [0.6, 0.4]
    factor_share = [0.7, 0.5]
    λ = [1.2, 0.8]
    labor_share = λ .* factor_share
    consumption_share_gross_output = [0.5, 0.5]
    grossy = [10.0, 8.0]
    value_added = factor_share .* grossy
    Data(io, Ω, consumption_share, factor_share, λ, labor_share,
        consumption_share_gross_output, grossy, value_added)
end

function show_cell(tag, model, sol)
    X = [sol.prices_raw; sol.quantities; sol.wages_raw; sol.external_transfer]
    can = external_balance_canary(model, X)
    @printf("%-24s resid %.2e  clearing %.2e\n", tag,
        maximum(abs, equilibrium_residuals(model, X)),
        maximum(abs, market_clearing_residuals(model, X)))
    @printf("%-24s p = %s\n", "", join(round.(sol.prices_raw; digits=16), ", "))
    @printf("%-24s q = %s\n", "", join(round.(sol.quantities; digits=16), ", "))
    @printf("%-24s w = %s\n", "", join(round.(sol.wages_raw; digits=16), ", "))
    @printf("%-24s F = %.16g  booked = %.16g  B_gov = %.16g  gap = %.2e  secgap = %.2e\n",
        "", sol.external_transfer, can.financing, can.programme_financing, can.diff,
        sectoral_labor_gap(model, sol.prices_raw, sol.quantities, sol.wages_raw))
    @printf("%-24s rgdp = %.16g  ngdp = %.16g\n", "", real_gdp(sol), nominal_gdp(sol))
end

data3 = three_sector_fixture()
sh0_3 = Shocks(ones(3), ones(3), zeros(3))
m0 = mobile_labor_model(data3, sh0_3, 0.5, 0.5, 0.9, 0.0; labor_bar=sum(data3.labor_share))
show_cell("3sec unshocked eta=0", m0, solve(m0))

supply = ones(3); supply[1] = 1.2
shA = Shocks(supply, ones(3); autonomous_demand=[0.1, 0.0, 0.0],
    investment_shock=[0.0, 0.05, 0.0])
mA = mobile_labor_model(data3, shA, 0.5, 0.5, 0.9, 0.0; labor_bar=sum(data3.labor_share))
show_cell("3sec additive eta=0", mA, solve(mA))

datat = tiny_fixture()
sh0_t = Shocks(ones(2), ones(2), zeros(2))
mt = mobile_labor_model(datat, sh0_t, 0.5, 0.5, 0.9, 0.0; labor_bar=sum(datat.labor_share))
show_cell("tiny unshocked eta=0", mt, solve(mt))

fx = v3_fixture()
shv = _v3_shocks()
for fin in (NoFinancing(), TaxFinanced(fx.g), ExternalDebt(fx.g))
    mv = mobile_labor_model(fx.data, shv, _V3_θ, _V3_ϵ, _V3_σ, 0.0; financing=fin)
    show_cell("v3 $(typeof(fin)) eta=0", mv, solve(mv))
end
