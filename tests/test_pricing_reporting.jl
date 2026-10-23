using Test, CSV, DataFrames
isdefined(Main, :consumption_unit_prices) ||
    include(joinpath(@__DIR__, "..", "experiments", "make_pricing_prices.jl"))

@testset "pricing reporting: invariant to the solver's nominal scale" begin
    p = [0.94, 1.08, 1.03]
    beta = [0.2, 0.5, 0.3]
    for sigma in (0.0001, 0.9, 1.0)
        a = consumption_unit_prices(p, beta, sigma, 0.99)
        b = consumption_unit_prices(3.7 .* p, beta, sigma, 3.7 * 0.99)
        @test isapprox(a.gdp_deflator_cpi, b.gdp_deflator_cpi; atol = 1e-12)
        @test isapprox(a.max_abs_price_dev_cpi, b.max_abs_price_dev_cpi; atol = 1e-12)
        normalized = consumption_unit_prices(p ./ a.cpi_raw, beta, sigma, a.gdp_deflator_cpi)
        @test isapprox(normalized.cpi_raw, 1.0; atol = 1e-12)
    end
    @test_throws ErrorException consumption_unit_prices([0.0, 1.0], [0.5, 0.5], 0.9, 1.0)
end

@testset "pricing reporting: derived table and source gates" begin
    root = default_root()
    stored = CSV.read(joinpath(root, PRICE_SIDECAR), DataFrame)
    @test length(unique(stored.run_id)) == nrow(stored)
    # A full reconstruction additionally checks real GDP invariance on every
    # wage-normalized cell. Run it when machine-local solution artifacts exist.
    available = isfile(joinpath(root, "data", RUN_IO_TABLE)) &&
        all(isfile(joinpath(root, "runs", id, "solution.csv")) for id in stored.run_id)
    if available
        got = pricing_price_rows(; root)
        @test got.run_id == stored.run_id
        @test got.manifest_sha256 == stored.manifest_sha256
        @test got.solution_sha256 == stored.solution_sha256
        @test isapprox(got.gdp_deflator_cpi, stored.gdp_deflator_cpi; atol = 1e-12)
        @test isapprox(got.max_abs_price_dev_cpi, stored.max_abs_price_dev_cpi; atol = 1e-12)
    else
        @test_skip "machine-local pricing solutions unavailable"
    end
    # The table builder must use the normalized prices and reject stale
    # provenance, including when a clean clone has no solution.csv files.
    code = """
import importlib.util
spec = importlib.util.spec_from_file_location('pricing', 'revision/tables/make_pricing_door.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
rid = 'pricing_capacity_rows-v2-DELTA-F3-dm05-s0'
v = m.values(rid)
assert abs(v['P'] - (-0.025285014)) < 1e-7
assert abs(v['P'] - 100*(m.manifest(rid)['metrics']['gdp_deflator']-1)) > 0.2
# The merged Table 7 keeps its no-door reference rows and W_D sidecar.
assert m.COLS_MK[3][0] == 'WD'
for row, _ in m.ROWS:
    for f in m.FIN:
        baseline = f'matrix_5x3-v11-{row}-{f}' + ('-etas05' if row == 'BETA' else '')
        got = m.values(baseline)
        assert abs(got['P'] - 100*(m.manifest(baseline)['metrics']['gdp_deflator']-1)) < 1e-6
        assert abs(got['WD'] - m.WD[baseline]) < 1e-12
assert abs(m.values('pricing_profit-v3-GAMMA-F3-k05-spis')['WD'] - 3.40) < 0.005
table = m.build_mk_table()
assert table.count('no door') == 5
assert 'W_D' in table
m.reported_prices.cache_clear()
m.PRICES[rid]['manifest_sha256'] = 'stale'
try:
    m.reported_prices(rid)
except AssertionError as e:
    assert 'stale prices' in str(e)
else:
    raise AssertionError('stale price sidecar accepted')
"""
    @test success(Cmd(`python3 -c $code`; dir = root))
end
