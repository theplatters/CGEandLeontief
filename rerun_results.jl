# rerun_results.jl — Reproduce the definitive-guide Part I results under the NEW
# §4.1 data (domestic Ω, basic-price gross output, standard Domar λ).
# Runs WITHOUT the heavy GLMakie/XLSX/Ipopt deps (includes only core model files),
# so it precompiles fast. Robust: warm-starts from the no-shock reference and
# catches solver stalls (η≥1 and :fixed currently STALL under the §4.1 data).
#
# Run from the project root:
#   cd /workspace/BFrep/'(3)BeyondHulten' && julia --project=. rerun_results.jl
# (On a machine with a cached precompile this takes a few minutes; in a cold
#  container it may be killed by the 10-min foreground cap — run in background.)
#
# Outputs: rerun_results.log (incremental), output/variance_decomposition_sobol.csv

using CSV, DataFrames, LinearAlgebra, Statistics, NonlinearSolve, ProgressMeter, ThreadsX
include("src/interface.jl"); include("src/solution.jl"); include("src/ces.jl")
include("src/mobile_labor.jl"); include("src/util.jl"); include("src/variance_decomposition.jl")

cd(@__DIR__)
data = Data("I-O_DE2019_formatiert.csv"); N = 71; SEC = 35
logio = open("rerun_results.log", "w")
log(msg) = (println(msg); println(logio, msg); flush(logio))
log("NEW §4.1 DATA: sum(λ)=$(round(sum(data.λ); digits=4))  sum(labor_share)=$(round(sum(data.labor_share); digits=4))")

# Reference (no shock) for Δ% normalisation + warm-start
ref = mobile_labor_model(data, Shocks(ones(N), ones(N), zeros(N)), 0.5, 0.5, 0.9, 0.5)
ref_sol = solve(ref); ref_gdp = real_gdp(ref_sol)
init_ref = [ref_sol.prices_raw; ref_sol.quantities; ref_sol.wages_raw[1]]
log("ref_gdp = $ref_gdp")

# ── Part I §3: binary wage regime (:fixed) vs continuous η (autonomous demand) ──
log("\n=== §3 wage-regime table (autonomous demand to construction) ===")
log("mult\tη=0(Δ%)\tη=1(Δ%)\tfixed(Δ%)\tη_bridge\tfix_bridge\tfixed_emp\tfixed_res")
for mult in [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
    sh = autonomous_shock(data; autonomous_mult = mult)
    g0 = g1 = gf = ef = mrf = NaN
    try; g0 = 100*(real_gdp(solve(mobile_labor_model(data, sh, 0.5,0.5,0.9,0.0; closure=:mobile); init=init_ref))/ref_gdp - 1); catch; end
    try; g1 = 100*(real_gdp(solve(mobile_labor_model(data, sh, 0.5,0.5,0.9,1.0; closure=:mobile); init=init_ref))/ref_gdp - 1); catch; end
    try
        sf = solve(mobile_labor_model(data, sh, 0.5,0.5,0.9,0.5; closure=:fixed); init=init_ref)
        gf = 100*(real_gdp(sf)/ref_gdp - 1); ef = sum(sectoral_labor_demand(sf.prices, sf.quantities, 1.0, sf.model))
        mrf = maximum(abs.(equilibrium_residuals(sf)))
    catch; end
    log("$(mult)\t$(round(g0;digits=6))\t$(round(g1;digits=6))\t$(round(gf;digits=6))\t$(round(g1-g0;digits=6))\t$(round(gf-g0;digits=6))\t$(round(ef;digits=4))\t$(mrf)")
end

# ── Part I §2: Sobol on (η,ε,θ,σ) — supply +30% to sector 35 ──
log("\n=== §2 Sobol decomposition (supply +30% to sector 35) ===")
ss = ones(N); ss[SEC] *= 1.30
sh = Shocks(ss, ones(N), zeros(N))
mkpath("output")
try
    vd = variance_decomposition(data, sh;
        η_values=[0.0,0.5,1.0,2.0], ϵ_values=[0.5,2.0], θ_values=[0.5,0.99], σ_values=[0.5,0.99],
        output=:real_gdp, verbose=true)
    summary_table(vd; save_csv=joinpath("output","variance_decomposition_sobol.csv"))
    log("η S_f=$(round(vd.S_f["η"];digits=4))  ST_f=$(round(vd.ST_f["η"];digits=4))")
    log("Sobol S_f: $([(f=>round(vd.S_f[f];digits=4)) for f in vd.factors])")
    log("Sobol ST_f: $([(f=>round(vd.ST_f[f];digits=4)) for f in vd.factors])")
    log("n_failed=$(vd.n_failed) / $(nrow(vd.grid)) grid points")
catch e
    log("Sobol FAILED: $e")
end
log("DONE"); close(logio)
