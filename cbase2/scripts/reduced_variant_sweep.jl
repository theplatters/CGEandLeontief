# ═══════════════════════════════════════════════════════════════════════════════
# reduced_variant_sweep.jl — the decisive self-loop experiment (continuation
# recipe). (cbase2, 2026-09-17, per docs/DOCS_ASSESSMENT.md "post-v4 findings")
#
# Q (dataset): does the "reduced" variant (71 + self-share>0.45 class dropped,
#     87.4% GO coverage) behave differently from "full"/"70s" under the SAME
#     solver recipe the acceptance gate uses?
#
# Recipe (mirrors verify_v3.jl section 1 — the standing gate): bisection for
# the continuation start exo_scale* (saving rate crosses 0), exo_scale
# continuation K=6, θ ladder 2.0 → 0.5 within each step, every solve
# warm-started from the previous solution. A cold solve at each θ is recorded
# alongside for contrast (that is the configuration that stalls at ~1.7e-3).
#
# Run:  julia --threads=4 --project=. cbase2/scripts/reduced_variant_sweep.jl
# Rows → cbase2/results_intermediate/reduced_variant_sweep.csv
# ═══════════════════════════════════════════════════════════════════════════════

using DataFrames, CSV
using LinearAlgebra, NonlinearSolve
const CB = normpath(joinpath(@__DIR__, ".."))
for f in ["interface.jl", "solution.jl", "ces.jl", "mobile_labor.jl", "leontief.jl", "util.jl"]
	include(joinpath(CB, "src", "core", f))
end
include(joinpath(CB, "src", "closures.jl"))
include(joinpath(CB, "src", "financing.jl"))
include(joinpath(CB, "src", "calibration.jl"))
include(joinpath(CB, "src", "solvers.jl"))

rows = DataFrame(variant = String[], exo_scale = Float64[], theta = Float64[],
	mode = String[], init = String[], retcode = String[],
	resid = Float64[], canary = Float64[], wage = Float64[],
	maxp = Float64[], secs = Float64[])

function record!(variant, esc, th, mode, iname, x, ret, t0, mdl)
	r = maximum(abs, equilibrium_residuals(mdl, x))
	(; nth_market = nth, w) = residual_canary(mdl, x)
	push!(rows, (variant, esc, th, mode, iname, ret, r, nth, w,
				 maximum(abs, x[1:length(mdl.data.factor_share)] .- 1),
				 round(time() - t0; digits = 2)))
	println(rpad(variant, 8), " exo=", rpad(round(esc; digits = 4), 6),
		" θ=", rpad(th, 5), " ", rpad(mode * "/" * iname, 18),
		" ret=", rpad(ret, 12), " resid=", Printf.@sprintf("%.3e", r),
		" canary=", Printf.@sprintf("%.3e", nth),
		" w*=", Printf.@sprintf("%.4f", w))
	return r
end

THETAS = [2.0, 1.5, 1.2, 1.0, 0.8, 0.65, 0.5]   # the gate's ladder
K = 6
raw = read_data(joinpath(CB, "data_raw", "I-O_DE2019_formatiert.csv"))

for (vname, drops) in [("full", DATASET_VARIANTS["full"]),
					   ("70s", DATASET_VARIANTS["70s"]),
					   ("reduced", DATASET_VARIANTS["reduced"])]
	data_v1 = retained_dataset(raw, drops)
	N = length(data_v1.factor_share)
	shocks = Shocks(ones(N), ones(N), zeros(N))
	cover = dataset_coverage(raw, drops)
	println("══ variant ", vname, ": N=", N, "  coverage GO=", cover.gross_share_kept,
		"% VA=", cover.va_share_kept, "% ══")

	# continuation start: smallest exo_scale with s >= 0 (bisection, as the gate)
	s_of(esc) = (d = recalibrate_open(data_v1; exo_scale = esc); d.saving_rate)
	lo, hi = 0.0, 1.0
	@assert s_of(hi) > 0
	for _ in 1:40
		mid = (lo + hi) / 2
		s_of(mid) < 0 ? (lo = mid) : (hi = mid)
	end
	esc0 = hi
	println("continuation start exo_scale* = ", round(esc0; digits = 4),
		" (s = ", round(s_of(esc0); digits = 4), ")")

	init_warm = nothing
	for k in 0:K
		esc = esc0 + (1.0 - esc0) * k / K
		es = recalibrate_open(data_v1; exo_scale = esc)
		for th in THETAS
			mdl = mobile_labor_model(es, shocks, th, 0.5, 0.9, 0.5)
			t0 = time()
			if init_warm === nothing
				# first solve ever: linear fixed-point warm start (gate recipe)
				init_warm = fixed_point_init(es)
			end
			prob = NonlinearSolve.NonlinearProblem(problem, Float64.(init_warm), mdl)
			res = NonlinearSolve.solve(prob; reltol = 1e-6, abstol = 1e-6, maxiters = 20_000)
			x = Float64.(res.u)
			r = record!(vname, esc, th, "cont", "warm", x, string(res.retcode), t0, mdl)
			# bounded LM polish when above the gate (as solve() does)
			if r > 1e-6
				xl, retl, _ = attempt_lm(mdl, x; tol = 1e-10, maxiters = 20_000)
				r = record!(vname, esc, th, "cont", "lm", xl, retl, t0, mdl)
			end
			# IPOPT rung from the stalled point (the new framework's finisher)
			if r > 1e-6
				xi, infoi, ri = ipopt_residual_solve(mdl, x; tol = 1e-8,
					max_iter = 2000, formulation = :ls)
				r = record!(vname, esc, th, "cont", "ipopt", xi, infoi, t0, mdl)
			end
			init_warm = [x[1:N]; x[N+1:2N]; x[2N+1]]   # continue from best Newton point
		end
	end

	# cold-solve contrast: the configuration that stalls (θ ladder, direct init)
	es = recalibrate_open(data_v1; exo_scale = 1.0)
	for th in THETAS
		mdl = mobile_labor_model(es, shocks, th, 0.5, 0.9, 0.5)
		t0 = time()
		x, ret, _ = attempt_newton(mdl, default_init(mdl); tol = 1e-10, maxiters = 20_000)
		record!(vname, 1.0, th, "cold", "default", x, ret, t0, mdl)
	end
end

out = joinpath(CB, "results_intermediate")
mkpath(out)
fname = joinpath(out, "reduced_variant_sweep.csv")
CSV.write(fname, rows)
println("\nrows: ", size(rows, 1), " → ", fname)

# summary per variant: worst residual along the continuation and cold stalls
println("\n── summary ──")
for vname in ("full", "70s", "reduced")
	sub = filter(r -> r.variant == vname && r.mode == "cont", rows)
	cold = filter(r -> r.variant == vname && r.mode == "cold", rows)
	ok = filter(r -> r.resid <= 1e-6, sub)
	println(rpad(vname, 8), ": continuation steps resid≤1e-6: ", nrow(ok), "/", nrow(sub),
		"  worst cont resid = ", Printf.@sprintf("%.2e", maximum(sub.resid)),
		"  |  cold θ=0.5 resid = ",
		Printf.@sprintf("%.2e", only(cold[cold.theta .== 0.5, :resid])))
end
