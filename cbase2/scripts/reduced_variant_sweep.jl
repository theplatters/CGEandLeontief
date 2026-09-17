# ═══════════════════════════════════════════════════════════════════════════════
# reduced_variant_sweep.jl — the decisive self-loop experiment + solver ladder
# (cbase2, 2026-09-17, per docs/DOCS_ASSESSMENT.md "post-v4 findings")
#
# Questions:
#  Q1 (dataset): does the "reduced" variant (71 + self-share>0.45 class
#     dropped, 87.4% GO coverage) converge where "full" and "70s" stall?
#     → falsifies / isolates the high-self-loop sectors as the cause.
#  Q2 (solver): does the new IPOPT rung (src/solvers.jl) crack the stalls
#     that Newton/LM leave at the 3.6e-4 floor?
#
# Run:  julia --threads=4 --project=. cbase2/scripts/reduced_variant_sweep.jl
# Rows appended to cbase2/results_intermediate/reduced_variant_sweep.csv
# ═══════════════════════════════════════════════════════════════════════════════

using DataFrames, CSV, Dates
const ROOT = "/workspace/git/BFRep/(3)BeyondHulten"
const CB = joinpath(ROOT, "cbase2")
for f in ["interface.jl", "solution.jl", "ces.jl", "mobile_labor.jl", "leontief.jl", "util.jl"]
	include(joinpath(CB, "src", "core", f))
end
include(joinpath(CB, "src", "closures.jl"))
include(joinpath(CB, "src", "financing.jl"))
include(joinpath(CB, "src", "calibration.jl"))
include(joinpath(CB, "src", "solvers.jl"))

rows = DataFrame(variant = String[], theta = Float64[], init = String[],
	stage = String[], method = String[], retcode = String[],
	resid = Float64[], canary = Float64[], secs = Float64[])

function record!(variant, th, iname, stage, meth, x, ret, t0)
	r = maximum(abs, equilibrium_residuals(mdl, x))  # mdl from enclosing scope
	(; nth) = residual_canary(mdl, x)
	push!(rows, (variant, th, iname, stage, meth, ret, r, nth,
				 round(time() - t0; digits = 2)))
	println(rpad(variant, 8), " θ=", rpad(th, 5), " init=", rpad(iname, 8),
		" ", rpad(stage * "/" * meth, 22), " ret=", rpad(ret, 28),
		" resid=", Printf.@sprintf("%.3e", r), " canary=", Printf.@sprintf("%.3e", nth),
		" (", round(time() - t0; digits = 1), "s)")
	return r, nth
end

raw = read_data(joinpath(CB, "data_raw", "I-O_DE2019_formatiert.csv"))
for (vname, drops) in [("full", DATASET_VARIANTS["full"]),
					   ("70s", DATASET_VARIANTS["70s"]),
					   ("reduced", DATASET_VARIANTS["reduced"])]
	data_v1 = drop_sectors(raw, drops)
	N = length(data_v1.factor_share)
	cover = dataset_coverage(data_v1, drops)
	println("── variant ", vname, ": N=", N, " coverage GO=", cover.gross_share_kept,
		"% VA=", cover.va_share_kept, "%")
	es = recalibrate_open(data_v1, CB; exo_scale = 1.0, drops = drops)
	shocks = Shocks(ones(N), ones(N), zeros(N))
	lin = fixed_point_init(es)
	for th in (2.0, 1.0, 0.5)
		global mdl = mobile_labor_model(es, shocks, th, 0.5, 0.9, 0.5)
		for (iname, x0) in [("default", [ones(N); es.λ; 1.0]), ("linear", lin)]
			t0 = time()
			# ── current setup: Newton + LM polish (standing `solve()` path) ──
			x, ret, _ = attempt_newton(mdl, x0; tol = 1e-10, maxiters = 20_000)
			r, nth = record!(vname, th, iname, "current", "newton", x, ret, t0)
			if r > 1e-6 || !canary_ok(nth)
				x2, ret2, _ = attempt_lm(mdl, x; tol = 1e-10, maxiters = 20_000)
				r2, nth2 = record!(vname, th, iname, "current", "lm-polish", x2, ret2, t0)
				if r2 > 1e-6 || !canary_ok(nth2)
					# ── new framework: IPOPT rung from the stalled point ──
					x3, info, r3 = ipopt_residual_solve(mdl, x2; tol = 1e-8,
						max_iter = 3000)
					record!(vname, th, iname, "ladder", info, x3, info, t0)
					if r3 > 1e-6 || !canary_ok(residual_canary(mdl, x3).nth_market)
						# ── full ladder (multistart + IPOPT) from scratch ──
						rr = solve_robust(mdl; init = nothing, tol = 1e-8,
							verbose = false, ladder = (:multistart, :ipopt))
						record!(vname, th, iname, "ladder", "solve_robust:" * rr.method,
							rr.x, "best", t0)
					end
				end
			end
		end
	end
end

out = joinpath(CB, "results_intermediate")
mkpath(out)
fname = joinpath(out, "reduced_variant_sweep.csv")
CSV.write(fname, rows)
println("\nrows: ", size(rows, 1), " → ", fname)

# summary: converged = resid ≤ 1e-6 AND canary ok
ok(row) = row.resid <= 1e-6 && abs(row.canary) < 1e-5
println("\n── convergence summary (resid ≤ 1e-6, canary ok) ──")
for vname in ("full", "70s", "reduced")
	sub = filter(r -> r.variant == vname && ok(r), rows)
	sol = filter(r -> r.variant == vname, rows)
	println(rpad(vname, 8), ": ", nrow(sub), "/", nrow(sol), " attempts converged; methods: ",
		join(unique(sub.method), ", "))
end
