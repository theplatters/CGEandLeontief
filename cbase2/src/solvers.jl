# ═══════════════════════════════════════════════════════════════════════════════
# cbase2 solver framework (Stage 1, 2026-09-17) — src/solvers.jl
# ═══════════════════════════════════════════════════════════════════════════════
#
# Robust equilibrium-solving harness, additive to the core (the standing
# `solve()` in `src/core/mobile_labor.jl` is untouched). Motivated by the
# residual-floor finding (FD-Newton stalls at ~3.6e-4 under some inits while
# the identical system converges to 2e-10 under another — init/path-sensitive,
# not a model-structure problem) and by the Baqaee–Farhi practice of solving
# equilibria as general NLPs (KNITRO/fmincon in MATLAB; SQP/interior-point
# globalization: line searches, trust regions, feasible-point maintenance).
#
# The ladder (each rung only runs if the previous one missed the gate):
#   A  default NonlinearSolve (TRF/LevenbergMarquardt as chosen by the library)
#      + bounded LM polish — exactly the standing `solve()` behaviour;
#   B  LM-primary with ForwardDiff, then LM with finite-difference AD;
#   C  multi-start battery (linear fixed-point init, default init, jittered
#      variants) through Newton and LM;
#   D  IPOPT: the square system as a constrained NLP —
#          minimize 0.5·‖r(x)‖²   s.t.  r(x) = 0,   x ≥ 0,
#      with dense central-difference Jacobians (no ForwardDiff, so the log/
#      floor clamps are never differentiated through) and the limited-memory
#      Hessian approximation. Positivity bounds keep the iterates strictly
#      inside the region where the residual's internal floors are inactive —
#      this is the feasible-point maintenance the FD-Newton lacks. IPOPT is
#      the free KNITRO analogue and the standard in the CGE world.
#
# Nothing in the core changes; the notebooks may adopt `solve_robust` per
# experiment after the acceptance gate (`verify_v3.jl`) passes on it.
#
# Public API:
#   fixed_point_init(data)                  warm start from the linear (θ=1)
#                                           fixed point [ones(N); ylin; 1]
#   solve_robust(model; init, tol, verbose) → RobustResult(x, resid, method,
#                                           attempts, report DataFrame)
#   solution_from_x(model, x)               Solution built from a raw vector
#                                           (post-solve consumption/Tornqvist)
# ═══════════════════════════════════════════════════════════════════════════════

using LinearAlgebra, Printf, Random, Statistics
using NonlinearSolve

# IPOPT rung D needs JuMP + Ipopt (both now in Project.toml); degrade
# gracefully when absent so scripts that only need rungs A–C still load.
const _IPOPT_AVAILABLE = Ref(false)
try
	@eval using JuMP, Ipopt
	_IPOPT_AVAILABLE[] = true
catch
	@warn "JuMP/Ipopt not available — solve_robust ladder rung D disabled"
end

# ──────────────────────────────────────────────────────────────────────────────
# Warm starts
# ──────────────────────────────────────────────────────────────────────────────

"""
	fixed_point_init(data::Data) -> Vector{Float64}

Linear (θ = 1) fixed-point warm start: solve the open-economy Leontief-style
linear demand system for gross output and append p = 1, w = 1. This is the
"linear" init of the solver sweep — cheap, and under a nearby init the
identical system has been observed to converge to 2e-10.
"""
function fixed_point_init(data::Data)
	N = length(data.factor_share)
	M0 = data.Ω_raw' * Diagonal(1.0 .- data.factor_share)
	Gk = M0 + Diagonal(1.0 .- data.import_margin) * data.consumption_share *
		 data.factor_share' * (1 - data.saving_rate) * (1 - sum(data.gov_demand))
	ylin = (I - Gk) \ ((1.0 .- data.import_margin) .* (data.gov_demand .+ data.exo_demand) .+
					   data.exports_demand)
	return [ones(N); ylin; 1.0]
end

"""
	default_init(model) -> Vector{Float64}

Standing default init: p = 1, y = λ, w = 1.
"""
default_init(model::Model{MobileLaborCES}) =
	[ones(length(model.data.factor_share)); model.data.λ; 1.0]

"""
	init_battery(model; seed = 20260917) -> Vector{Tuple{String,Vector{Float64}}}

Deterministic init battery for the multi-start rung: default, linear fixed
point, two sine-jittered and two seeded-random perturbations of each.
"""
function init_battery(model::Model{MobileLaborCES}; seed::Int = 20260917)
	N = length(model.data.factor_share)
	d = default_init(model)
	l = fixed_point_init(model.data)
	rng = MersenneTwister(seed)
	b = Tuple{String,Vector{Float64}}[
		("default", d), ("linear", l),
		("lin+sin", [l[1:N]; l[N+1:2N] .+ 0.02 .* sin.(1:N); l[2N+1]]),
		("def+sin", [d[1:N]; d[N+1:2N] .+ 0.02 .* sin.(1:N); d[2N+1]]),
	]
	for k in 1:2
		eps1 = 0.02 .* (rand(rng, N) .- 0.5)
		eps2 = 0.02 .* (rand(rng, N) .- 0.5)
		push!(b, ("lin+r$(k)", [l[1:N]; l[N+1:2N] .* (1.0 .+ eps1); l[2N+1]]))
		push!(b, ("def+r$(k)", [d[1:N]; d[N+1:2N] .* (1.0 .+ eps2); d[2N+1]]))
	end
	return b
end

# ──────────────────────────────────────────────────────────────────────────────
# Residual utilities
# ──────────────────────────────────────────────────────────────────────────────

"""
	residual_canary(model, x) -> (rmax, nth_market, w)

Max residual of the enforced system plus the ACCOUNTING CANARY: the
saving-identity residual S − (I + X − M), i.e. s·E − p′(I+X) + M(x) with M
the import content of gross household, government and investment demand.
Since the 2026-09-17 formulation repair ALL N markets are enforced and the
external balance closes residually, this residual vanishes with the system
residual at any true equilibrium — asserting it in the gate catches exactly
the formulation regression (dropped-market over-determination) that
previously certified non-equilibria with machine-zero residuals.
"""
function residual_canary(model::Model{MobileLaborCES}, x::AbstractVector)
	N = length(model.data.factor_share)
	r = equilibrium_residuals(model, x)
	(; data, options, shocks) = model
	p = _positive_floor(x[1:N]); y = _positive_floor(x[N+1:2N])
	w = max(x[2N+1], 1e-10)
	(; σ,) = options.elasticities
	L_i = sectoral_labor_demand(p, y, w, model)
	fin = model.financing
	ds_eff = preference_weights(fin, shocks.demand_shock)
	L_sum = sum(L_i)
	E = household_expenditure(fin, model, w * L_sum, p, L_sum)
	agg = sum(data.consumption_share .* ds_eff .* p .^ (1 - σ))
	cg_gross = (1 .- data.saving_rate) .*
			   (data.consumption_share .* ds_eff) .* E .* p .^ (-σ) ./ agg
	M = dot(p, data.import_margin .* (cg_gross .+ data.gov_demand .+ data.exo_demand))
	saving_identity = data.saving_rate * E - dot(p, data.exo_demand .+ data.exports_demand) + M
	return (rmax = maximum(abs, r), nth_market = saving_identity, w = w)
end

"Canary gate: the saving identity S = I + X − M must hold at the equilibrium."
canary_ok(nth) = abs(nth) < 1e-5

# ──────────────────────────────────────────────────────────────────────────────
# Report type
# ──────────────────────────────────────────────────────────────────────────────

"""
	RobustResult

Outcome of `solve_robust`: the best `x`, its max residual and canary value,
the winning method label, and the full attempt log (each entry: method, init,
retcode, residual, canary, seconds, message).
"""
struct RobustResult
	x::Vector{Float64}
	resid::Float64
	nth_market::Float64
	method::String
	attempts::Vector{NamedTuple}
end

Base.show(io::IO, r::RobustResult) = begin
	println(io, "RobustResult: method=", r.method, "  max|resid|=", r.resid,
			"  canary(N-th mkt)=", r.nth_market, "  attempts=", length(r.attempts))
end

_attempt!(rr, kw...) = push!(rr, (; kw...))

# ──────────────────────────────────────────────────────────────────────────────
# Rungs A–C: NonlinearSolve attempts
# ──────────────────────────────────────────────────────────────────────────────

function _nsolve(model::Model{MobileLaborCES}, x0, alg; tol, maxiters)
	prob = NonlinearSolve.NonlinearProblem(problem, Float64.(x0), model)
	res = alg === nothing ? NonlinearSolve.solve(prob; reltol = tol, abstol = tol, maxiters) :
		NonlinearSolve.solve(prob, alg; reltol = tol, abstol = tol, maxiters)
	x = Float64.(res.u)
	rmax = maximum(abs, equilibrium_residuals(model, x))
	return x, string(res.retcode), rmax
end

function attempt_newton(model, x0; tol = 1e-10, maxiters = 20_000)
	x, ret, rmax = _nsolve(model, x0, nothing; tol, maxiters)
	return x, ret, rmax
end

function attempt_lm(model, x0; tol = 1e-10, maxiters = 20_000, ad = :default)
	# 2026-09-17 (review §4 / LaForge): `ad` was an INVALID keyword for
	# NonlinearSolve v4 — the correct kwarg is `autodiff`. The old try/catch
	# silently swallowed the MethodError, so the :fd arm ran ForwardDiff all
	# along and every past FD-vs-AD comparison was fake. No try/catch: a bad
	# algorithm construction must fail loudly.
	alg = ad === :fd ? LevenbergMarquardt(; autodiff = AutoFiniteDiff()) :
					  LevenbergMarquardt()
	x, ret, rmax = _nsolve(model, x0, alg; tol, maxiters)
	return x, ret, rmax
end

"""
	attempt_projected_newton(model, x0; tol, maxiters, lo, hi) -> (x, ret, rmax)

Damped Newton with box projection and central-difference Jacobians — the
globalization the plain FD-Newton lacks. Motivation (measured, 2026-09-17):
the stall attractor is a path excursion into the y→0 cliff — at the gate's
failing first rung the stalled point has y₅ = 2.1e-6 (240× below λ₅), where
the residual's log-sensitivity makes the FD Jacobian catastrophically
ill-conditioned (σmax 2e15, cond 3.5e20) and every Newton/IPOPT step
diverges. Fencing the iterates into the economically meaningful region
(relative lower bounds on y, generous bounds on p and w) keeps the Jacobian
sane; the backtracking line search on ‖r‖₂ does the rest.

Default fence: p ≥ 0.05, y_i ≥ 1e-3·λ_i, w ≥ 0.05, no upper bounds except
p ≤ 50, y ≤ 200·λ_i (far outside any admissible equilibrium; blocks the
price-explosion branch).
"""
function attempt_projected_newton(model::Model{MobileLaborCES}, x0;
		tol::Float64 = 1e-10, maxiters::Int = 200,
		lo::Union{Nothing,Vector{Float64}} = nothing,
		hi::Union{Nothing,Vector{Float64}} = nothing)
	nv = length(x0)
	N = length(model.data.factor_share)
	lo = lo === nothing ? vcat(fill(0.05, N), 1e-3 .* max.(model.data.λ, 1e-6), [0.05]) : lo
	hi = hi === nothing ? vcat(fill(50.0, N), 200.0 .* max.(model.data.λ, 1e-6), [50.0]) : hi
	project!(x) = (x .= clamp.(x, lo, hi); x)
	r! = (out, xx) -> (problem(out, xx, model); out)
	nv == length(lo) == length(hi) || throw(DimensionMismatch("fence size mismatch"))
	r = Vector{Float64}(undef, nv); buf = Vector{Float64}(undef, nv)
	buf2 = Vector{Float64}(undef, nv); trial = Vector{Float64}(undef, nv)
	J = Matrix{Float64}(undef, nv, nv)
	x = project!(Float64.(x0))
	rmax = maximum(abs, r!(r, x))
	ret = "MaxIters"
	for it in 1:maxiters
		rmax = maximum(abs, r)
		rmax <= tol && (ret = string("converged(it=", it, ")"); break)
		_fd_jacobian!(J, r!, x, buf, buf2)
		# damped Newton direction: least-squares fallback when J is singular
		dx = try
			-(J \ r)
		catch
			-(J' * ((J * J' + 1e-8 * I) \ r))
		end
		all(isfinite, dx) || (dx = -(J' * ((J * J' + 1e-6 * I) \ r)))
		gn = sqrt(dot(r, r))
		# backtracking line search on ‖r‖₂ with projection after every trial
		t = 1.0
		improved = false
		for _ in 1:40
			trial .= x .+ t .* dx
			project!(trial)
			r!(r, trial)
			if sqrt(dot(r, r)) < (1 - 1e-4 * t) * gn
				improved = true
				break
			end
			t /= 2
		end
		improved || (ret = string("stalled(it=", it, ")"); break)
		x .= trial
	end
	rmax = maximum(abs, r!(r, x))    # r may hold a rejected trial's residual
	return x, ret, rmax
end

"""
	attempt_logspace_newton(model, x0; tol, maxiters) -> (x, ret, rmax)

Damped Newton in LOG coordinates: variables z = [log p; log y; log w], so
positivity holds by construction and the eps-floor cliff (the measured cause
of every stall: a y-component crossing zero under an O(0.1) step while λ
itself is 5e-4) is unreachable — in log space the residual is smooth along
the whole path (the floor at y = eps sits ~33 log-units below the data).
This is the standard CGE practice of solving in relative changes (MPSGE).

FD Jacobian in z-space; backtracking line search on ‖r‖₂. Returns the
x-space vector. This is the preferred opener rung; the fence of
`attempt_projected_newton` turned out to be the wrong cure (a projected
step at the fence boundary cannot descend, and any fence low enough to
permit convergence sits inside the cliff).
"""
function attempt_logspace_newton(model::Model{MobileLaborCES}, x0;
		tol::Float64 = 1e-10, maxiters::Int = 500)
	N = length(model.data.factor_share)
	nv = 2N + 1
	floorlog = log(1e-12)
	z0 = vcat(log.(max.(x0[1:2N], 1e-12)), [log(max(x0[2N+1], 1e-12))])
	xf(z) = vcat(exp.(z[1:2N]), [exp(z[2N+1])])
	rlog! = (out, z) -> (problem(out, xf(z), model); out)
	r = Vector{Float64}(undef, nv); buf = Vector{Float64}(undef, nv)
	buf2 = Vector{Float64}(undef, nv); trial = Vector{Float64}(undef, nv)
	J = Matrix{Float64}(undef, nv, nv)
	z = copy(z0)
	rmax = maximum(abs, rlog!(r, z))
	ret = "MaxIters"
	for it in 1:maxiters
		rmax = maximum(abs, r)
		rmax <= tol && (ret = string("converged(it=", it, ")"); break)
		_fd_jacobian!(J, rlog!, z, buf, buf2)
		dx = try
			-(J \ r)
		catch
			-(J' * ((J * J' + 1e-8 * I) \ r))
		end
		all(isfinite, dx) || (dx = -(J' * ((J * J' + 1e-6 * I) \ r)))
		gn = sqrt(dot(r, r))
		t = 1.0
		improved = false
		for _ in 1:60
			trial .= z .+ t .* dx
			all(isfinite, trial) || (t /= 2; continue)
			rlog!(r, trial)
			if sqrt(dot(r, r)) < (1 - 1e-4 * t) * gn
				improved = true
				break
			end
			t /= 2
		end
		improved || (ret = string("stalled(it=", it, ")"); break)
		z .= trial
		z[1:2N] = max.(z[1:2N], floorlog)   # keep exp() finite
	end
	x = xf(z)
	rmax = maximum(abs, rlog!(r, z))
	return x, ret, rmax
end

# ──────────────────────────────────────────────────────────────────────────────
# Rung D: IPOPT constrained-NLP formulation
# ──────────────────────────────────────────────────────────────────────────────

"""
Central-difference Jacobian of `r!` at `x`, written into the (n×m) matrix `J`
(rows = equations, columns = variables). Step h_i = 6e-6·max(1, |x_i|) — the
central-difference optimum for smooth functions; deliberate choice, NOT
ForwardDiff, so the internal log/floor clamps of the residual are evaluated,
never differentiated.
"""
function _fd_jacobian!(J::Matrix{Float64}, r!::Function, x::Vector{Float64},
		buf::Vector{Float64}, buf2::Vector{Float64})
	n = length(x); m = length(buf)
	for j in 1:n
		h = 6e-6 * max(1.0, abs(x[j]))
		xp = copy(x); xp[j] += h
		xm = copy(x); xm[j] -= h
		r!(buf, xp); r!(buf2, xm)
		for i in 1:m
			J[i, j] = (buf[i] - buf2[i]) / (2h)
		end
	end
	return J
end

"""
	ipopt_residual_solve(model, x0; tol, max_iter, verbose, pbnd, ubnd,
	                     formulation) -> (x, info::String, rmax)

IPOPT formulation (rung D): the equilibrium system as a bound-constrained
NLP with dense central-difference derivatives (one FD Jacobian cached per
trial point, shared by the objective gradient and the constraint rows) and
the limited-memory Hessian approximation (the KNITRO/fmincon analogue; free,
mature, standard in the CGE world). Implemented through JuMP nonlinear
operators with user-supplied gradients (IPOPT.jl v1.x removed the raw
CreateProblem API).

Two formulations:
  :square  minimize 0.5‖r(x)‖²  subject to  r(x) = 0,  pbnd ≤ x ≤ ubnd
           (classic square-system NLP; IPOPT reports the constraint
           violation directly).
  :ls      minimize 0.5‖r(x)‖²  subject to  pbnd ≤ x ≤ ubnd
           (bound-constrained nonlinear least squares; the line search
           cannot accept an objective increase, so the price-explosion
           branch is rejected outright rather than explored).

Bounds are load-bearing here: from the default init IPOPT's first full step
lands ON the documented price-explosion branch (objective 3.4e+20 — measured),
and plain filter logic then drifts near-feasible with a growing objective.
The upper bound ubnd (default 100, equilibrium prices/output live near 1)
makes that branch infeasible; `bound_push`/`bound_frac` are set to 1e-9 so
the start point is not shifted away from x0 (IPOPT's default 0.01 push
distorts small lambda components massively).

Returns the best-residual iterate among {IPOPT's final iterate, the warm
start} plus a status string. Requires JuMP + Ipopt in the project.
"""
function ipopt_residual_solve(model::Model{MobileLaborCES}, x0::Vector{Float64};
		tol::Float64 = 1e-8, max_iter::Int = 3000, verbose::Bool = false,
		pbnd::Float64 = 1e-8, ubnd::Float64 = 100.0,
		formulation::Symbol = :square)
	@assert _IPOPT_AVAILABLE[] "ipopt_residual_solve requires JuMP and Ipopt"
	formulation in (:square, :ls) ||
		throw(ArgumentError("formulation must be :square or :ls, got $formulation"))
	nv = length(x0)                    # 2N+1
	N = length(model.data.factor_share)
	# ── Scaling (percent-deviation variables for quantities) ──
	# The system mixes quantity scales of O(1) and ~5e-4 (small sectors); a
	# perturbation that drives a tiny y_i across zero hits the eps floor, log y
	# jumps by ~-28 and the wedge penalty (proportional to ratio²) explodes the
	# residual by 10 orders of magnitude (measured: r = 7.8e11 at a 0.01
	# absolute y shift). Scaling y by λ (the baseline gross output) puts every
	# decision variable at O(1) — the standard CGE practice (MPSGE solves in
	# relative deviations). p and w are already O(1) and stay in levels.
	scale = vcat(ones(N), max.(model.data.λ, 1e-6), [1.0])
	buf = Vector{Float64}(undef, nv)
	buf2 = Vector{Float64}(undef, nv)
	# per-point cache: r(x) and the FD Jacobian, computed once and shared by
	# the objective gradient and every constraint row. Evaluation happens in
	# SCALED coordinates v (x = scale .* v): MOI works in v, the residual in x.
	key_v = fill(NaN, nv)
	r_cache = fill(NaN, nv)
	J_cache = fill(NaN, nv, nv)     # ∂r/∂v (v-space Jacobian)
	lock_r = ReentrantLock()

	rv! = (out, v) -> (problem(out, scale .* v, model); out)
	function _ensure(v::Vector{Float64})
		lock(lock_r) do
			if key_v != v
				rv!(r_cache, v)
				_fd_jacobian!(J_cache, rv!, v, buf, buf2)
				key_v[:] = v
			end
		end
	end
	f_impl(x...) = (v = collect(Float64, x); _ensure(v); 0.5 * dot(r_cache, r_cache))
	# NOTE: MOI hands the gradient out as an `_UnsafeVectorView`, which BLAS
	# `mul!` cannot point into — so assign via broadcast (tiny per-call alloc,
	# negligible against the ~2n residual evaluations of the FD Jacobian).
	f_grad(gout, x...) = (v = collect(Float64, x); _ensure(v); gout .= J_cache' * r_cache; nothing)

	v0 = x0 ./ scale                  # warm start in scaled coordinates

	m = JuMP.Model(Ipopt.Optimizer)
	JuMP.set_silent(m)
	JuMP.set_attribute(m, "hessian_approximation", "limited-memory")
	JuMP.set_attribute(m, "mu_strategy", "adaptive")
	JuMP.set_attribute(m, "tol", tol)
	JuMP.set_attribute(m, "acceptable_tol", max(1e-6, 10 * tol))
	JuMP.set_attribute(m, "constr_viol_tol", tol)
	JuMP.set_attribute(m, "acceptable_constr_viol_tol", max(1e-5, 100 * tol))
	JuMP.set_attribute(m, "max_iter", max_iter)
	JuMP.set_attribute(m, "bound_relax_factor", 0.0)
	JuMP.set_attribute(m, "bound_push", 1e-9)
	JuMP.set_attribute(m, "bound_frac", 1e-9)
	JuMP.set_attribute(m, "sb", "yes")
	if verbose
		JuMP.unset_silent(m)
		JuMP.set_attribute(m, "print_level", 5)
	end

	JuMP.@variable(m, vv[i = 1:nv], start = v0[i])
	JuMP.set_lower_bound.(vv, pbnd ./ scale)   # bounds in x-units → v-units
	JuMP.set_upper_bound.(vv, ubnd ./ scale)   # excludes the explosion branch
	op_f = JuMP.add_nonlinear_operator(m, nv, f_impl, f_grad)
	JuMP.@objective(m, Min, op_f(vv...))
	if formulation === :square
		for i in 1:nv
			ri(x...) = (v = collect(Float64, x); _ensure(v); r_cache[i])
			rgi(gout, x...) = (v = collect(Float64, x); _ensure(v); gout .= J_cache[i, :]; nothing)
			# UNIQUE name is load-bearing: closures from the same source line
			# share one mangled name, and Symbol(f) collisions make the MOI
			# registry overwrite earlier registrations — all 141 constraints
			# would silently evaluate as the LAST row (the numeraire residual
			# ~1e-15), so IPOPT sees a feasible system and never solves it.
			op = JuMP.add_nonlinear_operator(m, nv, ri, rgi;
				name = Symbol("res_", i))
			JuMP.@constraint(m, op(vv...) == 0)
		end
	end
	JuMP.optimize!(m)
	st = JuMP.termination_status(m)
	v_final = try Float64.(JuMP.value.(vv)) catch; copy(v0) end
	x = scale .* v_final                  # back to x-space
	info = "ipopt($st)"
	# best-of-iterate: IPOPT's final iterate can drift off the best residual
	# seen along the path; compare it against the warm start.
	cands = [(x, maximum(abs, equilibrium_residuals(model, x))),
			 (Float64.(x0), maximum(abs, equilibrium_residuals(model, x0)))]
	best = argmin([c[2] for c in cands])
	return cands[best][1], info, cands[best][2]
end

# ──────────────────────────────────────────────────────────────────────────────
# The ladder
# ──────────────────────────────────────────────────────────────────────────────

"""
	solve_robust(model; init = nothing, tol = 1e-8, verbose = true,
	             ladder = (:newton, :lm, :multistart, :ipopt), seed = 20260917)
	             -> RobustResult

Run the robustness ladder on `model` and return the best raw vector found.
The gate for stopping early is `max|resid| ≤ tol` AND the Walras canary
(N-th clearing residual) below 1e-5. Every attempt is logged in the result's
`attempts` vector. Never throws on failure — return the best attempt so the
caller (sweep scripts, notebooks) can decide.
"""
function solve_robust(model::Model{MobileLaborCES};
		init = nothing, tol::Float64 = 1e-8, verbose::Bool = true,
		ladder::Tuple = (:pnewton, :newton, :lm, :multistart, :ipopt),
		seed::Int = 20260917, maxiters::Int = 20_000, ipopt_max_iter::Int = 3000)
	attempts = Vector{NamedTuple}()
	best_x = Float64[]; best_r = Inf; best_m = "none"
	gate(r, nth) = r <= tol && canary_ok(nth)

	# Accept a candidate if it beats the incumbent on (resid, canary).
	function offer!(x, meth, ret, t0, initname)
		r = maximum(abs, equilibrium_residuals(model, x))
		(; nth) = residual_canary(model, x)
		push!(attempts, (method = meth, init = initname, retcode = ret,
			resid = r, canary = nth, secs = round(time() - t0; digits = 2), msg = ""))
		verbose && println("  [", rpad(meth, 16), "init=", rpad(initname, 8),
			"] ret=", rpad(ret, 10), " resid=", Printf.@sprintf("%.3e", r),
			" canary=", Printf.@sprintf("%.3e", nth))
		better = isempty(best_x) || r < best_r ||
				 (r ≈ best_r && canary_ok(nth) && !canary_ok(residual_canary(model, best_x).nth_market))
		if better
			best_x = Float64.(x); best_r = r; best_m = meth
		end
		return r, nth
	end

	verbose && println("solve_robust: tol=", tol, " ladder=", ladder)

	t0 = time()
	# Rung A2 (first): projected damped Newton — fences the y→0 cliff that
	# every unprojected method wanders into (measured stall anatomy 2026-09-17).
	if :pnewton in ladder
		for (nm, x0) in inits
			x, ret, _ = attempt_projected_newton(model, x0; tol, maxiters = 500)
			r, nth = offer!(x, "pnewton", ret, t0, nm)
			gate(r, nth) && return RobustResult(best_x, best_r, nth, "pnewton($nm)", attempts)
		end
	end

	# Rung A/B/C share the init sequence; rung D starts from the best so far.
	inits = init === nothing ? init_battery(model; seed) :
			[("given", Float64.(init))]

	if :newton in ladder
		for (nm, x0) in inits
			x, ret, _ = attempt_newton(model, x0; tol, maxiters)
			r, nth = offer!(x, "newton", ret, t0, nm)
			gate(r, nth) && return RobustResult(best_x, best_r, nth, "newton($nm)", attempts)
		end
	end

	if :lm in ladder
		for adsym in (:default, :fd)
			for (nm, x0) in inits
				x, ret, _ = attempt_lm(model, x0; tol, maxiters, ad = adsym)
				r, nth = offer!(x, "lm/$adsym", ret, t0, nm)
				gate(r, nth) && return RobustResult(best_x, best_r, nth, "lm/$adsym($nm)", attempts)
			end
		end
	end

	# multistart: jittered re-launches of Newton/LM from the best point so far
	if :multistart in ladder
		rng = MersenneTwister(seed + 1)
		base = isempty(best_x) ? inits[1][2] : best_x
		for k in 1:4
			N = length(model.data.factor_share)
			jit = [base[1:N]; base[N+1:2N] .* (1.0 .+ 0.05 .* (rand(rng, N) .- 0.5));
				   base[2N+1]]
			for (meth, fn) in (("jitter+newton", attempt_newton), ("jitter+lm", attempt_lm))
				x, ret, _ = fn(model, jit; tol, maxiters)
				r, nth = offer!(x, meth, ret, t0, "j$(k)")
				gate(r, nth) && return RobustResult(best_x, best_r, nth, "$meth($(k))", attempts)
			end
		end
	end

	# Rung D: IPOPT from the best-so-far point. :ls first (the line search
	# rejects the price-explosion branch outright), :square as the alternative.
	if :ipopt in ladder
		start = isempty(best_x) ? inits[1][2] : best_x
		if isempty(start) || any(!isfinite, start)
			start = inits[1][2]
		end
		for f in (:ls, :square)
			x, info, r = ipopt_residual_solve(model, start; tol,
				max_iter = ipopt_max_iter, verbose = false, formulation = f)
			(; nth) = residual_canary(model, x)
			push!(attempts, (method = "ipopt/$f", init = "best-so-far", retcode = info,
				resid = r, canary = nth, secs = round(time() - t0; digits = 2), msg = ""))
			verbose && println("  [", rpad("ipopt/$f", 16), "init=best-sof ] ret=", rpad(info, 10),
				" resid=", Printf.@sprintf("%.3e", r), " canary=", Printf.@sprintf("%.3e", nth))
			if r < best_r || (r ≈ best_r && canary_ok(nth))
				best_x = x; best_r = r; best_m = "ipopt/$f"
			end
			best_r <= tol && break
		end
		# final polish: one AD-Newton/LM pass from the IPOPT point (smooth
		# region, so AD is safe here and drives to machine precision)
		if best_r > tol
			for (meth, fn) in (("ipopt+newton", attempt_newton), ("ipopt+lm", attempt_lm))
				x2, ret2, _ = fn(model, best_x; tol = 1e-12, maxiters)
				r2 = maximum(abs, equilibrium_residuals(model, x2))
				(; nth2) = residual_canary(model, x2)
				push!(attempts, (method = meth, init = "post-ipopt", retcode = ret2,
					resid = r2, canary = nth2, secs = round(time() - t0; digits = 2), msg = ""))
				verbose && println("  [", rpad(meth, 16), "init=post-ipopt] ret=", rpad(ret2, 10),
					" resid=", Printf.@sprintf("%.3e", r2), " canary=", Printf.@sprintf("%.3e", nth2))
				if r2 < best_r
					best_x = x2; best_r = r2; best_m = meth
				end
				gate(r2, nth2) && break
			end
		end
	end

	nth_final = isempty(best_x) ? NaN : residual_canary(model, best_x).nth_market
	return RobustResult(best_x, best_r, nth_final, best_m, attempts)
end

# ──────────────────────────────────────────────────────────────────────────────
# Solution construction from a raw vector
# ──────────────────────────────────────────────────────────────────────────────

"""
	solution_from_x(model, x) -> Solution

Build the full `Solution` (prices, quantities, wage, financed consumption,
Tornqvist real GDP) from a raw equilibrium vector `x`, mirroring the
post-solve construction in `solve()` so the harness output is
report-equivalent to a converged `solve()` call.
"""
function solution_from_x(model::Model{MobileLaborCES}, x::AbstractVector)
	(; data, options, shocks) = model
	N = length(data.factor_share)
	p = x[1:N]
	q = max.(x[N+1:2N], 0.0)
	w = x[2N+1]
	L_i = sectoral_labor_demand(p, q, w, model)
	wages = fill(w, N)
	numeraire = (data.consumption_share' * p .^ (1 - options.elasticities.σ)) ^
				(1 / (1 - options.elasticities.σ))
	fin = model.financing
	ds_eff = preference_weights(fin, shocks.demand_shock)
	L_sum = sum(L_i)
	total_income = w * L_sum
	E = household_expenditure(fin, model, total_income, p, L_sum)
	σ = options.elasticities.σ
	agg = sum(data.consumption_share .* ds_eff .* p .^ (1 - σ))
	consumption = (1 .- data.saving_rate) .*
				  (data.consumption_share .* ds_eff .* E .* p .^ (-σ)) ./ agg
	base_consumption = data.household_baseline
	real_gdp_index = all(>=(0), consumption) ?
		tornqvist_quantity_index(p, consumption, ones(N), base_consumption) : NaN
	nominal_gdp = w * sum(L_i)
	return Solution(p, q, wages, consumption, numeraire, real_gdp_index,
					nominal_gdp, model)
end

export fixed_point_init, default_init, init_battery, residual_canary, canary_ok,
	   solve_robust, ipopt_residual_solve, solution_from_x, RobustResult,
	   attempt_newton, attempt_lm, attempt_projected_newton
