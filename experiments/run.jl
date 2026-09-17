# experiments/run.jl — single experiment entry point (Phase 3, ADR-0006)
#
# Usage (from the repo root, with `julia --project=.`):
#
#   julia --project=. experiments/run.jl --list <design>
#   julia --project=. experiments/run.jl --preregister <design> [--actor NAME]
#   julia --project=. experiments/run.jl --design <design> [--cell <run_id>]
#       [--cells a,b,c] [--runs-dir DIR] [--budget-seconds N] [--actor NAME]
#
# - `--list`: print the design's cells with pinned parameters, exit.
# - `--preregister`: write/update `registry/preregistration.toml` from
#   `experiments/designs/<design>.toml`; this is the *only* way that file is
#   written.
# - `--design`: execute cells in file order. Refuses to start (before creating
#   any run dir) when the design file is not preregistered with a matching
#   SHA-256. The design's reference continuation is computed once per batch
#   and its final solution warm-starts all mobile cells and anchors `real_gdp`.
#
# Run layout: `runs/<run_id>/manifest.toml` (written `running` at cell start,
# `executed`/`failed` at the end), `log.txt` (appended progress lines),
# `solution.csv` (`sector,price,quantity`; only on success). An existing run
# dir is never overwritten. `runs/index.csv` (header
# `run_id,date,design,closures,status,gate_summary,headline_metrics,commit`)
# is rewritten on every status transition, sorted by `run_id`.
#
# The file can be `include`d by tests without executing `main()` (guard at
# the bottom). Every function that touches the disk takes an explicit `root`
# (repo root) and/or `runs_dir` keyword so tests can run against temporary
# directories. Public API: `load_design`, `cell_order`, `design_sha256`,
# `list_cells`, `preregister_design`, `preregistration_status`,
# `programme_vectors`, `build_reference`, `build_cell_model`,
# `evaluate_gates`, `execute_cell`, `run_design`, `main`.
#
# Headless-safe: never loads GLMakie (only `BeyondHulten` plus stdlibs and
# the kernel's CSV/DataFrames deps).

using BeyondHulten
using TOML
using SHA
using Dates
using Printf
using LinearAlgebra
using CSV
using DataFrames

const RUN_SEED = 1234
const RUN_SCHEMA_VERSION = 1
const RUN_IO_TABLE = "I-O_DE2019_formatiert.csv"
const RUN_INDEX_HEADER = "run_id,date,design,closures,status,gate_summary,headline_metrics,commit"

"""Default repo root (parent of `experiments/`)."""
default_root() = normpath(joinpath(@__DIR__, ".."))

"""Default runs directory for a root."""
default_runs_dir(root::AbstractString) = joinpath(root, "runs")

"""Default actor: the `USER` environment variable, else `"unknown"`."""
default_actor() = get(ENV, "USER", "unknown")

# ── Small utilities ────────────────────────────────────────────────────

"""SHA-256 hex digest of a file's bytes."""
function sha256_file(path::AbstractString)::String
    return bytes2hex(SHA.sha256(read(path)))
end

"""`true` when `path` exists; otherwise record `"absent"` (never throw)."""
sha256_or_absent(path::AbstractString) = isfile(path) ? sha256_file(path) : "absent"

function _git(root::AbstractString, args::String...)::Union{String,Nothing}
    try
        return strip(read(`git -C $root $args`, String))
    catch
        return nothing
    end
end

"""Recorded HEAD commit, or `"unknown"` outside git."""
git_commit(root::AbstractString) = something(_git(root, "rev-parse", "HEAD"), "unknown")

"""`true` when the tree is dirty or git is unavailable (conservative)."""
function git_dirty(root::AbstractString)::Bool
    porch = _git(root, "status", "--porcelain")
    porch === nothing && return true
    return !isempty(porch)
end

"""Today's ISO date (`yyyy-mm-dd`)."""
iso_date() = string(Dates.today())

"""Current UTC timestamp (`yyyy-mm-ddTHH:MM:SS`)."""
iso_timestamp() = Dates.format(Dates.now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS")

# ── Design files ───────────────────────────────────────────────────────

"""Absolute path of a design file."""
design_path(design::AbstractString; root::AbstractString = default_root()) =
    joinpath(root, "experiments", "designs", design * ".toml")

"""Parse `experiments/designs/<design>.toml` (throws on missing file)."""
function load_design(design::AbstractString; root::AbstractString = default_root())::Dict{String,Any}
    path = design_path(design; root = root)
    isfile(path) || throw(ArgumentError("design file not found: $path"))
    return TOML.parsefile(path)
end

"""Cell ids of a design file in *file order* (TOML tables are unordered)."""
function cell_order(design::AbstractString; root::AbstractString = default_root())::Vector{String}
    path = design_path(design; root = root)
    return cell_order_file(path)
end

"""Cell ids of a design file path in file order."""
function cell_order_file(path::AbstractString)::Vector{String}
    ids = String[]
    for line in eachline(path)
        m = match(r"^\s*\[cells\.([^\]]+)\]\s*$", line)
        m === nothing || push!(ids, strip(m.captures[1]))
    end
    return ids
end

"""SHA-256 hex of the design file's bytes."""
design_sha256(design::AbstractString; root::AbstractString = default_root()) =
    sha256_file(design_path(design; root = root))

"""Fetch a cell table (throws on unknown `run_id`)."""
function design_cell(design_d::Dict{String,Any}, run_id::AbstractString)::Dict{String,Any}
    cells = get(design_d, "cells", Dict{String,Any}())
    haskey(cells, run_id) || throw(ArgumentError("unknown cell \"$run_id\" for this design"))
    return cells[run_id]
end

"""One-line summary of a cell's pinned parameters."""
function cell_summary(run_id::AbstractString, cell::Dict{String,Any})::String
    parts = String["labor=$(cell["labor"])", "financing=$(cell["financing"])"]
    for k in ("eta", "eta_s", "theta", "epsilon", "sigma", "delta_epsilon")
        haskey(cell, k) && push!(parts, "$k=$(cell[k])")
    end
    return "$run_id: " * join(parts, " ") * " — " * get(cell, "note", "")
end

"""Print a design's cells with pinned parameters (implements `--list`)."""
function list_cells(design::AbstractString; root::AbstractString = default_root(),
        io::IO = stdout)::Nothing
    design_d = load_design(design; root = root)
    println(io, "design $(design_d["design"]) — $(get(design_d, "description", ""))")
    for id in cell_order(design; root = root)
        println(io, cell_summary(id, design_cell(design_d, id)))
    end
    return nothing
end

# ── Preregistration ────────────────────────────────────────────────────

"""Absolute path of the preregistration record."""
prereg_path(; root::AbstractString = default_root()) =
    joinpath(root, "registry", "preregistration.toml")

"""Load the preregistration record (`Dict()` when absent)."""
function load_prereg(; root::AbstractString = default_root())::Dict{String,Any}
    path = prereg_path(; root = root)
    isfile(path) || return Dict{String,Any}()
    return TOML.parsefile(path)
end

"""Key scalar parameters pinned by a design (recorded informationally)."""
function pinned_params(design_d::Dict{String,Any})::Dict{String,Any}
    prog = get(design_d, "programme", Dict{String,Any}())
    ref = get(design_d, "reference", Dict{String,Any}())
    return Dict{String,Any}(
        "theta" => get(ref, "theta", 0.0),
        "epsilon" => get(ref, "epsilon", 0.0),
        "sigma" => get(ref, "sigma", 0.0),
        "eta" => get(ref, "eta", 0.0),
        "exo_scale_steps" => get(ref, "exo_scale_steps", 0),
        "total_eur_m" => get(prog, "total_eur_m", 0.0),
        "programme_year" => get(prog, "year", 0),
        "f1_shift" => get(prog, "f1_shift", ""),
    )
end

"""
Write/update `registry/preregistration.toml` from the design file (the
*only* way that file is written; implements `--preregister`). Returns the
recorded SHA-256.
"""
function preregister_design(design::AbstractString; root::AbstractString = default_root(),
        actor::AbstractString = default_actor())::String
    design_d = load_design(design; root = root)
    sha = design_sha256(design; root = root)
    preg = load_prereg(; root = root)
    get!(preg, "schema_version", RUN_SCHEMA_VERSION)
    designs = get!(preg, "designs", Dict{String,Any}())
    designs[design] = Dict{String,Any}(
        "design_sha256" => sha,
        "registered_at" => iso_timestamp(),
        "git_commit" => git_commit(root),
        "actor" => actor,
        "pinned" => pinned_params(design_d),
        "note" => get(design_d, "description", ""),
    )
    open(prereg_path(; root = root), "w") do io
        TOML.print(io, preg)
    end
    return sha
end

"""
Check a design against the preregistration record.
Returns `(true, "")` on match, else `(false, reason)`.
"""
function preregistration_status(design::AbstractString;
        root::AbstractString = default_root())::Tuple{Bool,String}
    preg = load_prereg(; root = root)
    designs = get(preg, "designs", Dict{String,Any}())
    haskey(designs, design) || return (false,
        "design \"$design\" is not preregistered; run `experiments/run.jl --preregister $design` first")
    recorded = get(designs[design], "design_sha256", "")
    current = design_sha256(design; root = root)
    recorded == current || return (false,
        "design \"$design\" changed since preregistration (recorded $(recorded[1:min(8,end)]) vs current $(current[1:8])); re-run `experiments/run.jl --preregister $design`")
    return (true, "")
end

# ── Programme incidence ────────────────────────────────────────────────

"""
Programme incidence `(ψ, g)` for calibrated `data`: renormalized impulse
shares over the kept sectors and `g = (total_eur_m / GDP) .* ψ`.

Smoke designs (`smoke = true`) pin `[programme] explicit = [...]` and skip
the impulses file. `kept` selects kept sector indices for the raw 71-vector.
"""
function programme_vectors(design_d::Dict{String,Any}, data::Data, kept::AbstractVector{<:Integer};
        root::AbstractString = default_root())::Tuple{Vector{Float64},Vector{Float64}}
    prog = design_d["programme"]
    n = length(data.factor_share)
    if get(design_d, "smoke", false) && haskey(prog, "explicit")
        ψ = Float64.(prog["explicit"])
        length(ψ) == n || throw(ArgumentError(
            "explicit programme vector has length $(length(ψ)), expected $n"))
        return ψ, ψ
    end
    total_eur_m = Float64(prog["total_eur_m"])
    slice = prog["column_slice"]
    imp = CSV.read(joinpath(root, prog["source"]), DataFrame)
    col_lo, col_hi = Int(slice[1]), Int(slice[2])
    rows = imp[imp.year .== Int(prog["year"]), :]
    size(rows, 1) >= 1 || throw(ArgumentError("no programme row for year $(prog["year"])"))
    raw = Matrix{Float64}(rows[1:1, col_lo:col_hi])[:]
    v = raw[kept]
    sum(v) > 0 || throw(ArgumentError("programme incidence has zero mass over kept sectors"))
    ψ = v ./ sum(v)
    g = (total_eur_m / data.gdp_production) .* ψ
    return ψ, g
end

# ── Reference continuation (ports cbase2/scripts/verify_v3.jl 18–77) ───

"""
Reference continuation for a design: `read_data` → `retained_dataset` →
bisect `exo_scale` for the smallest scale with `saving_rate ≥ 0` → loop
`exo_scale` steps × the θ ladder with
`mobile_labor_model(data, shocks, θ, ϵ, σ, η)` (NoFinancing), warm-starting
each solve (first init from the linear fixed point
`y0 = (I − Gk) \\ b`), stopping a ladder early when `max|p−1| > 10`.

`design.theta` must equal the final ladder value. Returns a NamedTuple
with `data` (final, `exo_scale = 1`), `sol`, `init_warm = [p; q; w]`,
`resid`, `w_star`, `max_p_dev`, and the `S = I + X − M` canary.

Smoke path: with `data` given, the calibration (read/drop/bisect/loop) is
skipped and a single θ ladder runs on the passed `Data` (tests use this
with an explicit `[programme]` vector).
"""
function build_reference(design_d::Dict{String,Any}; root::AbstractString = default_root(),
        data::Union{Data,Nothing} = nothing)::NamedTuple
    ref = design_d["reference"]
    eta = Float64(ref["eta"])
    epsilon = Float64(ref["epsilon"])
    sigma = Float64(ref["sigma"])
    thetas = Float64.(ref["thetas"])
    Float64(ref["theta"]) == thetas[end] || throw(ArgumentError(
        "reference theta ($(ref["theta"])) must equal the final ladder value ($(thetas[end])))"))

    if data !== nothing
        N = length(data.factor_share)
        shocks = Shocks(ones(N), ones(N), zeros(N))
        init_warm, ref_sol = _theta_ladder(data, shocks, thetas, epsilon, sigma, eta, nothing)
        return _reference_result(data, ref_sol, init_warm)
    end

    dat = design_d["data"]
    drops = Int.(dat["drops"])
    K = Int(ref["exo_scale_steps"])

    data_full = read_data(RUN_IO_TABLE; datadir = root)
    data_v1 = retained_dataset(data_full, drops)
    N = length(data_v1.factor_share)
    shocks = Shocks(ones(N), ones(N), zeros(N))

    # Bisection: smallest exo_scale with saving_rate ≥ 0 (s is monotone).
    s_of(esc) = recalibrate_open(data_v1; exo_scale = esc).saving_rate
    s_of(1.0) > 0 || throw(ErrorException("reference continuation: saving rate at exo_scale = 1 is not positive"))
    lo, hi = 0.0, 1.0
    for _ in 1:40
        mid = (lo + hi) / 2
        s_of(mid) < 0 ? (lo = mid) : (hi = mid)
    end
    esc0 = hi

    init_warm = nothing
    data_cal = data_v1
    ref_sol = nothing
    t0 = time()
    for k in 0:K
        exo_scale = esc0 + (1.0 - esc0) * k / K
        data_cal = recalibrate_open(data_v1; exo_scale = exo_scale)
        init_warm, ref_sol = _theta_ladder(data_cal, shocks, thetas, epsilon, sigma, eta, init_warm)
    end
    @info "reference continuation done" seconds = round(time() - t0; digits = 1)
    return _reference_result(data_cal, ref_sol, init_warm)
end

"""One θ ladder over fixed `(ϵ, σ, η)`; returns `(init_warm, last_sol)`."""
function _theta_ladder(data::Data, shocks::Shocks, thetas::Vector{Float64},
        epsilon::Float64, sigma::Float64, eta::Float64,
        init_warm::Union{Vector{Float64},Nothing})::Tuple{Vector{Float64},Solution}
    N = length(data.factor_share)
    ref_sol = nothing
    for θ in thetas
        ref = mobile_labor_model(data, shocks, θ, epsilon, sigma, eta)
        if init_warm === nothing
            # Very first solve: warm start from the linear fixed point (exact
            # on real calibrations; singular toy fixtures fall back to λ).
            (; Ω_raw, factor_share, consumption_share, import_margin,
                gov_demand, exo_demand, exports_demand) = data
            init_warm = try
                M0 = Ω_raw' * Diagonal(1.0 .- factor_share)
                Gk = M0 + Diagonal(1.0 .- import_margin) * consumption_share *
                    factor_share' * (1 - data.saving_rate) * (1 - sum(gov_demand))
                y0 = (I - Gk) \ ((1.0 .- import_margin) .* (gov_demand .+ exo_demand) .+ exports_demand)
                [ones(N); y0; 1.0]
            catch
                [ones(N); data.λ; 1.0]
            end
        end
        ref_sol = solve(ref; init = init_warm)
        mxp = maximum(abs.(ref_sol.prices_raw .- 1))
        init_warm = [ref_sol.prices_raw; ref_sol.quantities; ref_sol.wages_raw[1]]
        @info "reference rung" theta = θ resid =
            maximum(abs, equilibrium_residuals(ref, init_warm)) w_star = ref_sol.wages_raw[1] max_p_dev = mxp
        mxp > 10 && (@warn "exploded branch; stopping the θ ladder"; break)
    end
    return init_warm, ref_sol
end

"""
    assert_external_canary(model, sol)

Assert the review-2.1 canary identity at a mobile η = 1 solution: the omitted
N-th market residual (`market_clearing_residuals`) must equal the
external-account imbalance `S − (I+X−M)` (`external_balance_canary`). Fixed
closures enforce all N clearings and η = 0 additionally carries the
fixed-allocation gap, so both are exempt.
"""
function assert_external_canary(model::Model, sol::Solution)
    labor_closure(model.options) isa FixedWageClosure && return nothing
    model.options.elasticities.η == 1.0 || return nothing
    X = [sol.prices_raw; sol.quantities; sol.wages_raw[1]]
    market = dot(sol.prices_raw, market_clearing_residuals(model, X))
    canary = external_balance_canary(model, X)
    isapprox(market, canary.diff; atol = 1e-6) || error(
        "external-account canary mismatch: omitted market residual = $market, " *
        "S − (I+X−M) = $(canary.diff) (review finding 2.1, ADR-0010)")
    return nothing
end

"""Final reference solution plus the asserted `S = I + X − M` canary."""
function _reference_result(data::Data, ref_sol::Solution, init_warm::Vector{Float64})::NamedTuple
    p = ref_sol.prices_raw
    w = ref_sol.wages_raw[1]
    L = sum(sectoral_labor_demand(p, ref_sol.quantities, w, ref_sol.model))
    X = [p; ref_sol.quantities; w]
    canary = external_balance_canary(ref_sol.model, X)
    assert_external_canary(ref_sol.model, ref_sol)
    return (data = data, sol = ref_sol, init_warm = init_warm,
        resid = maximum(abs, equilibrium_residuals(ref_sol.model, X)),
        w_star = w, employment = L, max_p_dev = maximum(abs.(p .- 1)),
        canary_s = canary.S, canary_ixm = canary.IX - canary.M,
        canary_diff = canary.diff)
end

# ── Cell construction ──────────────────────────────────────────────────

"""F1 preference-tilt weights from verify_v3.jl (defect 9 fix, commit `0f33ad6`).

Restrict programme shares `ψ` to sectors with positive baseline household
demand `c0` (a zero category cannot be tilted), renormalise to sum 1, then
`d = 1 .+ G0 .* ψ1 ./ max.(c0, 1e-12)` with `G0 = sum(g)`. Throws an
`ArgumentError` when no baseline entry is positive.
"""
function f1_tilt_weights(baseline::AbstractVector, ψ::AbstractVector, g::AbstractVector)::Vector{Float64}
    n = length(baseline)
    length(ψ) == n && length(g) == n || throw(ArgumentError(
        "f1_tilt_weights: length mismatch (baseline=$n, ψ=$(length(ψ)), g=$(length(g)))"))
    all(isfinite, baseline) || throw(ArgumentError("f1_tilt_weights: baseline must be finite"))
    all(isfinite, ψ) || throw(ArgumentError("f1_tilt_weights: ψ must be finite"))
    all(isfinite, g) || throw(ArgumentError("f1_tilt_weights: g must be finite"))
    all(>=(0), ψ) || throw(ArgumentError("f1_tilt_weights: ψ entries must be nonnegative"))
    all(>=(0), g) || throw(ArgumentError("f1_tilt_weights: g entries must be nonnegative"))
    pos = baseline .> 0
    any(pos) || throw(ArgumentError("f1_tilt_weights: no baseline entry is positive"))
    ψ1 = ψ .* pos
    s = sum(ψ1)
    s > 0 || throw(ArgumentError("f1_tilt_weights: programme has zero mass over positive-baseline sectors"))
    ψ1 = ψ1 ./ s
    G0 = sum(g)
    return 1.0 .+ G0 .* ψ1 ./ max.(baseline, 1e-12)
end

"""Financing closure for a cell from `(ψ, g)` (`f1_shift = "tilt_g0_over_c0"`)."""
function cell_financing(fin::AbstractString, ψ::Vector{Float64}, g::Vector{Float64},
        f1_shift::AbstractString, data::Data)::AbstractFinancing
    fin == "F1" || return fin == "F2" ? TaxFinanced(g) : fin == "F3" ? ExternalDebt(g) :
        throw(ArgumentError("unknown financing id $fin"))
    if f1_shift == "one_plus_psi"
        throw(ArgumentError("retired f1_shift \"one_plus_psi\": the non-reference 1 .+ ψ tilt adds ψ directly instead of G0 .* ψ1 ./ c0; use \"tilt_g0_over_c0\""))
    end
    f1_shift == "tilt_g0_over_c0" || throw(ArgumentError("unknown f1_shift \"$f1_shift\""))
    return PreferenceReallocation(f1_tilt_weights(data.household_baseline, ψ, g))
end

"""
Equilibrium model for a design cell (no solve). BF/ALPHA: mobile at η;
BETA: `:beta` via `mobile_labor_model` (solved with `solve_beta`);
GAMMA: `:fixed`; DELTA: `delta_model` at `delta_epsilon`. Shocks are always
`Shocks(ones(N), ones(N), zeros(N))`.
"""
function build_cell_model(cell::Dict{String,Any}, design_d::Dict{String,Any},
        data::Data, ψ::Vector{Float64}, g::Vector{Float64})::Model
    N = length(data.factor_share)
    shocks = Shocks(ones(N), ones(N), zeros(N))
    labor, fin_id = cell["labor"], cell["financing"]
    fin = cell_financing(fin_id, ψ, g, design_d["programme"]["f1_shift"], data)
    if labor == "BF" || labor == "ALPHA"
        return mobile_labor_model(data, shocks, Float64(cell["theta"]),
            Float64(cell["epsilon"]), Float64(cell["sigma"]), Float64(cell["eta"]);
            financing = fin)
    elseif labor == "BETA"
        return mobile_labor_model(data, shocks, Float64(cell["theta"]),
            Float64(cell["epsilon"]), Float64(cell["sigma"]), Float64(cell["eta"]);
            financing = fin, eta_s = Float64(cell["eta_s"]))
    elseif labor == "GAMMA"
        return mobile_labor_model(data, shocks, Float64(cell["theta"]),
            Float64(cell["epsilon"]), Float64(cell["sigma"]), Float64(cell["eta"]);
            closure = :fixed, financing = fin)
    elseif labor == "DELTA"
        return delta_model(data, shocks; ε = Float64(cell["delta_epsilon"]), financing = fin)
    end
    throw(ArgumentError("unknown labor id $labor"))
end

"""Solve a cell model (BETA via `solve_beta`, else `solve`), warm-started."""
function solve_cell(cell::Dict{String,Any}, design_d::Dict{String,Any}, data::Data,
        ψ::Vector{Float64}, g::Vector{Float64}, init_warm::Vector{Float64})::Solution
    model = build_cell_model(cell, design_d, data, ψ, g)
    if cell["labor"] == "BETA"
        N = length(data.factor_share)
        shocks = Shocks(ones(N), ones(N), zeros(N))
        return solve_beta(data, shocks, Float64(cell["theta"]), Float64(cell["epsilon"]),
            Float64(cell["sigma"]), Float64(cell["eta"]); eta_s = Float64(cell["eta_s"]),
            financing = cell_financing(cell["financing"], ψ, g, design_d["programme"]["f1_shift"], data),
            init = init_warm)
    end
    return solve(model; init = init_warm)
end

# ── Gates / metrics / diagnostics ──────────────────────────────────────

"""Format a gate summary fragment `name value<=tol ok|FAIL`."""
gate_frag(name::AbstractString, value::Real, tol::Real, pass::Bool) =
    @sprintf("%s %.3g<=%.0e %s", name, value, tol, pass ? "ok" : "FAIL")

"""
Evaluate a solved cell: residual / budget / labour-or-wage gates, headline
metrics (against the reference solution), and the `S = I + X − M` canary
plus `external_balance` / `public_budget` diagnostics (never gates).
"""
function evaluate_gates(cell::Dict{String,Any}, design_d::Dict{String,Any},
        model::Model, sol::Solution, ref_sol::Solution)::NamedTuple
    gates = design_d["gates"]
    residual_tol = Float64(gates["residual_tol"])
    budget_tol = Float64(gates["budget_tol"])
    labour_tol = Float64(gates["labour_tol"])
    wage_tol = Float64(gates["wage_tol"])
    data = model.data
    fin = model.financing
    fixed = labor_closure(model.options) isa FixedWageClosure

    p, q = sol.prices_raw, sol.quantities
    w = sol.wages_raw[1]
    X = fixed ? [p; q] : [p; q; w]
    resid = maximum(abs, equilibrium_residuals(model, X))
    L_sum = sum(sectoral_labor_demand(p, q, w, model))
    E = household_expenditure(fin, model, w * L_sum, p, L_sum)
    budget = abs(dot(p, sol.consumption) - (1 - data.saving_rate) * E)
    resid_pass = resid < residual_tol
    budget_pass = budget < budget_tol

    third_name, third_value, third_tol, third_pass =
        # Fixed closures hard-pin w = 1 (wages_raw); the normalized sol.wages
        # folds in the CPI numeraire and must not enter the gate.
        fixed ? ("wage", maximum(abs.(sol.wages_raw .- 1)), wage_tol,
                 maximum(abs.(sol.wages_raw .- 1)) < wage_tol) :
                ("labour", abs(labor_market_residual(labor_closure(model.options), model, L_sum, w)),
                 labour_tol, abs(labor_market_residual(labor_closure(model.options), model, L_sum, w)) < labour_tol)

    rgdp = real_gdp(sol)
    rgdp_ref = real_gdp(ref_sol)
    overall = (resid_pass && budget_pass && third_pass) ? "pass" : "fail"
    summary = gate_frag("resid", resid, residual_tol, resid_pass) * "; " *
        gate_frag("budget", budget, budget_tol, budget_pass) * "; " *
        gate_frag(third_name, third_value, third_tol, third_pass)

    # Value-consistent external-account canary (review finding 2.1, ADR-0010):
    # recorded for every cell and asserted for mobile η = 1 cells, where the
    # omitted N-th market residual must equal S − (I+X−M).
    assert_external_canary(model, sol)
    canary = external_balance_canary(model, [p; q; fixed ? 1.0 : w])
    return (
        gates = Dict{String,Any}(
            "residual" => Dict{String,Any}("value" => resid, "tolerance" => residual_tol, "pass" => resid_pass),
            "budget" => Dict{String,Any}("value" => budget, "tolerance" => budget_tol, "pass" => budget_pass),
            third_name => Dict{String,Any}("value" => third_value, "tolerance" => third_tol, "pass" => third_pass),
            "overall" => overall),
        gate_summary = summary,
        metrics = Dict{String,Any}(
            "real_gdp" => rgdp, "real_gdp_ref" => rgdp_ref,
            "real_gdp_rel" => rgdp / rgdp_ref - 1,
            "employment" => L_sum, "wage" => w,
            "nominal_gdp" => nominal_gdp(sol),
            "max_abs_price_dev" => maximum(abs.(p .- 1))),
        diagnostics = Dict{String,Any}(
            "canary_s" => canary.S, "canary_ixm" => canary.IX - canary.M,
            "canary_diff" => canary.diff,
            "external_balance" => external_balance(fin, model, p),
            "public_budget" => public_budget(fin, model, p)),
    )
end

# ── Manifest / index / scenarios ───────────────────────────────────────

"""Provenance block shared by all manifests of a batch."""
function batch_provenance(design::AbstractString, design_d::Dict{String,Any};
        root::AbstractString = default_root())::Dict{String,Any}
    dat = design_d["data"]
    prog = design_d["programme"]
    HEAD = git_commit(root)
    # Smoke designs pin an explicit programme vector instead of a source
    # file; absent keys simply contribute no hashes.
    paths = String[
        "experiments/designs/" * design * ".toml",
        "data/" * RUN_IO_TABLE,
        String.(dat["calibration_artifacts"])...,
    ]
    haskey(prog, "source") && push!(paths, prog["source"])
    data_sha = Dict{String,Any}(p => sha256_or_absent(joinpath(root, p)) for p in paths)
    return Dict{String,Any}(
        "git_commit" => HEAD,
        "git_dirty" => git_dirty(root),
        "julia_version" => string(VERSION),
        "manifest_sha256" => sha256_or_absent(joinpath(root, "Manifest.toml")),
        "design_sha256" => design_sha256(design; root = root),
        "data_sha256" => data_sha,
        "seed" => RUN_SEED,
    )
end

"""Append a line to a run's `log.txt` (and echo to stdout)."""
function run_log(rundir::AbstractString, line::AbstractString)::Nothing
    open(joinpath(rundir, "log.txt"), "a") do io
        println(io, iso_timestamp() * " " * line)
    end
    println(line)
    return nothing
end

"""Write `solution.csv` (`sector,price,quantity`; raw prices)."""
function write_solution(rundir::AbstractString, data::Data, sol::Solution)::Nothing
    open(joinpath(rundir, "solution.csv"), "w") do io
        println(io, "sector,price,quantity")
        for i in 1:length(sol.quantities)
            println(io, "\"", replace(string(data.io.Sektoren[i]), "\"" => "\"\""), "\"",
                ",", sol.prices_raw[i], ",", sol.quantities[i])
        end
    end
    return nothing
end

"""Minimal CSV row writer (quotes fields containing `,\"` or newlines)."""
csv_field(s::AbstractString) =
    (occursin(r"[,\"\n]", s) ? "\"" * replace(s, "\"" => "\"\"") * "\"" : s)

"""Rewrite `runs/index.csv` with one row per run dir, sorted by `run_id`."""
function rewrite_index(; runs_dir::AbstractString)::Nothing
    rows = Dict{String,Vector{String}}()
    isdir(runs_dir) || mkpath(runs_dir)
    for entry in readdir(runs_dir)
        mandir = joinpath(runs_dir, entry)
        manpath = joinpath(mandir, "manifest.toml")
        isdir(mandir) && isfile(manpath) || continue
        man = TOML.parsefile(manpath)
        scen = get(man, "scenario", Dict{String,Any}())
        gt = get(man, "gates", Dict{String,Any}())
        me = get(man, "metrics", Dict{String,Any}())
        closures = string(get(scen, "labour", ""), "+", get(scen, "financing", ""))
        gs = if haskey(gt, "overall")
            parts = String[]
            for (k, v) in gt
                k == "overall" && continue
                v isa Dict || continue
                push!(parts, gate_frag(k == "residual" ? "resid" : k,
                    v["value"], v["tolerance"], v["pass"]))
            end
            isempty(parts) ? "overall=$(gt["overall"])" : join(sort(parts), "; ")
        else
            "overall=$(get(man, "status", "?"))"
        end
        hm = haskey(me, "real_gdp_rel") ?
            @sprintf("gdp_rel=%.6f; L=%.6f; w=%.6f", me["real_gdp_rel"], me["employment"], me["wage"]) : ""
        rows[string(man["run_id"])] = String[
            string(man["run_id"]), string(get(man, "date", "")),
            string(get(man, "design", "")), closures, string(get(man, "status", "")),
            gs, hm, string(get(get(man, "provenance", Dict()), "git_commit", ""))]
    end
    open(joinpath(runs_dir, "index.csv"), "w") do io
        println(io, RUN_INDEX_HEADER)
        for id in sort(collect(keys(rows)))
            println(io, join(csv_field.(rows[id]), ","))
        end
    end
    return nothing
end

"""Update one `registry/scenarios.csv` row from a finished (or started) cell."""
function update_scenario_row(run_id::AbstractString, design::AbstractString,
        cell::Dict{String,Any}, status::AbstractString, commit::AbstractString;
        root::AbstractString = default_root(), note_suffix::AbstractString = "")::Nothing
    scenpath = joinpath(root, "registry", "scenarios.csv")
    # Normalize every column to String: without `stringtype = String`, CSV
    # infers narrow InlineString widths (e.g. String7) from the current cell
    # values and the longer status/evidence/commit assignments below throw
    # `ArgumentError: string too large`. (Found by the smoke-manifest tests.)
    raw = DataFrame(CSV.File(scenpath; stringtype = String))
    df = DataFrame([c => string.(coalesce.(raw[!, c], "")) for c in names(raw)])
    r = findfirst(==(run_id), string.(coalesce.(df.run_id, "")))
    r === nothing && throw(ArgumentError("no scenarios.csv row for run_id \"$run_id\""))
    df[r, :design] = design
    df[r, :status] = status
    df[r, :labor] = string(cell["labor"])
    df[r, :financing] = string(cell["financing"])
    df[r, :eta] = string(cell["eta"])
    df[r, :eta_s] = string(get(cell, "eta_s", 0.0))
    df[r, :theta] = string(cell["theta"])
    df[r, :epsilon] = string(cell["epsilon"])
    df[r, :sigma] = string(cell["sigma"])
    df[r, :shock] = "impulses.csv"
    df[r, :magnitude] = "1.0"
    df[r, :data_vintage] = "cbase2-v3"
    df[r, :evidence] = "runs/$run_id/manifest.toml; runs/$run_id/log.txt"
    df[r, :commit] = commit
    isempty(note_suffix) || (df[r, :notes] = string(coalesce(df[r, :notes], ""), note_suffix))
    CSV.write(scenpath, df)
    return nothing
end

# ── Cell execution ─────────────────────────────────────────────────────

"""
Execute one cell: refuse when the run dir exists (operator must register a
`-v2` variant), else write `manifest.toml` (`running`), solve, evaluate,
and update the manifest, `solution.csv`, `runs/index.csv`, and
`scenarios.csv`. Any per-cell failure becomes a `failed` manifest with an
`[error]` table; the batch continues.
Returns the final status (`"executed"`, `"failed"`, or `"refused"`).
"""
function execute_cell(run_id::AbstractString, design::AbstractString,
        design_d::Dict{String,Any}, data::Data, ψ::Vector{Float64}, g::Vector{Float64},
        ref_sol::Solution, init_warm::Vector{Float64}, prov::Dict{String,Any};
        root::AbstractString = default_root(), runs_dir::AbstractString = default_runs_dir(root),
        actor::AbstractString = default_actor())::String
    cell = design_cell(design_d, run_id)
    rundir = joinpath(runs_dir, run_id)
    if ispath(rundir)
        println("refusing to overwrite existing run dir runs/$run_id/ — " *
            "register a variant run_id ($run_id-v2) for a rerun")
        return "refused"
    end
    mkpath(rundir)
    man = Dict{String,Any}(
        "schema_version" => RUN_SCHEMA_VERSION,
        "run_id" => run_id, "design" => design, "cell" => run_id,
        "status" => "running", "date" => iso_date(), "actor" => actor,
        "provenance" => prov,
        "scenario" => Dict{String,Any}(
            "labour" => string(cell["labor"]), "financing" => string(cell["financing"]),
            "eta" => Float64(cell["eta"]), "eta_s" => Float64(get(cell, "eta_s", 0.0)),
            "theta" => Float64(cell["theta"]), "epsilon" => Float64(cell["epsilon"]),
            "sigma" => Float64(cell["sigma"]),
            "shock" => "impulses.csv", "magnitude" => 1.0),
        "solver" => Dict{String,Any}(
            "init" => "warm:reference-final", "reference" => "batch-continuation",
            "algorithm" => "kernel solve / solve_beta (Newton + residual-gated LM polish)"),
        "artifacts" => Dict{String,Any}("log" => "log.txt", "solution" => ""),
    )
    if haskey(cell, "delta_epsilon")
        man["scenario"]["delta_epsilon"] = Float64(cell["delta_epsilon"])
    end
    open(joinpath(rundir, "manifest.toml"), "w") do io
        TOML.print(io, man)
    end
    run_log(rundir, "start $run_id design=$design actor=$actor commit=$(prov["git_commit"])")
    run_log(rundir, "reference real_gdp_ref=$(real_gdp(ref_sol))")
    update_scenario_row(run_id, design, cell, "running", prov["git_commit"]; root = root)
    try
        sol = solve_cell(cell, design_d, data, ψ, g, init_warm)
        ev = evaluate_gates(cell, design_d, sol.model, sol, ref_sol)
        for (k, v) in ev.metrics
            isfinite(v) || throw(ErrorException("non-finite metric $k in $run_id"))
        end
        status = ev.gates["overall"] == "pass" ? "executed" : "failed"
        man["status"] = status
        man["gates"] = ev.gates
        man["metrics"] = ev.metrics
        man["diagnostics"] = ev.diagnostics
        if status == "executed"
            write_solution(rundir, data, sol)
            man["artifacts"]["solution"] = "solution.csv"
        end
        open(joinpath(rundir, "manifest.toml"), "w") do io
            TOML.print(io, man)
        end
        run_log(rundir, "$(status): $(ev.gate_summary)")
        rewrite_index(; runs_dir = runs_dir)
        update_scenario_row(run_id, design, cell, status, prov["git_commit"]; root = root,
            note_suffix = " | $(status) $(iso_date()) (see runs/$run_id/manifest.toml)")
        return status
    catch e
        # Any per-cell failure — solve, gates, metrics, or writes — becomes
        # a `failed` manifest; the batch continues.
        man["status"] = "failed"
        man["gates"] = Dict{String,Any}("overall" => "fail")
        man["error"] = Dict{String,Any}("type" => string(typeof(e)), "message" => sprint(showerror, e))
        try
            open(joinpath(rundir, "manifest.toml"), "w") do io
                TOML.print(io, man)
            end
        catch io_e
            println("could not write failed manifest for $run_id: $(sprint(showerror, io_e))")
        end
        run_log(rundir, "failed: $(typeof(e)): $(sprint(showerror, e))")
        rewrite_index(; runs_dir = runs_dir)
        update_scenario_row(run_id, design, cell, "failed", prov["git_commit"]; root = root,
            note_suffix = " | failed $(iso_date()) (see runs/$run_id/manifest.toml)")
        return "failed"
    end
end

# ── Batch driver ───────────────────────────────────────────────────────

"""
Execute a design's cells in file order (implements `--design`).

Refuses to start — before creating any run dir — when the design file is
not preregistered with a matching SHA-256. The reference continuation is
computed once per batch; its final solution warm-starts every cell and
anchors `real_gdp`. `--cell`/`--cells` select a subset; `--budget-seconds`
stops the batch cleanly before a cell whose start would exceed the budget.
Per-cell exceptions are caught by `execute_cell`; the batch continues.
"""
function run_design(design::AbstractString; root::AbstractString = default_root(),
        runs_dir::AbstractString = default_runs_dir(root),
        cell::Union{AbstractString,Nothing} = nothing,
        cells::Union{AbstractString,Vector{<:AbstractString},Nothing} = nothing,
        budget_seconds::Union{Real,Nothing} = nothing,
        actor::AbstractString = default_actor(),
        data::Union{Data,Nothing} = nothing)::Dict{String,String}
    ok, reason = preregistration_status(design; root = root)
    ok || throw(ErrorException("refusing to run: $reason"))
    design_d = load_design(design; root = root)
    order = cell_order(design; root = root)
    wanted = if cell !== nothing && cells !== nothing
        throw(ArgumentError("pass at most one of --cell and --cells"))
    elseif cell !== nothing
        [cell]
    elseif cells !== nothing
        cells isa AbstractString ? String.(split(cells, ",")) : String.(cells)
    else
        order
    end
    for id in wanted
        id in order || throw(ArgumentError("unknown cell \"$id\" for design \"$design\""))
    end

    ref = build_reference(design_d; root = root, data = data)
    drops = Int.(design_d["data"]["drops"])
    kept = if data !== nothing
        collect(1:length(ref.data.factor_share))
    else
        n_full = length(ref.data.factor_share) + length(drops)
        sort(setdiff(1:n_full, drops))
    end
    ψ, g = programme_vectors(design_d, ref.data, kept; root = root)
    prov = batch_provenance(design, design_d; root = root)

    t0 = time()
    results = Dict{String,String}()
    for id in order
        id in wanted || continue
        if budget_seconds !== nothing && (time() - t0) > Float64(budget_seconds)
            println("budget exceeded ($(round(time()-t0; digits=1))s > $budget_seconds s); " *
                "stopping cleanly before $id")
            break
        end
        results[id] = try
            execute_cell(id, design, design_d, ref.data, ψ, g,
                ref.sol, ref.init_warm, prov; root = root, runs_dir = runs_dir, actor = actor)
        catch e
            # Infrastructure failures outside the cell itself (already
            # manifest-backed inside execute_cell) must not abort the batch.
            println("cell $id errored outside its manifest: $(sprint(showerror, e))")
            "error"
        end
    end
    return results
end

# ── CLI ────────────────────────────────────────────────────────────────

"""Print usage and exit with code 1."""
function usage()::Nothing
    println(stderr, """usage:
  julia --project=. experiments/run.jl --list <design>
  julia --project=. experiments/run.jl --preregister <design> [--actor NAME]
  julia --project=. experiments/run.jl --design <design> [--cell <run_id>] [--cells a,b,c] [--runs-dir DIR] [--budget-seconds N] [--actor NAME]""")
    exit(1)
end

"""CLI entry point (only via the `PROGRAM_FILE` guard below)."""
function main(args::Vector{String} = ARGS; root::AbstractString = default_root())::Nothing
    isempty(args) && usage()
    mode, design = args[1], length(args) >= 2 && !startswith(args[2], "--") ? args[2] : nothing
    opts = Dict{String,String}()
    i = 3
    while i <= length(args)
        a = args[i]
        startswith(a, "--") || usage()
        key = a[3:end]
        key in ("cell", "cells", "runs-dir", "budget-seconds", "actor") ||
            (println(stderr, "unknown option $a"); usage())
        i + 1 > length(args) && usage()
        opts[key] = args[i+1]
        i += 2
    end
    if mode == "--list"
        design === nothing && usage()
        list_cells(design; root = root)
    elseif mode == "--preregister"
        design === nothing && usage()
        sha = preregister_design(design; root = root, actor = get(opts, "actor", default_actor()))
        println("preregistered $design (sha256=$sha)")
    elseif mode == "--design"
        design === nothing && usage()
        budget = haskey(opts, "budget-seconds") ? parse(Float64, opts["budget-seconds"]) : nothing
        results = run_design(design; root = root,
            runs_dir = get(opts, "runs-dir", default_runs_dir(root)),
            cell = get(opts, "cell", nothing), cells = get(opts, "cells", nothing),
            budget_seconds = budget, actor = get(opts, "actor", default_actor()))
        for id in sort(collect(keys(results)))
            println("$id: $(results[id])")
        end
    else
        usage()
    end
    return nothing
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
