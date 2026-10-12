# scripts/reproduce.jl — reproduction baseline for the paper pipeline (WP0)
#
# Usage (from the repository root):
#
#   julia --project=. scripts/reproduce.jl --write-baseline
#   julia --project=. scripts/reproduce.jl --check
#
# What it does, in both modes: re-solves EVERY cell of the pinned paper
# design `matrix_5x3_v10` (33 cells) with the harness (`solve_cell` +
# `evaluate_gates` from `experiments/run.jl`) and hard-validates each cell
# against its recorded run artifacts (`runs/<id>/manifest.toml` +
# `solution.csv`) via `BeyondHulten.validate_cell`, reporting per-cell
# PASS/FAIL with the equilibrium residual and the max abs deviation from
# the stored artifacts.
#
# --write-baseline additionally:
#   * dumps full-precision per-cell metrics to `baseline/v10_metrics.csv`
#     (17 significant digits, stable row order: run_id, then metric /
#     manifest-metric / gate blocks sorted by key, then the price /
#     quantity / wage vectors in sector order);
#   * regenerates the tidy output CSVs with the documented generator
#     (`experiments/plot_matrix.jl --design matrix_5x3_v10 --no-figures
#     --data-dir output`, invoked as `plot_main`) and copies them into
#     `baseline/output/` (the generator writes round-trip-precision floats;
#     full-precision copies are written alongside iff a round-trip check
#     finds rounded values);
#   * copies `output/variance_decomposition_sobol.csv` into
#     `baseline/output/` as-is (its generator is the legacy root-level
#     `rerun_results.jl`; see `baseline/README.md` for its provenance).
#
# --check repeats the re-solve and the generator and compares the result
# against `baseline/` within 1e-8 (abs + rel, NaN == NaN), listing every
# mismatch and exiting non-zero on drift or on any failed cell.
#
# The tool is read-only outside `baseline/` and the gitignored `output/`
# scratch area: it never writes `runs/`, `registry/`, or any tracked state.
# Reused machinery: `experiments/plot_matrix.jl` (which includes
# `experiments/run.jl`).

using BeyondHulten
using TOML
using CSV
using DataFrames
using Printf

# run.jl is include-safe (PROGRAM_FILE guard); plot_matrix.jl includes it when
# it is not already loaded and only runs `plot_main()` when it is the program
# file, so a plain include here is safe.
include(joinpath(@__DIR__, "..", "experiments", "plot_matrix.jl"))

const REPRO_DESIGN = "matrix_5x3_v10"
const REPRO_TOL = 1e-8
const METRIC_COLUMNS = ["run_id", "series", "key", "idx", "label", "value"]

"""Full-precision Float64 text: 17 significant digits guarantee round-trip."""
fmt17(x::Real)::String = @sprintf("%.17g", Float64(x))

function reproduce_usage()::Int
    println(stderr, """usage:
  julia --project=. scripts/reproduce.jl --write-baseline
  julia --project=. scripts/reproduce.jl --check""")
    return 2
end

# ── Re-solve + hard validation of every cell ────────────────────────────

"""
Re-solve every cell of `matrix_5x3_v10` in sorted run-id order and
hard-validate against `runs/<id>/manifest.toml` + `solution.csv`.
Returns `(design_d, sectors, recs)` where `recs` is a vector of
`Dict{String,Any}` records with the keys `id`, `status`, `verdict`,
`reason`, `deltas`, `skipped`, `metrics`, `man_metrics`, `gates`,
`man_gates`, `prices`, `quantities`, `wages`, `resid`. A solve/gate/read
failure is recorded as a failing cell and never aborts the batch.
"""
function resolve_cells(root::AbstractString)
    ok, reason = preregistration_status(REPRO_DESIGN; root = root)
    ok || error("refusing to run: $reason")
    design_d = load_design(REPRO_DESIGN; root = root)
    ids = sort(String.(cell_order(REPRO_DESIGN; root = root)))

    println("building reference continuation for $REPRO_DESIGN ...")
    ref = build_reference(design_d; root = root)
    drops = Vector{Int}(design_d["data"]["drops"])
    n_full = length(ref.data.factor_share) + length(drops)
    kept = sort(setdiff(1:n_full, drops))
    ψ, g = programme_vectors(design_d, ref.data, kept; root = root)
    n = length(ref.data.factor_share)
    sectors = String[string(x) for x in ref.data.io.Sektoren[1:n]]

    recs = Dict{String,Any}[]
    for id in ids
        rec = Dict{String,Any}(
            "id" => id, "status" => "n/a", "verdict" => "FAIL", "reason" => "",
            "deltas" => Dict{String,Float64}(), "skipped" => String[],
            "metrics" => Dict{String,Any}(), "man_metrics" => Dict{String,Any}(),
            "gates" => Dict{String,Any}(), "man_gates" => Dict{String,Any}(),
            "prices" => Float64[], "quantities" => Float64[],
            "wages" => Float64[], "resid" => NaN)
        try
            println("re-solving $id ...")
            cell = design_cell(design_d, id)
            sol = solve_cell(cell, design_d, ref.data, ψ, g, ref.init_warm)
            ev = evaluate_gates(cell, design_d, sol.model, sol, ref.sol)
            ev.gates["overall"] == "pass" || throw(ErrorException(
                "re-solved gates report \"$(ev.gates["overall"])\""))
            manpath = joinpath(root, "runs", id, "manifest.toml")
            isfile(manpath) || throw(ErrorException(
                "recorded manifest not found: runs/$id/manifest.toml"))
            manifest = TOML.parsefile(manpath)
            stored_prices, stored_quantities =
                read_stored_solution(joinpath(root, "runs", id, "solution.csv"))
            vc = validate_cell(manifest, ev.metrics, sol.prices_raw,
                sol.quantities, stored_prices, stored_quantities)
            rec["status"] = string(get(manifest, "status", "?"))
            rec["deltas"] = vc.deltas
            rec["skipped"] = collect(String, vc.skipped)
            rec["metrics"] = ev.metrics
            rec["man_metrics"] = get(manifest, "metrics", Dict{String,Any}())
            rec["gates"] = ev.gates
            rec["man_gates"] = get(manifest, "gates", Dict{String,Any}())
            rec["prices"] = Float64.(collect(sol.prices_raw))
            rec["quantities"] = Float64.(collect(sol.quantities))
            rec["wages"] = Float64.(collect(sol.wages_raw))
            rec["resid"] = Float64(ev.gates["residual"]["value"])
            if vc.pass
                rec["verdict"] = "PASS"
            else
                rec["reason"] = join(vc.failures, "; ")
            end
        catch e
            rec["reason"] = sprint(showerror, e)
        end
        push!(recs, rec)
    end
    return design_d, sectors, recs
end

"""Per-cell validation report: PASS/FAIL with residual and max abs deviations."""
function report_cells(recs)::Int
    println("validation report for $REPRO_DESIGN (vs runs/<id> artifacts):")
    println("run_id | status | resid | max|dp| | max|dq| | max|dmetric| | verdict")
    n_pass = 0
    for rec in recs
        d = rec["deltas"]
        dp = get(d, "price", NaN)
        dq = get(d, "quantity", NaN)
        dm = plot_max_metric_delta(d)
        tail = rec["verdict"] == "PASS" ? "PASS" : "FAIL ($(rec["reason"]))"
        rec["verdict"] == "PASS" && (n_pass += 1)
        println("$(rec["id"]) | $(rec["status"]) | $(@sprintf("%.3g", rec["resid"])) | " *
            "$(@sprintf("%.3g", dp)) | $(@sprintf("%.3g", dq)) | " *
            "$(@sprintf("%.3g", dm)) | $tail")
    end
    println("$n_pass/$(length(recs)) cells validated")
    return n_pass
end

# ── Full-precision metrics rows ─────────────────────────────────────────

"""
One long-format metrics record per row, in stable order. Columns:
`run_id,series,key,idx,label,value`; `series` in `metric` (re-solved
headline metrics), `manifest_metric` (the manifest `[metrics]` block),
`gate` / `manifest_gate` (gate values), `price`, `quantity`, `wage`
(indexed vectors in sector order, `label` = sector name). Scalars carry
`idx = 0` and an empty label. Order: metric blocks sorted by key, then
the vectors in sector order. Values are formatted with 17 significant
digits (`%.17g`), which round-trips every Float64.
"""
function metric_rows(rec::Dict{String,Any}, sectors::Vector{String})::Vector{Vector{String}}
    rows = Vector{String}[]
    id = rec["id"]
    scalar(name, key, v) = push!(rows,
        [id, name, key, "0", "", fmt17(v)])
    for k in sort!(collect(keys(rec["metrics"])))
        scalar("metric", String(k), Float64(rec["metrics"][k]))
    end
    for k in sort!(collect(keys(rec["man_metrics"])))
        scalar("manifest_metric", String(k), Float64(rec["man_metrics"][k]))
    end
    for k in sort!(collect(keys(rec["gates"])))
        String(k) == "overall" && continue
        gv = rec["gates"][k]
        scalar("gate", String(k), Float64(gv isa AbstractDict ? gv["value"] : gv))
    end
    for k in sort!(collect(keys(rec["man_gates"])))
        String(k) == "overall" && continue
        gv = rec["man_gates"][k]
        gv isa AbstractDict || continue
        scalar("manifest_gate", String(k), Float64(gv["value"]))
    end
    for (series, vec) in (("price", rec["prices"]), ("quantity", rec["quantities"]),
            ("wage", rec["wages"]))
        for (i, v) in enumerate(vec)
            lab = length(vec) == length(sectors) ? sectors[i] : ""
            push!(rows, [id, series, "", string(i), lab, fmt17(v)])
        end
    end
    return rows
end

"""All metrics rows of all records, in record order."""
function all_metric_rows(recs, sectors::Vector{String})::Vector{Vector{String}}
    rows = Vector{String}[]
    for rec in recs
        append!(rows, metric_rows(rec, sectors))
    end
    return rows
end

"""Write the metrics rows to `path` (CSV, stable header)."""
function write_metrics_csv(path::AbstractString, rows::Vector{Vector{String}})
    df = DataFrame(
        run_id = [r[1] for r in rows], series = [r[2] for r in rows],
        key = [r[3] for r in rows], idx = [r[4] for r in rows],
        label = [r[5] for r in rows], value = [r[6] for r in rows])
    CSV.write(path, df)
    println("wrote $path ($(length(rows)) rows)")
end

"""`true` when the baseline value `b` and the recomputed value `r` agree
within `REPRO_TOL` abs + rel (NaN == NaN, Inf == Inf)."""
function value_match(b::Float64, r::Float64)::Bool
    (isnan(b) && isnan(r)) && return true
    (isinf(b) || isinf(r)) && return b == r
    return abs(r - b) <= REPRO_TOL + REPRO_TOL * abs(b)
end

"""
Compare freshly built metrics rows against `baseline/v10_metrics.csv`
(row order and keys exact, values within 1e-8 abs + rel). Returns the
list of mismatch descriptions (empty = consistent).
"""
function check_metrics_csv(path::AbstractString, rows::Vector{Vector{String}})::Vector{String}
    problems = String[]
    isfile(path) || return ["baseline metrics file missing: $path"]
    df = CSV.read(path, DataFrame; types = String)
    names(df) == METRIC_COLUMNS || push!(problems,
        "header mismatch: baseline=$(names(df)) expected=$METRIC_COLUMNS")
    if nrow(df) != length(rows)
        push!(problems, "row count mismatch: baseline=$(nrow(df)) recomputed=$(length(rows))")
    end
    n = min(nrow(df), length(rows))
    for i in 1:n
        r = rows[i]
        bkey = (String(df.run_id[i]), String(df.series[i]), String(df.key[i]),
            coalesce(df.idx[i], ""), coalesce(df.label[i], ""))
        rkey = (r[1], r[2], r[3], r[4], r[5])
        if bkey != rkey
            push!(problems, "row $i key mismatch: baseline=$bkey recomputed=$rkey")
            continue
        end
        bv = parse(Float64, df.value[i])
        rv = parse(Float64, r[6])
        value_match(bv, rv) || push!(problems,
            "row $i $(r[1]) $(r[2]):$(r[3]) differs: baseline=$bv recomputed=$rv " *
            "(|d|=$(abs(rv - bv)))")
    end
    return problems
end

# ── Output CSVs (documented generator) ──────────────────────────────────

"""Regenerate the tidy output CSVs with the documented generator
(`experiments/plot_matrix.jl --design matrix_5x3_v10 --no-figures
--data-dir output`); returns the generator's exit code."""
function run_generator(root::AbstractString)::Int
    args = ["--design", REPRO_DESIGN, "--no-figures", "--data-dir",
        joinpath(root, "output")]
    println("running documented generator: experiments/plot_matrix.jl " * join(args, " "))
    return plot_main(args; root = root)
end

"""Names of the generator artifacts captured in `baseline/output/`."""
function output_artifact_names()::Vector{String}
    return [REPRO_DESIGN * "_summary.csv", REPRO_DESIGN * "_sectoral.csv",
        "variance_decomposition_sobol.csv"]
end

"""
Copy the generated CSVs (and the Sobol decomposition CSV, whose generator
is the legacy root-level `rerun_results.jl` — see `baseline/README.md`)
from `output/` into `baseline/output/`.
"""
function copy_outputs(root::AbstractString)
    for name in output_artifact_names()
        src = joinpath(root, "output", name)
        isfile(src) || error("missing generated artifact: output/$name")
        dst = joinpath(root, "baseline", "output", name)
        cp(src, dst; force = true)
        println("copied output/$name -> baseline/output/$name")
    end
end

"""
Round-trip check: every numeric value of the regenerated summary CSV must
parse back to the exact Float64 of the re-solved metric (`CSV.write` writes
shortest round-trip floats). Returns the max abs deviation (0.0 = the
generator writes full precision). Deviations > 0 trigger full-precision
copies `<name>.full-precision.csv` alongside the plain ones.
"""
function summary_roundtrip_check(root::AbstractString, recs)::Float64
    by_id = Dict{String,Dict{String,Any}}(String(r["id"]) => r["metrics"] for r in recs)
    path = joinpath(root, "output", REPRO_DESIGN * "_summary.csv")
    df = CSV.read(path, DataFrame; types = String)
    maxdev = 0.0
    rounded = false
    for i in 1:nrow(df)
        id = String(df.run_id[i])
        haskey(by_id, id) || continue
        m = by_id[id]
        for col in names(df)
            (col in METRIC_COLUMNS || !haskey(m, col)) && continue
            v = parse(Float64, df[i, col])
            dev = v == Float64(m[col]) ? 0.0 :
                (isnan(v) && isnan(Float64(m[col]))) ? 0.0 : abs(v - Float64(m[col]))
            dev > 0 && (rounded = true)
            maxdev = max(maxdev, dev)
        end
    end
    if rounded
        for name in (REPRO_DESIGN * "_summary.csv", REPRO_DESIGN * "_sectoral.csv")
            src = joinpath(root, "output", name)
            df2 = CSV.read(src, DataFrame; types = String)
            for col in names(df2), i in 1:nrow(df2)
                v = tryparse(Float64, df2[i, col])
                v === nothing && continue
                df2[i, col] = fmt17(v)
            end
            dst = joinpath(root, "baseline", "output",
                replace(name, ".csv" => ".full-precision.csv"))
            CSV.write(dst, df2)
            println("wrote $dst (generator output was rounded)")
        end
    end
    println("summary CSV round-trip check: max|dev| = $maxdev" *
        (rounded ? " (rounded values found)" : " (full precision)"))
    return maxdev
end

"""
Compare a regenerated tidy CSV against its `baseline/output/` copy
(headers exact, cells exact for text and within 1e-8 abs + rel for
numbers, NaN == NaN). Returns mismatch descriptions (empty = consistent).
"""
function compare_tidy_csv(new_path::AbstractString, base_path::AbstractString)::Vector{String}
    problems = String[]
    isfile(new_path) && isfile(base_path) ||
        return ["missing tidy CSV: $new_path / $base_path"]
    a = CSV.read(new_path, DataFrame; types = String)
    b = CSV.read(base_path, DataFrame; types = String)
    names(a) == names(b) || return ["$new_path: header mismatch: $(names(a)) vs $(names(b))"]
    if (nrow(a), ncol(a)) != (nrow(b), ncol(b))
        push!(problems, "$new_path: shape mismatch: $(size(a)) vs $(size(b))")
    end
    nr = min(nrow(a), nrow(b))
    for col in names(a), i in 1:nr
        sa = coalesce(a[i, col], "")
        sb = coalesce(b[i, col], "")
        va = tryparse(Float64, sa)
        vb = tryparse(Float64, sb)
        ok = if va !== nothing && vb !== nothing
            value_match(vb, va)
        else
            sa == sb
        end
        ok || push!(problems,
            "$new_path: cell ($i,$col) differs: baseline=\"$sb\" regenerated=\"$sa\"")
    end
    return problems
end

# ── Driver ──────────────────────────────────────────────────────────────

function reproduce_main(args::Vector{String} = ARGS)::Int
    mode = nothing
    for a in args
        if a == "--write-baseline" || a == "--check"
            mode === nothing || return reproduce_usage()
            mode = a
        else
            println(stderr, "unknown option \"$a\"")
            return reproduce_usage()
        end
    end
    mode === nothing && return reproduce_usage()
    root = default_root()

    design_d, sectors, recs = resolve_cells(root)
    n_pass = report_cells(recs)
    rows = all_metric_rows(recs, sectors)

    problems = String[]
    if mode == "--write-baseline"
        mkpath(joinpath(root, "baseline", "output"))
        write_metrics_csv(joinpath(root, "baseline", "v10_metrics.csv"), rows)
    else
        append!(problems,
            check_metrics_csv(joinpath(root, "baseline", "v10_metrics.csv"), rows))
    end

    rc = run_generator(root)
    rc == 0 || push!(problems,
        "documented generator failed (exit $rc): a cell of the headline set did not validate")

    if mode == "--write-baseline"
        copy_outputs(root)
        summary_roundtrip_check(root, recs)
    else
        for name in output_artifact_names()
            name == "variance_decomposition_sobol.csv" && continue  # not re-generated here
            append!(problems, compare_tidy_csv(
                joinpath(root, "output", name),
                joinpath(root, "baseline", "output", name)))
        end
    end

    n_pass == length(recs) || push!(problems,
        "$(length(recs) - n_pass)/$(length(recs)) cells failed re-validation")
    if isempty(problems)
        println(mode == "--write-baseline" ?
            "baseline written under baseline/; all $(length(recs)) cells validated" :
            "check PASSED: baseline consistent within $REPRO_TOL")
        return 0
    end
    println(stderr, "$(length(problems)) problem(s):")
    for p in problems
        println(stderr, "  " * p)
    end
    return 1
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    exit(reproduce_main())
end
