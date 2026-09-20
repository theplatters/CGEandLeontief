# experiments/plot_matrix.jl — comparative figures for the executed 5x3 matrix
# (see `experiments/README.md`, "Plotting the matrix").
#
# Usage (from the repo root, with `julia --project=.`):
#
#   julia --project=. experiments/plot_matrix.jl --design <design> [--cells a,b,c]
#       [--outdir DIR] [--data-dir DIR] [--no-figures] [--quiet]
#
# Validation contract: every selected cell is RE-SOLVED with the harness
# (`solve_cell` + `evaluate_gates` from `run.jl`) and HARD-VALIDATED against
# its recorded run artifacts (`runs/<id>/manifest.toml` + `solution.csv`)
# via `BeyondHulten.validate_cell`. Figures and tidy CSVs cover only the
# cells that validate; invalid cells are excluded from the dataset but
# listed in the report. The driver never creates run dirs and never modifies
# `runs/`, `registry/`, or any other tracked state — it only reads them
# (writing just the tidy CSVs under `<data-dir>/` and the figures under
# `<outdir>/`, both gitignored areas). The design must be preregistered with
# a matching SHA-256, tying the figures to the pinned design.
#
# Headless-safe with `--no-figures`: only `BeyondHulten`, stdlibs, and the
# kernel's CSV/DataFrames/TOML deps are loaded; GLMakie is touched only on
# the figure path.
#
# The file can be `include`d without running `plot_main()` (guard at the
# bottom). Entry point: `plot_main` (not `main`, which `run.jl` defines).

using BeyondHulten
using TOML
using CSV
using DataFrames
using Printf

# run.jl is include-safe (PROGRAM_FILE guard); skip the re-include when it is
# already loaded (as `tests/test_run_manifest.jl` does).
isdefined(Main, :run_design) || include(joinpath(@__DIR__, "run.jl"))

"""Print usage to stderr and return exit code 1."""
function plot_usage()::Int
    println(stderr, """usage:
  julia --project=. experiments/plot_matrix.jl --design <design> [--cells a,b,c] [--outdir DIR] [--data-dir DIR] [--no-figures] [--quiet]""")
    return 1
end

"""
Parse the driver CLI args into a NamedTuple (`design`, `cells`, `outdir`,
`data_dir`, `no_figures`, `quiet`). `--design` is required (checked by
`plot_main`); `--cells` is a comma-separated run-id list. Throws an
`ArgumentError` on a missing value or an unknown option.
"""
function parse_plot_args(args::Vector{String})::NamedTuple
    design::Union{String,Nothing} = nothing
    cells::Union{Vector{String},Nothing} = nothing
    outdir::Union{String,Nothing} = nothing
    data_dir::Union{String,Nothing} = nothing
    no_figures = false
    quiet = false
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--design" || a == "--cells" || a == "--outdir" || a == "--data-dir"
            i + 1 > length(args) && throw(ArgumentError("missing value for $a"))
            key = a[3:end]
            val = args[i+1]
            if key == "design"
                design = val
            elseif key == "cells"
                cells = unique(String.(strip.(split(val, ","))))
            elseif key == "outdir"
                outdir = val
            else
                data_dir = val
            end
            i += 2
        elseif a == "--no-figures"
            no_figures = true
            i += 1
        elseif a == "--quiet"
            quiet = true
            i += 1
        else
            throw(ArgumentError("unknown option \"$a\""))
        end
    end
    return (; design = design, cells = cells, outdir = outdir,
        data_dir = data_dir, no_figures = no_figures, quiet = quiet)
end

"""
`true` when the GLMakie package is installed (findable), without loading it.
The figure path `@eval using GLMakie` only after this check passes.
"""
function glmakie_available()::Bool
    return Base.find_package("GLMakie") !== nothing
end

"""
Short sector labels for the plots: the `N` names in columns `3:(2+N)` of the
design programme's `source` header row (`N` sectors). Falls back to
`"1".."N"` when the programme has no `source` key, the file is missing, or
the header is too short. Never hardcodes the impulses path beyond the
design's own `source` key.
"""
function plot_sector_labels(design_d::Dict{String,Any}, n::Integer;
        root::AbstractString = default_root())::Vector{String}
    try
        prog = design_d["programme"]
        haskey(prog, "source") || return string.(1:n)
        path = joinpath(root, prog["source"])
        isfile(path) || return string.(1:n)
        cols = names(DataFrame(CSV.File(path; limit = 1)))
        length(cols) >= 2 + n || return string.(1:n)
        return String.(cols[3:(2+n)])
    catch
        return string.(1:n)
    end
end

"""
Read a recorded `solution.csv` (`sector,price,quantity`; `stringtype =
String`) into `(prices, quantities)` float vectors. Returns `(nothing,
nothing)` when the file is absent (e.g. a failed cell, whose manifest
carries no solution); a present-but-unparseable file throws.
"""
function read_stored_solution(path::AbstractString)::Tuple
    isfile(path) || return (nothing, nothing)
    df = CSV.read(path, DataFrame; stringtype = String)
    return (Float64.(df.price), Float64.(df.quantity))
end

"""
Maximum over the `metric:*` entries of a `validate_cell` deltas dict;
`NaN` when no metric delta is finite.
"""
function plot_max_metric_delta(deltas::AbstractDict)::Float64
    best = NaN
    for (k, v) in deltas
        startswith(String(k), "metric:") || continue
        v isa Real && isfinite(v) || continue
        best = isnan(best) ? Float64(v) : max(best, Float64(v))
    end
    return best
end

"""
Driver entry point: preregistration gate, reference continuation, per-cell
re-solve + hard validation against the recorded artifacts, validation
report, tidy CSV export, and (unless `--no-figures`) the GLMakie figures.
Returns the process exit code (0 only when every selected cell validated
and figures were written or `--no-figures` was passed).
"""
function plot_main(args::Vector{String} = ARGS;
        root::AbstractString = default_root())::Int
    opts = try
        parse_plot_args(args)
    catch e
        println(stderr, sprint(showerror, e))
        return plot_usage()
    end
    if opts.design === nothing
        println(stderr, "missing required --design")
        return plot_usage()
    end
    design = opts.design
    quiet = opts.quiet

    # 1. The figures are tied to the pinned design: abort (before solving)
    # when the design file is not preregistered with a matching SHA-256.
    ok, reason = preregistration_status(design; root = root)
    if !ok
        println(stderr, "refusing to plot: $reason")
        return 1
    end

    # 2. Design, file order, and the wanted subset (every --cells id must be
    # a cell of this design; the default is the design's matrix cells).
    design_d = load_design(design; root = root)
    order = cell_order(design; root = root)
    wanted = if opts.cells === nothing
        matrix_cell_ids(order, design)
    else
        for id in opts.cells
            if !(id in order)
                println(stderr, "unknown cell \"$id\" for design \"$design\"")
                return 1
            end
        end
        String.(opts.cells)
    end
    if isempty(wanted)
        println(stderr, "no matrix cells selected for design \"$design\"")
        return 1
    end
    wanted_set = Set(wanted)
    outdir = opts.outdir === nothing ? joinpath(root, "plots") : opts.outdir
    data_dir = opts.data_dir === nothing ? joinpath(root, "output") : opts.data_dir

    # 3. Reference continuation and programme incidence, exactly as
    # `run_design` in `run.jl` computes them. Typed conversion (see the note
    # in `build_reference`): broadcasting over an empty `drops = []` yields
    # `Vector{Any}` once BeyondHulten is loaded, so `Vector{Int}(...)` is
    # inference-independent.
    quiet || println("building reference continuation for $design ...")
    ref = build_reference(design_d; root = root)
    drops = Vector{Int}(design_d["data"]["drops"])
    n_full = length(ref.data.factor_share) + length(drops)
    kept = sort(setdiff(1:n_full, drops))
    ψ, g = programme_vectors(design_d, ref.data, kept; root = root)

    # 4. Short sector labels from the design's own programme source.
    n = length(ref.data.factor_share)
    labels = plot_sector_labels(design_d, n; root = root)
    baseline = matrix_baseline(ref.data, ref.sol)

    # 5-6. Per cell, in file (canonical) order: re-solve, re-evaluate the
    # gates, and hard-validate against the recorded artifacts. A
    # solve/gate/read failure is recorded as an invalid cell and never
    # aborts the batch.
    cells = MatrixCellData[]
    validation = Dict{String,String}()
    man_commits = Dict{String,Any}()
    report = Tuple{String,String,Dict{String,Float64},String,String,Int}[]
    for id in order
        id in wanted_set || continue
        status = "n/a"
        deltas = Dict{String,Float64}()
        verdict = "FAIL"
        reason = ""
        n_skipped = 0
        try
            quiet || println("re-solving $id ...")
            cell = design_cell(design_d, id)
            sol = solve_cell(cell, design_d, ref.data, ψ, g, ref.init_warm)
            ev = evaluate_gates(cell, design_d, sol.model, sol, ref.sol)
            ev.gates["overall"] == "pass" ||
                throw(ErrorException("re-solved gates report \"$(ev.gates["overall"])\""))
            mandir = joinpath(root, "runs", id)
            manpath = joinpath(mandir, "manifest.toml")
            isfile(manpath) ||
                throw(ErrorException("recorded manifest not found: runs/$id/manifest.toml"))
            manifest = TOML.parsefile(manpath)
            status = string(get(manifest, "status", "?"))
            man_commits[id] = string(get(get(manifest, "provenance",
                Dict{String,Any}()), "git_commit", "unknown"))
            stored_prices, stored_quantities =
                read_stored_solution(joinpath(mandir, "solution.csv"))
            vc = validate_cell(manifest, ev.metrics, sol.prices_raw,
                sol.quantities, stored_prices, stored_quantities)
            deltas = vc.deltas
            n_skipped = length(vc.skipped)
            if vc.pass
                push!(cells, matrix_cell_data(id, string(cell["labor"]),
                    string(cell["financing"]), Float64(cell["eta"]),
                    Float64(get(cell, "eta_s", 0.0)), status,
                    ref.data, sol, ref.sol; metrics = ev.metrics,
                    diagnostics = ev.diagnostics, labels = labels))
                verdict = "ok"
            else
                reason = join(vc.failures, "; ")
            end
        catch e
            reason = sprint(showerror, e)
        end
        validation[id] = verdict == "ok" ?
            (n_skipped > 0 ? "ok ($n_skipped metrics skipped)" : "ok") :
            (isempty(reason) ? "failed" : reason)
        push!(report, (id, status, deltas, verdict, reason, n_skipped))
    end

    # 7. Assemble the dataset: the shared baseline, the valid cells, the
    # per-run validation map, and the provenance block.
    prov = Dict{String,Any}(
        "design" => design,
        "design_sha256" => design_sha256(design; root = root),
        "git_commit" => git_commit(root),
        "julia_version" => string(VERSION),
        "plotted_at" => iso_timestamp(),
        "preregistration" => "ok",
        "run_commits" => man_commits,
    )
    ds = matrix_dataset(design, baseline, cells; provenance = prov,
        validation = validation)

    # 8. Validation report: one line per selected cell (never silenced by
    # --quiet), with the deltas from `validate_cell`.
    println("validation report for $design:")
    println("run_id | status | max|Δp| | max|Δq| | max|Δmetric| | verdict | skipped")
    for (id, status, deltas, verdict, reason, n_skipped) in report
        dp = get(deltas, "price", NaN)
        dq = get(deltas, "quantity", NaN)
        dm = plot_max_metric_delta(deltas)
        tail = verdict == "ok" ? "ok" : "FAIL ($reason)"
        println("$id | $status | $(@sprintf("%.3g", dp)) | " *
            "$(@sprintf("%.3g", dq)) | $(@sprintf("%.3g", dm)) | $tail | skipped=$n_skipped")
    end

    # 9. Tidy data for ad-hoc exploration (valid cells only).
    mkpath(data_dir)
    sum_path = joinpath(data_dir, design * "_summary.csv")
    sec_path = joinpath(data_dir, design * "_sectoral.csv")
    CSV.write(sum_path, matrix_summary_frame(ds))
    CSV.write(sec_path, matrix_sectoral_frame(ds))
    println("wrote $sum_path")
    println("wrote $sec_path")

    # 10. Figures (skipped headless-safe with --no-figures, which never
    # loads GLMakie).
    if !opts.no_figures
        if !glmakie_available()
            println(stderr, "`BeyondHulten.save_matrix_figures` requires " *
                "the GLMakie extension, which is not loaded. " *
                "Run `using GLMakie` together with `using BeyondHulten` " *
                "(install it first with `import Pkg; Pkg.add(\"GLMakie\")` " *
                "if necessary).")
            return 1
        end
        try
            @eval using GLMakie
        catch e
            println(stderr, "could not load GLMakie: $(sprint(showerror, e))")
            return 1
        end
        mkpath(outdir)
        # `invokelatest`: the extension method is born when `using GLMakie`
        # runs above — after `plot_main` was already compiled — so a plain
        # call would dispatch in the caller's (older) world age and hit the
        # headless stub. `invokelatest` dispatches in the latest world.
        written = try
            Base.invokelatest(save_matrix_figures, ds;
                outdir = outdir, prefix = design)
        catch e
            println(stderr,
                "save_matrix_figures failed: $(sprint(showerror, e))")
            return 1
        end
        paths = written isa AbstractVector ? written :
            (written === nothing ? String[] : [written])
        for p in paths
            println("wrote $p")
        end
    end

    # 11. Exit code: 0 only when every selected cell validated (figures were
    # written above, or --no-figures skipped them).
    n_ok = count(v -> startswith(v, "ok"), values(validation))
    n_tot = length(wanted)
    println("$n_ok/$n_tot cells validated")
    return n_ok == n_tot ? 0 : 1
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    exit(plot_main())
end
