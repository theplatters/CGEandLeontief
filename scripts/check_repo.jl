# scripts/check_repo.jl — repository gate (Phase 3, ADR-0007)
#
# Usage:
#   julia --project=. scripts/check_repo.jl [--root DIR] [--quiet]
#
# Runs ten checks over the working tree and exits 1 on any violation:
#   1. board ......... docs/status.md is up to date and the registry reports 0 warnings
#   2. closures ...... registry/closures.toml is bidirectional with the kernel registry
#   3. closure-subtypes  every AbstractLabor/FinancingClosure subtype is a registered symbol
#   4. kernel-reachability  every src/**/*.jl is reachable from src/BeyondHulten.jl
#   5. runs .......... every runs/<id>/ has a schema-valid manifest, index, and scenario row
#   6. runs-visible .. executed/provisional/failed rows are manifest-backed or pre-manifest history
#   7. tracked-artifacts  no generated artifacts are tracked by git
#   8. kernel-copies . no living file includes a cbase2/src path
#   9. preregistration  every preregistration record matches its design file SHA-256
#  10. wip-limit ..... at most one scenario row is `running`
#
# The file can be `include`d by tests without running the CLI (guard at the
# bottom). Every check takes an explicit `root` so tests can point at
# fixtures. Public API: `check_repo(root)` (all violations),
# `check_board`, `check_closure_registry`, `check_closure_subtypes`,
# `check_kernel_reachability`, `check_runs`, `check_failed_visible`,
# `check_tracked_artifacts`, `check_kernel_copies`, `check_preregistration`,
# `check_wip_limit`, `check_main`.
#
# Headless-safe: only BeyondHulten plus stdlibs and the kernel's CSV/DataFrames deps.

using BeyondHulten
using TOML
using Dates
using SHA
using CSV
using DataFrames

include(joinpath(@__DIR__, "status.jl"))  # build_board, normalize (include-safe via its own PROGRAM_FILE guard)

"""One gate violation: which check failed and why."""
struct RepoViolation
    check::String
    message::String
end
Base.show(io::IO, v::RepoViolation) = print(io, "[", v.check, "] ", v.message)

const RUN_STATUSES = Set(["running", "executed", "provisional", "failed"])
const EXPECTED_INDEX_HEADER = "run_id,date,design,closures,status,gate_summary,headline_metrics,commit"

"""Pre-manifest historical rows (ADR-0004): the six cbase2-v3 verification runs
predate manifests and are never retrofitted with invented ones. Every other
executed/provisional/failed scenario row must be manifest-backed."""
const PRE_MANIFEST_ROWS = Set([
    "cbase2-v3-ALPHA-F1-mobile",
    "cbase2-v3-ALPHA-F2-mobile",
    "cbase2-v3-ALPHA-F3-mobile",
    "cbase2-v3-DELTA-analytic-F2",
    "cbase2-v3-DELTA-analytic-F3",
    "cbase2-v3-BETA-mobile",
])

"""Concrete closure-description types that are NOT taxonomy closures and hence
need no registry symbol: `ExogenousLaborClosure` (the legacy CES exogenous-labor
callback wrapper, ADR-0005) and `NoFinancing` (the unfinanced baseline reference;
an unfinanced experiment is explicitly inadmissible, ROADMAP §4.3)."""
const NON_CLOSURE_SUBTYPES = Set(["ExogenousLaborClosure", "NoFinancing"])

"""SHA-256 hex digest of a file's bytes."""
file_sha256(path::AbstractString) = bytes2hex(SHA.sha256(read(path)))

"""Read-only/frozen zone paths from registry/freeze.toml."""
function frozen_paths(root::AbstractString)::Set{String}
    out = Set{String}()
    fp = joinpath(root, "registry", "freeze.toml")
    if isfile(fp)
        fr = TOML.parsefile(fp)
        for section in ("frozen", "read_only")
            sec = get(fr, section, Dict())
            sec isa AbstractDict || continue
            for (_, v) in sec
                (v isa AbstractDict && haskey(v, "path")) || continue
                push!(out, strip(string(v["path"])))
            end
        end
    end
    return out
end

"""True when repo-relative `path` lies inside one of the frozen zones."""
function is_frozen(path::AbstractString, frozen::Set{String})::Bool
    p = replace(string(path), "\\" => "/")
    for z in frozen
        z == p && return true
        startswith(p, z * "/") && return true
    end
    return false
end

"""Tracked files (`git ls-files`), or `String[]` when git is unavailable."""
function tracked_files(root::AbstractString)::Vector{String}
    try
        out = read(`git -C $root ls-files`, String)
        return [l for l in split(out, '\n') if !isempty(strip(l))]
    catch
        return String[]
    end
end

"""All `.jl` files under a directory, sorted repo-relatively for stable output."""
function all_jl_files(dir::AbstractString)::Vector{String}
    out = String[]
    isdir(dir) || return out
    for (d, _, files) in walkdir(dir)
        for f in files
            endswith(f, ".jl") && push!(out, joinpath(d, f))
        end
    end
    sort!(out)
    return out
end

"""Closure entries as id => (axis, table) from registry/closures.toml."""
function registry_entries(root::AbstractString)::Dict{String,Tuple{String,Dict{String,Any}}}
    clos = TOML.parsefile(joinpath(root, "registry", "closures.toml"))
    out = Dict{String,Tuple{String,Dict{String,Any}}}()
    for axis in ("labor", "financing")
        sec = get(clos, axis, Dict())
        sec isa AbstractDict || continue
        for (key, v) in sec
            v isa AbstractDict || continue
            t = Dict{String,Any}(string(k) => val for (k, val) in v)
            id = haskey(t, "id") ? string(t["id"]) : string(key)
            out[id] = (axis, t)
        end
    end
    return out
end

"""All code symbols claimed by the closure registry."""
function registry_symbols(entries::Dict{String,Tuple{String,Dict{String,Any}}})::Set{String}
    syms = Set{String}()
    for (_, (_, t)) in entries
        v = get(t, "symbols", String[])
        v === nothing && continue
        for s in (v isa AbstractVector ? v : [v])
            st = strip(string(s))
            (isempty(st) || st == "TBD") || push!(syms, st)
        end
    end
    return syms
end

# ── Check 1: board up to date + zero registry warnings ───────────────────

"""Rebuild the status board, require zero warnings, and require docs/status.md to match."""
function check_board(root::AbstractString)::Vector{RepoViolation}
    out = RepoViolation[]
    board, nwarn = build_board(root)
    nwarn == 0 || push!(out, RepoViolation("board",
        "registry reports $nwarn warning(s); fix the references (see the Warnings section of the generated board)"))
    outpath = joinpath(root, "docs", "status.md")
    if !isfile(outpath)
        push!(out, RepoViolation("board", "docs/status.md is missing; run `julia --project=. scripts/status.jl` first"))
    elseif normalize(read(outpath, String)) != normalize(board)
        push!(out, RepoViolation("board", "docs/status.md is stale; run `julia --project=. scripts/status.jl` to regenerate"))
    end
    return out
end

# ── Check 2: closure registry is bidirectional ───────────────────────────

"""Every registry id resolves in the kernel and vice versa, axes agree,
non-idea closures have constructors, and ZETA has none."""
function check_closure_registry(root::AbstractString)::Vector{RepoViolation}
    out = RepoViolation[]
    entries = registry_entries(root)
    code_ids = Set(string.(closure_ids()))
    reg_ids = Set(keys(entries))
    reg_ids == code_ids || push!(out, RepoViolation("closures",
        "closure id mismatch: registry has $(sort(collect(reg_ids))), code has $(sort(collect(code_ids)))"))
    for id in sort(collect(reg_ids))
        axis, t = entries[id]
        status = strip(string(get(t, "status", "")))
        ok = try
            string(closure_axis(Symbol(id))) == axis
        catch
            false
        end
        ok || push!(out, RepoViolation("closures",
            "closure $id: axis mismatch or unknown to the kernel (registry says \"$axis\")"))
        ctor = try
            closure_constructor(Symbol(id))
        catch
            :__unknown__
        end
        if id == "ZETA"
            (ctor === nothing && status == "idea") || push!(out, RepoViolation("closures",
                "closure ZETA must have no constructor and status `idea`"))
        elseif status != "idea"
            (ctor !== nothing && ctor !== :__unknown__) || push!(out, RepoViolation("closures",
                "closure $id: status `$status` requires a kernel constructor"))
        end
    end
    return out
end

# ── Check 3: no unregistered labour-closure subtypes ─────────────────────

"""Every `struct X <: AbstractLaborClosure` / `AbstractFinancing` in src/ must
appear in some closure entry's symbols (or be listed in NON_CLOSURE_SUBTYPES)."""
function check_closure_subtypes(root::AbstractString)::Vector{RepoViolation}
    out = RepoViolation[]
    syms = registry_symbols(registry_entries(root))
    for jl in all_jl_files(joinpath(root, "src"))
        contents = try
            read(jl, String)
        catch
            continue
        end
        for m in eachmatch(r"struct\s+([A-Za-z_]\w*)\s*<:\s*Abstract(Labor|Financing)Closure", contents)
            t = string(m.captures[1])
            (t in syms || t in NON_CLOSURE_SUBTYPES) || push!(out, RepoViolation("closure-subtypes",
                "unregistered closure subtype `$t` in $(relpath(jl, root)): add it to a closure entry's symbols"))
        end
    end
    return out
end

# ── Check 4: no zombie kernel files ──────────────────────────────────────

"""Every src/**/*.jl must be reachable from src/BeyondHulten.jl via include()."""
function check_kernel_reachability(root::AbstractString)::Vector{RepoViolation}
    out = RepoViolation[]
    entry = joinpath(root, "src", "BeyondHulten.jl")
    if !isfile(entry)
        return [RepoViolation("kernel-reachability", "src/BeyondHulten.jl is missing")]
    end
    seen = Set{String}()
    stack = String[normpath(entry)]
    while !isempty(stack)
        f = pop!(stack)
        f in seen && continue
        push!(seen, f)
        contents = try
            read(f, String)
        catch
            push!(out, RepoViolation("kernel-reachability", "included file missing: $(relpath(f, root))"))
            continue
        end
        for m in eachmatch(r"include\s*\(\s*\"([^\"]+)\"", contents)
            push!(stack, normpath(joinpath(dirname(f), string(m.captures[1]))))
        end
    end
    for jl in all_jl_files(joinpath(root, "src"))
        normpath(jl) in seen || push!(out, RepoViolation("kernel-reachability",
            "zombie kernel file $(relpath(jl, root)): not reachable from src/BeyondHulten.jl (promote it, or archive it with a DE record)"))
    end
    return out
end

# ── Check 5: runs are manifest-backed ────────────────────────────────────

"""Load scenarios.csv rows keyed by run_id (all values as strings)."""
function scenario_rows(root::AbstractString)::Dict{String,Dict{String,String}}
    scenpath = joinpath(root, "registry", "scenarios.csv")
    isfile(scenpath) || return Dict{String,Dict{String,String}}()
    raw = DataFrame(CSV.File(scenpath))
    out = Dict{String,Dict{String,String}}()
    for r in eachrow(raw)
        d = Dict{String,String}(string(c) => string(coalesce(r[c], "")) for c in names(raw))
        out[d["run_id"]] = d
    end
    return out
end

"""Every runs/<id>/ directory has a schema-valid manifest; the index and the
scenario rows agree with the manifests (ADR-0004/ADR-0006)."""
function check_runs(root::AbstractString)::Vector{RepoViolation}
    out = RepoViolation[]
    runs_dir = joinpath(root, "runs")
    isdir(runs_dir) || return out
    manifests = Dict{String,Dict{String,Any}}()
    for entry in sort(readdir(runs_dir))
        full = joinpath(runs_dir, entry)
        isdir(full) || continue
        manpath = joinpath(full, "manifest.toml")
        if !isfile(manpath)
            push!(out, RepoViolation("runs", "runs/$entry/ has no manifest.toml"))
            continue
        end
        man = try
            TOML.parsefile(manpath)
        catch e
            push!(out, RepoViolation("runs", "runs/$entry/manifest.toml does not parse: $(sprint(showerror, e))"))
            continue
        end
        for k in ("schema_version", "run_id", "design", "status", "date", "actor")
            haskey(man, k) || push!(out, RepoViolation("runs", "runs/$entry/manifest.toml is missing key `$k`"))
        end
        for t in ("provenance", "scenario", "solver", "artifacts")
            (haskey(man, t) && man[t] isa AbstractDict) ||
                push!(out, RepoViolation("runs", "runs/$entry/manifest.toml is missing table `[$t]`"))
        end
        if haskey(man, "run_id") && string(man["run_id"]) != entry
            push!(out, RepoViolation("runs", "runs/$entry/manifest.toml run_id \"$(man["run_id"])\" does not match the directory"))
        end
        if haskey(man, "status") && !(string(man["status"]) in RUN_STATUSES)
            push!(out, RepoViolation("runs", "runs/$entry/manifest.toml has invalid status \"$(man["status"])\""))
        end
        # The run journal: run.jl appends progress lines to log.txt (ADR-0004 raw results/logs).
        isfile(joinpath(full, "log.txt")) ||
            push!(out, RepoViolation("runs", "runs/$entry/ has no log.txt (run journal)"))
        manifests[entry] = man
    end
    # runs/index.csv: exact header, exactly one row per manifest and vice versa.
    idxpath = joinpath(runs_dir, "index.csv")
    if !isfile(idxpath)
        push!(out, RepoViolation("runs", "runs/index.csv is missing"))
    else
        raw = read(idxpath, String)
        lines = split(raw, '\n')
        isempty(lines) || lines[1] == EXPECTED_INDEX_HEADER ||
            push!(out, RepoViolation("runs", "runs/index.csv has wrong header: \"$(lines[1])\""))
        idx_ids = Set{String}()
        if filesize(idxpath) > 0
            try
                df = DataFrame(CSV.File(idxpath; stringtype = String))
                hasproperty(df, :run_id) || push!(out, RepoViolation("runs", "runs/index.csv has no run_id column"))
                for r in eachrow(df)
                    push!(idx_ids, string(coalesce(r.run_id, "")))
                end
            catch e
                push!(out, RepoViolation("runs", "runs/index.csv does not parse: $(sprint(showerror, e))"))
            end
        end
        man_ids = Set(keys(manifests))
        for id in sort(collect(setdiff(man_ids, idx_ids)))
            push!(out, RepoViolation("runs", "runs/$id/ has a manifest but no runs/index.csv row"))
        end
        for id in sort(collect(setdiff(idx_ids, man_ids)))
            push!(out, RepoViolation("runs", "runs/index.csv row `$id` has no runs/$id/manifest.toml"))
        end
    end
    # Every manifest has a scenarios.csv row with matching labour/financing.
    rows = scenario_rows(root)
    for id in sort(collect(keys(manifests)))
        haskey(rows, id) || push!(out, RepoViolation("runs",
            "runs/$id/manifest.toml has no registry/scenarios.csv row `$id`"))
        continue_check = haskey(rows, id)
        continue_check || continue
        scen = get(manifests[id], "scenario", Dict{String,Any}())
        for (mkey, ckey) in (("labour", "labor"), ("financing", "financing"))
            haskey(scen, mkey) || continue
            string(scen[mkey]) == rows[id][ckey] || push!(out, RepoViolation("runs",
                "runs/$id/manifest.toml $mkey \"$(scen[mkey])\" does not match scenarios.csv $ckey \"$(rows[id][ckey])\""))
        end
    end
    return out
end

# ── Check 6: failed/provisional runs stay visible ────────────────────────

"""Executed/provisional/failed scenario rows are manifest-backed, unless they
are pre-manifest history (PRE_MANIFEST_ROWS, ADR-0004)."""
function check_failed_visible(root::AbstractString)::Vector{RepoViolation}
    out = RepoViolation[]
    for (id, row) in scenario_rows(root)
        row["status"] in ("executed", "provisional", "failed") || continue
        isfile(joinpath(root, "runs", id, "manifest.toml")) && continue
        id in PRE_MANIFEST_ROWS && continue
        push!(out, RepoViolation("runs-visible",
            "scenario `$id` is `$(row["status"])` but has no runs/$id/manifest.toml (failed/provisional runs stay visible; never drop them)"))
    end
    return out
end

# ── Check 7: no tracked generated artifacts ─────────────────────────────

"""No tracked, non-frozen file may be a generated artifact: Manifest.toml,
*.log, .ipynb_checkpoints, plots/*.png, output/ contents, or runs/ contents
other than runs/index.csv and runs/*/manifest.toml."""
function check_tracked_artifacts(root::AbstractString)::Vector{RepoViolation}
    out = RepoViolation[]
    frozen = frozen_paths(root)
    for f in tracked_files(root)
        is_frozen(f, frozen) && continue
        rel = "tracked artifact `$f`"
        if basename(f) == "Manifest.toml"
            push!(out, RepoViolation("tracked-artifacts", "$rel: Manifest.toml is generated; do not commit it"))
            continue
        end
        if endswith(f, ".log")
            push!(out, RepoViolation("tracked-artifacts", "$rel: *.log files are transient; do not commit them"))
            continue
        end
        if ".ipynb_checkpoints" in split(replace(f, "\\" => "/"), "/")
            push!(out, RepoViolation("tracked-artifacts", "$rel: notebook checkpoints are transient; do not commit them"))
            continue
        end
        if occursin(r"(^|/)plots/[^/]*\.png$", replace(f, "\\" => "/"))
            push!(out, RepoViolation("tracked-artifacts", "$rel: generated figures are not committed (regenerate them)"))
            continue
        end
        if f == "output" || startswith(f, "output/")
            push!(out, RepoViolation("tracked-artifacts", "$rel: output/ is generated diagnostic output; do not commit it"))
            continue
        end
        if f == "runs" || startswith(f, "runs/")
            f == "runs/index.csv" && continue
            match(r"^runs/[^/]+/manifest\.toml$", f) !== nothing && continue
            push!(out, RepoViolation("tracked-artifacts", "$rel: only runs/index.csv and runs/*/manifest.toml are committed (ADR-0004)"))
            continue
        end
    end
    return out
end

# ── Check 8: no living kernel copies ─────────────────────────────────────

"""No tracked, non-frozen file may include a cbase2/src path (ADR-0001)."""
function check_kernel_copies(root::AbstractString)::Vector{RepoViolation}
    out = RepoViolation[]
    frozen = frozen_paths(root)
    for f in tracked_files(root)
        is_frozen(f, frozen) && continue
        full = joinpath(root, f)
        (isfile(full) && filesize(full) < 5_000_000) || continue
        contents = try
            read(full, String)
        catch
            continue
        end
        occursin(r"include\s*\(.*cbase2/src", contents) || continue
        push!(out, RepoViolation("kernel-copies",
            "`$f` includes a cbase2/src path: the root src/ is the one kernel (ADR-0001); read cbase2, never wire it in"))
    end
    return out
end

# ── Check 9: preregistration integrity ───────────────────────────────────

"""Every preregistration record matches the SHA-256 of its design file."""
function check_preregistration(root::AbstractString)::Vector{RepoViolation}
    out = RepoViolation[]
    preg_path = joinpath(root, "registry", "preregistration.toml")
    isfile(preg_path) || return out
    preg = try
        TOML.parsefile(preg_path)
    catch e
        return [RepoViolation("preregistration", "registry/preregistration.toml does not parse: $(sprint(showerror, e))")]
    end
    designs = get(preg, "designs", Dict())
    for design in sort(collect(keys(designs)))
        rec = designs[design]
        recorded = rec isa AbstractDict ? strip(string(get(rec, "design_sha256", ""))) : ""
        dpath = joinpath(root, "experiments", "designs", design * ".toml")
        if !isfile(dpath)
            push!(out, RepoViolation("preregistration", "preregistered design `$design` is missing experiments/designs/$design.toml"))
            continue
        end
        current = file_sha256(dpath)
        recorded == current || push!(out, RepoViolation("preregistration",
            "design `$design` changed since preregistration (recorded $(recorded[1:min(8, end)]) vs current $(current[1:8])); re-run `experiments/run.jl --preregister $design`"))
    end
    return out
end

# ── Check 10: WIP limit ──────────────────────────────────────────────────

"""At most one scenario row may be `running` at a time (ADR-0007)."""
function check_wip_limit(root::AbstractString)::Vector{RepoViolation}
    out = RepoViolation[]
    running = sort([id for (id, row) in scenario_rows(root) if row["status"] == "running"])
    length(running) > 1 && push!(out, RepoViolation("wip-limit",
        "WIP limit exceeded: $(length(running)) scenarios `running` ($(join(running, ", "))); at most one at a time (ADR-0007)"))
    return out
end

# ── Driver ───────────────────────────────────────────────────────────────

const CHECKS = [
    ("board", check_board),
    ("closures", check_closure_registry),
    ("closure-subtypes", check_closure_subtypes),
    ("kernel-reachability", check_kernel_reachability),
    ("runs", check_runs),
    ("runs-visible", check_failed_visible),
    ("tracked-artifacts", check_tracked_artifacts),
    ("kernel-copies", check_kernel_copies),
    ("preregistration", check_preregistration),
    ("wip-limit", check_wip_limit),
]

"""Run all checks over `root`; returns the violation list (empty means clean)."""
function check_repo(root::AbstractString)::Vector{RepoViolation}
    out = RepoViolation[]
    for (name, fn) in CHECKS
        try
            append!(out, fn(root))
        catch e
            push!(out, RepoViolation(name, "check errored: $(typeof(e)): $(sprint(showerror, e))"))
        end
    end
    return out
end

"""Print usage and exit with code 1."""
function check_usage()::Nothing
    println(stderr, "usage:\n  julia --project=. scripts/check_repo.jl [--root DIR] [--quiet]")
    exit(1)
end

"""CLI entry point (only via the `PROGRAM_FILE` guard below)."""
function check_main(args::Vector{String} = ARGS)::Int
    root = rootdir()
    quiet = false
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--root"
            i += 1
            i > length(args) && check_usage()
            root = args[i]
        elseif a == "--quiet"
            quiet = true
        else
            check_usage()
        end
        i += 1
    end
    vs = check_repo(root)
    for (name, _) in CHECKS
        n = count(v -> v.check == name, vs)
        quiet || println(n == 0 ? "[ok] $name" : "[FAIL] $name ($n)")
    end
    if isempty(vs)
        quiet || println("check_repo: 0 violations")
    else
        println("check_repo: $(length(vs)) violation(s):")
        for v in sort(vs; by = v -> (v.check, v.message))
            println("- ", v)
        end
    end
    return isempty(vs) ? 0 : 1
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    exit(check_main())
end
