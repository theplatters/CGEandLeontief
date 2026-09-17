# scripts/status.jl — generate docs/status.md from registry/ (ADR-0003)
#
# Usage:
#   julia --project=. scripts/status.jl          # write docs/status.md
#   julia --project=. scripts/status.jl --check  # exit 1 if the board is stale

using TOML
using Dates
using Printf
using CSV
using DataFrames

const CLOSURE_STATES = ("idea", "spec", "implemented", "tested", "validated", "rejected", "superseded")
const SCENARIO_STATES = ("planned", "running", "executed", "provisional", "failed", "superseded", "cancelled")
const PARAM_COLS = ("eta", "eta_s", "theta", "epsilon", "sigma")

"""Repository root derived from this file's location."""
rootdir() = normpath(joinpath(@__DIR__, ".."))

"""Run a git command, returning stripped stdout or `nothing` on failure."""
function gitout(root::AbstractString, args::AbstractString...)::Union{String,Nothing}
    try
        return strip(read(`git -C $root $args`, String))
    catch
        return nothing
    end
end

"""Escape a value for a Markdown table cell."""
function esc(x)::String
    s = x === missing ? "" : string(x)
    s = replace(s, r"\s+" => " ")
    return replace(strip(s), "|" => "\\|")
end

"""Fetch an optional key from a TOML table with a default."""
getk(d, k, default) = isa(d, AbstractDict) && haskey(d, k) ? d[k] : default

"""Coerce a TOML value to a vector of strings (missing -> empty)."""
function asvec(x)::Vector{String}
    x === nothing && return String[]
    x === missing && return String[]
    isa(x, AbstractVector) || return [string(x)]
    return [string(v) for v in x]
end

"""True when a registry string cell carries no information."""
isblank(s)::Bool = (t = strip(string(s === missing ? "" : s)); isempty(t) || t == "TBD")

"""Shorten a hash for display."""
short_hash(s)::String = (t = strip(string(s)); length(t) > 7 ? t[1:7] : (isempty(t) ? "—" : t))

"""Title of a DE/ADR file: first line starting with `# `, or the file name."""
function doctitle(path::String)::String
    try
        for line in eachline(path)
            startswith(line, "# ") && return strip(line[3:end])
        end
    catch
    end
    return basename(path)
end

"""List `(stem, filename, title)` triples for `prefix`-`*.md` files, sorted by file."""
function docindex(dir::String, prefix::String)::Vector{Tuple{String,String,String}}
    out = Tuple{String,String,String}[]
    isdir(dir) || return out
    for f in sort(readdir(dir))
        (startswith(f, prefix) && endswith(f, ".md")) || continue
        push!(out, (splitext(f)[1], f, doctitle(joinpath(dir, f))))
    end
    return out
end

"""Check whether `id` has a matching `<id>-*.md` file in `dir`."""
function hasdoc(dir::String, id::AbstractString)::Bool
    isdir(dir) || return false
    return any(f -> (startswith(f, id * "-") || startswith(f, id * "_")) && endswith(f, ".md"), readdir(dir))
end

"""Load one closure axis (`labor`/`financing`) as id-sorted dicts."""
function load_axis(clos::Dict, axis::String)::Vector{Dict{String,Any}}
    raw = get(clos, axis, Dict())
    rows = Dict{String,Any}[]
    for (key, v) in raw
        d = Dict{String,Any}(string(k) => val for (k, val) in v)
        haskey(d, "id") || (d["id"] = string(key))
        push!(rows, d)
    end
    sort!(rows, by = r -> string(getk(r, "id", "")))
    return rows
end

"""Validate closure path lists and symbols; append to `warns`. Returns file cache."""
function check_closure(row::Dict{String,Any}, root::String, warns::Vector{String})::Dict{String,String}
    id = string(getk(row, "id", "?"))
    status = string(getk(row, "status", ""))
    isempty(status) || status in CLOSURE_STATES ||
        push!(warns, "closure $id: unknown status \"$status\"")
    cache = Dict{String,String}()
    for field in ("files", "references", "tests", "evidence")
        for p in asvec(getk(row, field, String[]))
            isblank(p) && continue
            full = joinpath(root, strip(p))
            (isfile(full) || isdir(full)) ||
                push!(warns, "closure $id: $field not found: $p")
        end
    end
    # Read each existing listed file once for the symbol grep.
    for p in asvec(getk(row, "files", String[]))
        isblank(p) && continue
        full = joinpath(root, strip(p))
        isfile(full) || continue
        haskey(cache, full) || (cache[full] = try read(full, String) catch; "" end)
    end
    for sym in asvec(getk(row, "symbols", String[]))
        isblank(sym) && continue
        any(contents -> occursin(sym, contents), values(cache)) ||
            push!(warns, "closure $id: symbol \"$sym\" not found in listed files")
    end
    for de in asvec(getk(row, "dead_ends", String[]))
        isblank(de) && continue
        hasdoc(joinpath(root, "docs", "dead-ends"), strip(de)) ||
            push!(warns, "closure $id: dead_ends record $de has no file in docs/dead-ends/")
    end
    for a in asvec(getk(row, "adrs", String[]))
        isblank(a) && continue
        hasdoc(joinpath(root, "docs", "decisions"), strip(a)) ||
            push!(warns, "closure $id: ADR record $a has no file in docs/decisions/")
    end
    return cache
end

"""Format the `parameters` cell of a scenario row."""
function params(row)::String
    parts = String[]
    for c in PARAM_COLS
        v = hasproperty(row, Symbol(c)) ? row[Symbol(c)] : ""
        s = v === missing ? "" : strip(string(v))
        (isempty(s) || s == "TBD") || push!(parts, "$c=$s")
    end
    return isempty(parts) ? "TBD" : join(parts, ", ")
end

"""Print one Markdown table row from raw cell values."""
function trow(buf::IOBuffer, cells::Vector)::Nothing
    println(buf, "| ", join(esc.(cells), " | "), " |")
    return nothing
end

"""Render a closures table (labour or financing)."""
function closure_table(rows::Vector{Dict{String,Any}})::String
    buf = IOBuffer()
    trow(buf, ["ID", "Status", "Formulation", "Implementation", "Tests", "Dead ends", "Open gates"])
    println(buf, "| --- | --- | --- | --- | --- | --- | --- |")
    for r in rows
        trow(buf, [getk(r, "id", ""), getk(r, "status", ""), getk(r, "formulation", ""),
            join(asvec(getk(r, "files", String[])), ", "), join(asvec(getk(r, "tests", String[])), ", "),
            join(asvec(getk(r, "dead_ends", String[])), ", "), join(asvec(getk(r, "open_gates", String[])), "; ")])
    end
    return String(take!(buf))
end

"""Render the labour x financing scenario matrix."""
function matrix_table(df::DataFrame, labids::Vector{String}, finids::Vector{String})::String
    cell = Dict{Tuple{String,String},Vector{String}}()
    for r in eachrow(df)
        l = string(coalesce(r.labor, ""))
        f = string(coalesce(r.financing, ""))
        push!(get!(cell, (l, f), String[]), string(coalesce(r.status, "")))
    end
    buf = IOBuffer()
    trow(buf, vcat(["Labour \\ Financing"], finids))
    println(buf, "| --- |", join(fill(" ---", length(finids)), " |"), " |")
    for l in labids
        trow(buf, vcat([l], [haskey(cell, (l, f)) ? join(cell[(l, f)], ", ") : "—" for f in finids]))
    end
    return String(take!(buf))
end

"""Render one design's run table."""
function design_table(df::DataFrame)::String
    buf = IOBuffer()
    trow(buf, ["run_id", "labour", "financing", "status", "parameters", "evidence", "notes"])
    println(buf, "| --- | --- | --- | --- | --- | --- | --- |")
    for r in eachrow(sort(df, :run_id))
        trow(buf, [coalesce(r.run_id, ""), coalesce(r.labor, ""), coalesce(r.financing, ""),
            coalesce(r.status, ""), params(r), coalesce(r.evidence, ""), coalesce(r.notes, "")])
    end
    return String(take!(buf))
end

"""Build the full board text (including the volatile header block)."""
function build_board(root::String)::Tuple{String,Int}
    warns = String[]
    # --- inputs (missing registry file is a hard error) ---
    clospath = joinpath(root, "registry", "closures.toml")
    scenpath = joinpath(root, "registry", "scenarios.csv")
    freezepath = joinpath(root, "registry", "freeze.toml")
    isfile(clospath) || (println(stderr, "missing registry file: registry/closures.toml"); exit(1))
    isfile(scenpath) || (println(stderr, "missing registry file: registry/scenarios.csv"); exit(1))
    isfile(freezepath) || (println(stderr, "missing registry file: registry/freeze.toml"); exit(1))

    clos = TOML.parsefile(clospath)
    lab = load_axis(clos, "labor")
    fin = load_axis(clos, "financing")
    labids = [string(getk(r, "id", "")) for r in lab]
    finids = [string(getk(r, "id", "")) for r in fin]
    for r in vcat(lab, fin)
        check_closure(r, root, warns)
    end

    df = DataFrame(CSV.File(scenpath; stringtype=String, silencewarnings=true))
    for c in ("run_id", "design", "status", "labor", "financing", "evidence", "notes")
        hasproperty(df, Symbol(c)) || (df[!, Symbol(c)] = fill("", nrow(df)))
    end
    sort!(df, :run_id)
    seen = Set{String}()
    for r in eachrow(df)
        rid = string(coalesce(r.run_id, ""))
        if rid in seen
            push!(warns, "scenario $rid: duplicate run_id")
        else
            push!(seen, rid)
        end
        l, f, st = string(coalesce(r.labor, "")), string(coalesce(r.financing, "")),
            string(coalesce(r.status, ""))
        (isempty(l) || l in labids) || push!(warns, "scenario $rid: unknown labor id \"$l\"")
        (isempty(f) || f in finids) || push!(warns, "scenario $rid: unknown financing id \"$f\"")
        st in SCENARIO_STATES || push!(warns, "scenario $rid: invalid status \"$st\"")
        ev = string(coalesce(r.evidence, ""))
        for p in split(ev, ';')
            q = strip(p)
            (isempty(q) || q == "TBD") && continue
            ispath(joinpath(root, q)) || push!(warns, "scenario $rid: evidence not found: $q")
        end
    end

    fr = TOML.parsefile(freezepath)
    zones = Tuple{String,Dict{String,Any}}[]  # (key, table)
    for section in ("frozen", "read_only")
        sec = get(fr, section, Dict())
        for (k, v) in sec
            push!(zones, (string(k), Dict{String,Any}(string(kk) => vv for (kk, vv) in v)))
        end
    end
    sort!(zones, by = z -> z[1])
    fallback = strip(string(get(fr, "frozen_at_commit", "")))
    for (key, z) in zones
        p = strip(string(getk(z, "path", "")))
        isempty(p) && (push!(warns, "freeze zone $key: no path recorded"); continue)
        full = joinpath(root, p)
        ispath(full) || (push!(warns, "freeze zone $key ($p): path missing on disk"); continue)
        baseline = strip(string(getk(z, "recorded_commit", "")))
        isempty(baseline) && (baseline = fallback)
        if isempty(baseline)
            push!(warns, "freeze zone $key ($p): no baseline commit found")
            continue
        end
        recorded = strip(string(getk(z, "tree_hash", "")))
        if !isempty(recorded)
            actual = gitout(root, "rev-parse", "$baseline:$p")
            if actual !== nothing && actual != recorded
                push!(warns, "freeze zone $key ($p): tree_hash differs from $(short_hash(baseline)):$p " *
                             "($(short_hash(recorded)) vs $(short_hash(actual))) — modified since freeze")
            end
        end
        changed = String[]
        for out in (gitout(root, "diff", "--name-only", baseline, "--", p),
                    gitout(root, "ls-files", "--others", "--exclude-standard", "--", p))
            out === nothing && continue
            append!(changed, [l for l in split(out, '\n') if !isempty(strip(l))])
        end
        # The FROZEN.md marker is freeze metadata, not a zone modification.
        filter!(l -> strip(l) != string(p) * "/FROZEN.md", changed)
        isempty(changed) || push!(warns, "freeze zone $key ($p): dirty working tree — modified since freeze")
    end

    de_docs = docindex(joinpath(root, "docs", "dead-ends"), "DE-")
    adr_docs = docindex(joinpath(root, "docs", "decisions"), "ADR-")

    # --- volatile header (date, branch, HEAD, clean/dirty) ---
    branch = something(gitout(root, "rev-parse", "--abbrev-ref", "HEAD"), "unknown")
    head = something(gitout(root, "rev-parse", "--short", "HEAD"), "unknown")
    porch = gitout(root, "status", "--porcelain")
    clean = porch === nothing ? "unknown" : (isempty(porch) ? "clean" : "dirty")
    volatile = "<!-- volatile:start -->\nGenerated $(Dates.today()) by `scripts/status.jl` " *
               "· branch `$branch` · HEAD `$head` ($clean)\n<!-- volatile:end -->"

    # --- summary counts ---
    countby(rows, st) = count(r -> string(getk(r, "status", "")) == st, rows)
    scount(st) = count(r -> string(coalesce(r.status, "")) == st, eachrow(df))
    designs = sort(unique([string(coalesce(r.design, "")) for r in eachrow(df) if !isblank(coalesce(r.design, ""))]))

    buf = IOBuffer()
    println(buf, "# Status Board\n")
    println(buf, volatile, "\n")
    println(buf, "> Single source of truth: `registry/` (closures.toml, scenarios.csv, freeze.toml).")
    println(buf, "> Generated file — do not edit by hand (ADR-0003). See `registry/README.md`.\n")
    println(buf, "## Summary\n")
    println(buf, "| Axis | idea | spec | implemented | tested | validated | total |")
    println(buf, "| --- | ---: | ---: | ---: | ---: | ---: | ---: |")
    @printf(buf, "| Labour closures | %d | %d | %d | %d | %d | %d |\n",
        countby(lab, "idea"), countby(lab, "spec"), countby(lab, "implemented"),
        countby(lab, "tested"), countby(lab, "validated"), length(lab))
    @printf(buf, "| Financing closures | %d | %d | %d | %d | %d | %d |\n\n",
        countby(fin, "idea"), countby(fin, "spec"), countby(fin, "implemented"),
        countby(fin, "tested"), countby(fin, "validated"), length(fin))
    scenbits = join(["$s: $(scount(s))" for s in SCENARIO_STATES if scount(s) > 0], ", ")
    println(buf, "Scenarios: **$(nrow(df))** rows — $scenbits. " *
                 "Designs: $(join(["`$d`" for d in designs], ", ")).\n")
    println(buf, "## Labour closures\n")
    println(buf, closure_table(lab))
    println(buf, "## Financing closures\n")
    println(buf, closure_table(fin))
    println(buf, "## Scenario matrix\n")
    println(buf, matrix_table(df, labids, finids))
    println(buf, "### Runs by design\n")
    for d in designs
        println(buf, "#### `$d`\n")
        println(buf, design_table(filter(r -> string(coalesce(r.design, "")) == d, df)))
    end
    println(buf, "## Freeze board\n")
    trow(buf, ["Zone", "Kind", "Status", "Frozen at", "Files", "Tree hash", "Reason"])
    println(buf, "| --- | --- | --- | --- | --- | --- | --- |")
    for (key, z) in zones
        frozenat = getk(z, "recorded_commit_short", getk(z, "recorded_commit", getk(z, "last_touch_commit", "")))
        frozenat = isempty(strip(string(frozenat))) ? "—" : short_hash(frozenat)
        files = getk(z, "tracked_files", "")
        files = (files === "" || files === nothing) ? "—" : string(files)
        trow(buf, [key, getk(z, "kind", ""), getk(z, "status", ""), frozenat, files,
            short_hash(getk(z, "tree_hash", "")), getk(z, "reason", "")])
        extras = ["open: $item" for item in asvec(getk(z, "open_items", String[]))]
        for extra in ("successor", "notes")
            v = strip(string(getk(z, extra, "")))
            isempty(v) || push!(extras, "$extra: $v")
        end
        for e in extras
            println(buf, "  - ", esc(e))
        end
    end
    println(buf)
    for (heading, docs, prefix) in (("Dead ends", de_docs, "dead-ends"), ("Decisions", adr_docs, "decisions"))
        println(buf, "## $heading\n")
        if isempty(docs)
            println(buf, "None recorded.\n")
        else
            for (id, f, title) in docs
                println(buf, "- [$title]($prefix/$f)")
            end
            println(buf)
        end
    end
    println(buf, "## Warnings\n")
    if isempty(warns)
        println(buf, "None.\n")
    else
        for w in sort(warns)
            println(buf, "- ", esc(w))
        end
        println(buf)
    end
    return String(take!(buf)), length(warns)
end

"""Strip the volatile block and trailing whitespace per line for `--check`."""
function normalize(s::String)::String
    t = replace(s, r"<!-- volatile:start -->.*?<!-- volatile:end -->"s => "<!-- volatile -->")
    return join([rstrip(l) for l in split(t, '\n')], "\n")
end

"""Entry point: write the board, or compare it in `--check` mode."""
function main()::Nothing
    root = rootdir()
    check = "--check" in ARGS
    board, nwarn = build_board(root)
    outpath = joinpath(root, "docs", "status.md")
    if check
        isfile(outpath) || (println(stderr, "docs/status.md is missing; run scripts/status.jl first."); exit(1))
        if normalize(read(outpath, String)) != normalize(board)
            println(stderr, "docs/status.md is stale; run `julia --project=. scripts/status.jl` to regenerate.")
            exit(1)
        end
        println("docs/status.md is up to date.")
        return nothing
    end
    write(outpath, board)
    clos = TOML.parsefile(joinpath(root, "registry", "closures.toml"))
    nlab = length(load_axis(clos, "labor"))
    nfin = length(load_axis(clos, "financing"))
    nscen = nrow(DataFrame(CSV.File(joinpath(root, "registry", "scenarios.csv"); stringtype=String, silencewarnings=true)))
    println("docs/status.md written ($nlab labour, $nfin financing, $nscen scenarios, $nwarn warnings)")
    return nothing
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
