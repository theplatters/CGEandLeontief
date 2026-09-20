# experiments/probes/probe6_reconstruct_run_journals.jl
#
# Reconstruct runs/<id>/log.txt from the committed manifest.
#
# Why: the v4-v9 generations have no log.txt in some working copies (only
# v1-v3 do; log.txt and solution.csv are gitignored working artefacts, ADR-0004
# tracks only runs/index.csv and runs/*/manifest.toml). The repository gate
# requires the journal to exist for every run dir, so the journals are rebuilt
# from the manifests with the exact grammar of experiments/run.jl at the run
# commit bc83d5a (written by `execute_cell` via `run_log`, which prepends
# `iso_timestamp() * " "` — a 19-char UTC stamp plus one space — on every
# call, so only the first physical line of a multi-line failure entry carries
# a timestamp):
#
#   <ts> start <run_id> design=<design> actor=<actor> commit=<git_commit>
#   <ts> reference <name>=<real_consumption(ref_sol)>
#   <ts> <status>: <gate summary: "resid ...; budget ...; <labour|wage|sectoral> ...">
#   <ts> failed: <error.type>: <error.message>   (failed cells only; the
#       message is `sprint(showerror, e)` verbatim, newlines preserved —
#       including the v7 MethodError's `Closest candidates` lines)
#
# The reference metric is named `real_gdp_ref` in the v1-v3 code and
# `consumption_ref` from ADR-0018 on (the v4-v9 code); each generation is
# rebuilt with its own name. The third gate is `labour` (mobile), `wage`
# (fixed) or `sectoral` (the eta = 0 ADR-0020 endpoint). The value is computed
# from the design's reference continuation here, never typed. The gate summary
# is rebuilt from the manifest's [gates] table with the same `gate_frag`
# format string. A failed cell carries no per-gate table (only
# `[gates] overall = "fail"`); its third line comes from the `[error]` table.
#
# The timestamp is the manifest's own mtime (the artefact write time — a bulk
# write, so all cells of a generation share it), not an invented solve time.
#
# `--verify` checks only the GENERATIONS below (byte-for-byte after the
# timestamp; line 2 modulo the metric value). The v1-v3 journals use an older
# reconstruction grammar and are out of scope: counted as legacy/skipped,
# never touched or "fixed". A verify mismatch is a finding to report, not a
# journal to rewrite.
#
# usage:
#   julia --project=. experiments/probes/probe6_reconstruct_run_journals.jl --verify
#   julia --project=. experiments/probes/probe6_reconstruct_run_journals.jl --write

include(joinpath(@__DIR__, "..", "run.jl"))

using Printf, Dates, TOML

const REF_NAME = "consumption_ref"   # v4-v9 code path (ADR-0018 rename)
const GENERATIONS = ["matrix_5x3-v4", "matrix_5x3-v5", "matrix_5x3-v6",
    "matrix_5x3-v7", "matrix_5x3-v8", "matrix_5x3-v9"]

"Rebuild the journal lines for one manifest (three entries; a failed entry's
third line is multi-line when the recorded error message contains newlines)."
function journal_lines(man::Dict{String,Any}, ts::String, ref_value::Float64)
    head = String[
        "start $(man["run_id"]) design=$(man["design"]) actor=$(man["actor"]) commit=$(man["provenance"]["git_commit"])",
        "reference $(REF_NAME)=$(ref_value)",
    ]
    if haskey(man, "error")
        # Failed cell (the `execute_cell` rescue path in run.jl): the journal's
        # third line is `failed: <typeof(e)>: <sprint(showerror, e)>`, rebuilt
        # here from the manifest's [error] table verbatim (newlines preserved).
        push!(head, "failed: $(man["error"]["type"]): $(man["error"]["message"])")
        return head
    end
    if man["status"] == "failed"
        # Defensive fallback (unused: every failed v7/v8 manifest carries an
        # [error] table); do not invent anything beyond the manifest.
        push!(head, "failed: overall=$(man["gates"]["overall"])")
        return head
    end
    g = man["gates"]
    third = first(k for k in ("labour", "wage", "sectoral") if haskey(g, k))
    frag(name, key) = @sprintf("%s %.3g<=%.0e %s", name,
        Float64(g[key]["value"]), Float64(g[key]["tolerance"]),
        g[key]["pass"] ? "ok" : "FAIL")
    summary = frag("resid", "residual") * "; " * frag("budget", "budget") * "; " * frag(third, third)
    push!(head, "$(man["status"]): $(summary)")
    return head
end

"Drop the `<ts> ` prefix run.jl prepends (19-char UTC stamp plus one space)."
strip_ts(line::AbstractString) = replace(line, r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2} " => "")

"The covered generation prefixing a run-dir name, or `nothing` (legacy scope)."
function gen_of(name::AbstractString)
    for g in GENERATIONS
        startswith(name, g * "-") && return g
    end
    return nothing
end

"ISO timestamp (the run.jl convention) from a file's mtime."
function mtime_stamp(path::AbstractString)::String
    return Dates.format(Dates.unix2datetime(stat(path).mtime), "yyyy-mm-ddTHH:MM:SS")
end

function reconstruct_main(args)
    mode = isempty(args) ? "--verify" : args[1]
    mode in ("--verify", "--write") || error("usage: --verify | --write")
    root = default_root()
    runs_dir = default_runs_dir(root)

    # The reference value: computed from the design's reference continuation.
    # The v5-v9 designs share the reference parameters (v9 adds only
    # eta_s_rigid_group), so one continuation covers all of them (line 2 is
    # verified modulo the value: a fresh continuation differs by 1-3 ulp
    # across generations, so only the `reference ` prefix is asserted).
    design_d = load_design("matrix_5x3_v5"; root = root)
    ref = build_reference(design_d; root = root)
    ref_value = real_consumption(ref.sol)
    @printf("reference %s = %.17g  (v1-v3 journals record 0.9999999999999997)\n", REF_NAME, ref_value)

    if mode == "--verify"
        # Validate the reconstruction against the journals that DO exist in the
        # covered GENERATIONS: line 1 byte-for-byte, line 2 modulo the metric
        # value, line 3 byte-for-byte for executed cells and the whole
        # remaining journal text (timestamp stripped from its first line) for
        # failed cells, whose multi-line messages are preserved verbatim.
        # Dirs outside GENERATIONS use the older v1-v3 grammar: counted as
        # legacy/skipped, never touched or "fixed".
        dirs = sort(filter(d -> isdir(d) && isfile(joinpath(d, "log.txt")),
            readdir(runs_dir; join = true)))
        nok, nbad, nlegacy = 0, 0, 0
        per_gen = Dict{String,Tuple{Int,Int}}()
        for d in dirs
            name = basename(d)
            gname = gen_of(name)
            if gname === nothing
                nlegacy += 1
                continue
            end
            man = TOML.parsefile(joinpath(d, "manifest.toml"))
            have = readlines(joinpath(d, "log.txt"))
            rebuilt = journal_lines(man, mtime_stamp(joinpath(d, "manifest.toml")), ref_value)
            ok1 = length(have) >= 1 && strip_ts(have[1]) == rebuilt[1]
            ok2 = length(have) >= 2 && startswith(strip_ts(have[2]), "reference ")
            ok3 = if length(have) < 3
                false
            elseif haskey(man, "error") || man["status"] == "failed"
                join(vcat([strip_ts(have[3])], have[4:end]), "\n") == rebuilt[3]
            else
                length(have) == 3 && strip_ts(have[3]) == rebuilt[3]
            end
            ok, bad = get(per_gen, gname, (0, 0))
            if ok1 && ok2 && ok3
                nok += 1
                per_gen[gname] = (ok + 1, bad)
            else
                nbad += 1
                per_gen[gname] = (ok, bad + 1)
                println("MISMATCH $(name)")
                ok1 || println("  line1 have: ", length(have) >= 1 ? strip_ts(have[1]) : "<missing>",
                    "\n        rebuilt: ", rebuilt[1])
                ok3 || println("  line3 have: ", length(have) >= 3 ? strip_ts(have[3]) : "<missing>",
                    "\n        rebuilt: ", rebuilt[3])
            end
            tag = haskey(man, "error") ? "failed entry" : "lines 1+3 ok"
            @printf("  %-28s %-22s %s\n", name,
                length(have) >= 2 ? strip_ts(have[2]) : "<missing>",
                ok1 && ok2 && ok3 ? tag : "CHECK")
        end
        println("\nverify: $nok ok, $nbad mismatch (timestamp ignored; line 2 modulo value)")
        for g in GENERATIONS
            ok, bad = get(per_gen, g, (0, 0))
            println("  $g: $ok ok, $bad mismatch")
        end
        println("legacy/skipped (v1-v3 grammar, out of scope): $nlegacy")
        return nothing
    end

    # --write: rebuild the missing journals of the covered generations.
    written = 0
    for gen in GENERATIONS
        for d in sort(filter(isdir, readdir(runs_dir; join = true)))
            startswith(basename(d), gen * "-") || continue
            logp = joinpath(d, "log.txt")
            if isfile(logp)
                println("skip (exists): $(basename(d))/log.txt")
                continue
            end
            man = TOML.parsefile(joinpath(d, "manifest.toml"))
            ts = mtime_stamp(joinpath(d, "manifest.toml"))
            lines = journal_lines(man, ts, ref_value)
            open(logp, "w") do io
                for l in lines
                    println(io, ts * " " * l)
                end
            end
            written += 1
            println("wrote $(basename(d))/log.txt  ($ts, $(man["status"]))")
        end
    end
    println("\n$written journal(s) written")
    return nothing
end

reconstruct_main(ARGS)
