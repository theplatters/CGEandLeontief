# experiments/probes/probe6_reconstruct_run_journals.jl
#
# Reconstruct runs/<id>/log.txt from the committed manifest.
#
# Why: the v4 and v5 generations have no log.txt in the working copy (only
# v1-v3 do; log.txt and solution.csv are gitignored working artefacts, ADR-0004
# tracks only runs/index.csv and runs/*/manifest.toml). The repository gate
# requires the journal to exist for every run dir, so the journals are rebuilt
# from the manifests with the exact grammar of experiments/run.jl at the run's
# commit (three lines, written by `execute_cell` via `run_log`):
#
#   <ts> start <run_id> design=<design> actor=<actor> commit=<git_commit>
#   <ts> reference <name>=<real_consumption(ref_sol)>
#   <ts> <status>: <gate summary: "resid ...; budget ...; <labour|wage> ...">
#
# The reference metric is named `real_gdp_ref` in the v1-v3 code and
# `consumption_ref` from ADR-0018 on (the v4/v5 code); each generation is
# rebuilt with its own name. The value is computed from the design's reference
# continuation here, never typed. The gate summary is rebuilt from the
# manifest's [gates] table with the same `gate_frag` format string.
#
# The timestamp is the manifest's own mtime (the artefact write time — a bulk
# write, so all cells of a generation share it), not an invented solve time.
#
# usage:
#   julia --project=. experiments/probes/probe6_reconstruct_run_journals.jl --verify
#   julia --project=. experiments/probes/probe6_reconstruct_run_journals.jl --write

include(joinpath(@__DIR__, "..", "run.jl"))

using Printf, Dates, TOML

const REF_NAME = "consumption_ref"   # v4/v5 code path (ADR-0018 rename)
const GENERATIONS = ["matrix_5x3-v4", "matrix_5x3-v5"]

"Rebuild the three journal lines for one manifest."
function journal_lines(man::Dict{String,Any}, ts::String, ref_value::Float64)
    g = man["gates"]
    third = first(k for k in ("labour", "wage") if haskey(g, k))
    frag(name, key) = @sprintf("%s %.3g<=%.0e %s", name,
        Float64(g[key]["value"]), Float64(g[key]["tolerance"]),
        g[key]["pass"] ? "ok" : "FAIL")
    summary = frag("resid", "residual") * "; " * frag("budget", "budget") * "; " * frag(third, third)
    return String[
        "start $(man["run_id"]) design=$(man["design"]) actor=$(man["actor"]) commit=$(man["provenance"]["git_commit"])",
        "reference $(REF_NAME)=$(ref_value)",
        "$(man["status"]): $(summary)",
    ]
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
    # The two generations share the design parameters, so one continuation
    # covers both (verified against the value recorded in the v1-v3 journals).
    design_d = load_design("matrix_5x3_v5"; root = root)
    ref = build_reference(design_d; root = root)
    ref_value = real_consumption(ref.sol)
    @printf("reference %s = %.17g  (v1-v3 journals record 0.9999999999999997)\n", REF_NAME, ref_value)

    if mode == "--verify"
        # Validate the reconstruction against journals that DO exist (v1-v3):
        # lines 1 and 3 must match byte-for-byte, line 2 modulo the metric name.
        dirs = sort(filter(d -> isdir(d) && isfile(joinpath(d, "log.txt")),
            readdir(runs_dir; join = true)))
        nok, nbad = 0, 0
        for d in dirs
            man = TOML.parsefile(joinpath(d, "manifest.toml"))
            have = readlines(joinpath(d, "log.txt"))
            rebuilt = journal_lines(man, mtime_stamp(joinpath(d, "manifest.toml")), ref_value)
            ok1 = have[1][19:end] == rebuilt[1]                 # drop the timestamp
            ok3 = have[3][19:end] == rebuilt[3]
            name_ok = startswith(have[2][19:end], "reference ")
            if ok1 && ok3 && name_ok
                nok += 1
            else
                nbad += 1
                println("MISMATCH $(basename(d))")
                ok1 || println("  line1 have: ", have[1][19:end], "\n        rebuilt: ", rebuilt[1])
                ok3 || println("  line3 have: ", have[3][19:end], "\n        rebuilt: ", rebuilt[3])
            end
            @printf("  %-28s %-22s %s\n", basename(d), have[2][19:end], ok1 && ok3 ? "lines 1+3 ok" : "CHECK")
        end
        println("\nverify: $nok ok, $nbad mismatch (lines 1 and 3, timestamp ignored)")
        return nothing
    end

    # --write: rebuild the missing journals of the two generations.
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
