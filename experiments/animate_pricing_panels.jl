# experiments/animate_pricing_panels.jl — animate the manuscript panels across
# the markup slope kappa (ROADMAP Phase 7, ADR-0029/ADR-0030; see
# `experiments/README.md`, "Animating the pricing panels").
#
# Usage (from the repo root; needs a Makie backend, see below):
#
#   julia --project=. experiments/animate_pricing_panels.jl
#       [--design pricing_markup_v1] [--cells a,b,c] [--kappa-max 0.5]
#       [--steps 40] [--fps 8] [--formats mp4,gif] [--outdir output/animations]
#
# One animation per selected cell (default: the design's cells at the largest
# kappa, i.e. one per financing column). Each frame is the manuscript panel
# layout of `experiments/plot_matrix.jl` — grouped bars over the programme
# sectors (left) and the quantity/price scatter over all sectors (right) — at
# one kappa on the grid `kappa_max * (0:steps)/steps`, with axis limits FIXED
# across frames so the motion is the model's, not the autoscaler's.
#
# Validation contract (ADR-0025: figures are not runs). Frames are re-solves by
# continuation in kappa (each warm-started from the previous frame) from the
# NoPricing twin; they are never written to `runs/` or the registry. Every
# grid kappa that coincides with an executed run of the same labour, financing
# and s_pi in the design is HARD-VALIDATED against that run's `solution.csv`
# (max |dp|, |dq| < 1e-6, the design's uniqueness tolerance) and its manifest
# status; a mismatch aborts that animation. The validated run ids are printed
# on the frames. Frames between them are interpolating re-solves, not
# citable numbers.
#
# Backend: CairoMakie (file output, no GL context). It is not a project
# dependency (the committed Manifest pins the kernel generation); add it on
# demand, e.g. in a separate environment stacked onto this one:
#   julia --project=<env> -e 'using Pkg; Pkg.add("CairoMakie")'
#   JULIA_LOAD_PATH="@:<env>:@stdlib" julia --project=. experiments/animate_pricing_panels.jl
# Writes only under --outdir (default `output/animations/`, gitignored).

using OpenCGE
using CSV
using DataFrames
using Printf

isdefined(Main, :plot_sector_labels) || include(joinpath(@__DIR__, "plot_matrix.jl"))

"""Parse the CLI into a NamedTuple; throws `ArgumentError` on bad input."""
function parse_anim_args(args::Vector{String})::NamedTuple
    o = Dict{String,String}("design" => "pricing_markup_v1", "kappa-max" => "0.5",
        "steps" => "40", "fps" => "8", "formats" => "mp4,gif",
        "outdir" => joinpath(default_root(), "output", "animations"))
    i = 1
    while i <= length(args)
        key = startswith(args[i], "--") ? args[i][3:end] : ""
        key in ("design", "cells", "kappa-max", "steps", "fps", "formats", "outdir") ||
            throw(ArgumentError("unknown option $(args[i])"))
        i + 1 <= length(args) || throw(ArgumentError("missing value for $(args[i])"))
        o[key] = args[i+1]
        i += 2
    end
    steps = parse(Int, o["steps"])
    steps >= 1 || throw(ArgumentError("--steps must be >= 1"))
    formats = String.(split(o["formats"], ","))
    all(f -> f in ("mp4", "gif"), formats) || throw(ArgumentError("--formats takes mp4 and/or gif"))
    return (design = o["design"], cells = haskey(o, "cells") ? String.(split(o["cells"], ",")) : nothing,
        kappa_max = parse(Float64, o["kappa-max"]), steps = steps, fps = parse(Int, o["fps"]),
        formats = formats, outdir = o["outdir"])
end

"""One frame's data: relative sectoral vectors and the aggregate readouts."""
struct AnimFrame
    kappa::Float64
    cell::MatrixCellData
    defl::Float64
    employment::Float64
    consumption::Float64
    profit_share::Float64
    mu_min::Float64
    mu_max::Float64
end

"""
Solve the frames of one cell along the kappa grid by continuation and
validate them against the design's executed runs. Returns
`(frames, anchors)` with `anchors` the validated `kappa => run_id` pairs.
"""
function solve_frames(cell::Dict{String,Any}, design::AbstractString, design_d::Dict{String,Any},
        ref, ψ::Vector{Float64}, g::Vector{Float64}, labels::Vector{String},
        kappas::Vector{Float64}; root::AbstractString = default_root())
    pr0, _ = cell_pricing(design_d, cell)
    pr0 isa MarkupPricing || throw(ArgumentError("the animation needs a markup (P2) cell"))
    s_pi = pr0.s_pi
    # Anchors: executed runs of the same labour/financing/s_pi at a grid kappa.
    anchors = Dict{Float64,String}()
    for id in cell_order(design; root = root)
        c = design_cell(design_d, id)
        (c["labor"] == cell["labor"] && c["financing"] == cell["financing"]) || continue
        pc, _ = cell_pricing(design_d, c)
        (pc isa MarkupPricing && pc.s_pi == s_pi) || continue
        k = findfirst(x -> isapprox(x, pc.κ; atol = 1e-12), kappas)
        k === nothing || (anchors[kappas[k]] = id)
    end
    twin = build_cell_model(cell, design_d, ref.data, ψ, g; pricing = NoPricing())
    x = OpenCGE._solution_vector(twin, solve(twin; init = ref.init_warm))
    frames = AnimFrame[]
    for κ in kappas
        m = build_cell_model(cell, design_d, ref.data, ψ, g; pricing = MarkupPricing(κ; s_pi = s_pi))
        sol = solve(m; init = x)
        x = OpenCGE._solution_vector(m, sol)
        if haskey(anchors, κ)
            id = anchors[κ]
            man = TOML.parsefile(joinpath(root, "runs", id, "manifest.toml"))
            man["status"] == "executed" || error("anchor $id is not executed")
            p, q = read_stored_solution(joinpath(root, "runs", id, "solution.csv"))
            dev = max(maximum(abs.(sol.prices_raw .- p)), maximum(abs.(sol.quantities .- q)))
            dev < 1e-6 || error("frame kappa = $κ does not reproduce $id (max dev $dev)")
            println("  validated kappa = $κ against $id (max dev $(@sprintf("%.1e", dev)))")
        end
        mcd = matrix_cell_data(cell_id_label(cell, κ), string(cell["labor"]), string(cell["financing"]),
            Float64(cell["eta"]), Float64(get(cell, "eta_s", 0.0)), "re-solve",
            ref.data, sol, ref.sol; labels = labels)
        μ = pricing_markup(m.pricing, sol.quantities, ref.data.λ)
        L = sum(sectoral_labor_demand(sol.prices_raw, sol.quantities, sol.wages_raw[1], m))
        push!(frames, AnimFrame(κ, mcd, gdp_deflator(sol, ref.sol), L,
            real_consumption(sol) / real_consumption(ref.sol), profit_share(m, sol),
            minimum(μ), maximum(μ)))
    end
    return frames, anchors
end

cell_id_label(cell, κ) = @sprintf("%s-%s-kappa%.4f", cell["labor"], cell["financing"], κ)

"""Padded fixed limits over every frame's values."""
function fixed_limits(vals; pad = 0.08, min_span = 0.05)
    f = filter(isfinite, vals)
    lo, hi = minimum(f), maximum(f)
    span = max(hi - lo, min_span)
    return (lo - pad * span, hi + pad * span)
end

"""Record one animation (every requested format) for a solved frame sequence."""
function record_cell(M, frames::Vector{AnimFrame}, anchors::Dict{Float64,String},
        programme::Vector{Float64}, baseline_y::Vector{Float64}, title_stub::String,
        paths_base::String; fps::Int, formats::Vector{String})
    M = M.Makie   # the backend re-exports Makie; utilities like wong_colors live there
    n = length(frames[1].cell.quantities)
    idx = [j for j in sortperm(programme; rev = true) if programme[j] > 0]
    labels = frames[1].cell.labels
    colors = M.wong_colors()
    pct(v) = 100 .* (v .- 1)
    # Left panel geometry: four series over the programme sectors.
    ns, width = 4, 0.7 / 4
    xs = [k + (si - 2.5) * width for si in 1:ns for k in eachindex(idx)]
    bar_colors = [colors[si] for si in 1:ns for _ in eachindex(idx)]
    prog_pct = [100 * programme[j] / baseline_y[j] for j in idx]
    heights(fr) = vcat(prog_pct, pct(fr.cell.quantities[idx]), pct(fr.cell.consumption[idx]),
        pct(fr.cell.prices[idx]))
    left_lim = fixed_limits(vcat(0.0, (heights(fr) for fr in frames)...))
    # Right panel: fixed limits over all frames and sectors.
    allq = vcat((pct(fr.cell.quantities) for fr in frames)...)
    allp = vcat((pct(fr.cell.prices) for fr in frames)...)
    xlim, ylim = fixed_limits(vcat(0.0, allq)), fixed_limits(vcat(0.0, allp))
    inprog = Set(idx)
    others = [j for j in 1:n if !(j in inprog)]
    # Labelled movers, fixed by their motion at the last frame.
    last = frames[end].cell
    move(j) = max(abs(pct(last.quantities[j])), abs(pct(last.prices[j])))
    tagged = vcat(sort(idx; by = move, rev = true)[1:min(3, length(idx))],
        sort(others; by = move, rev = true)[1:min(3, length(others))])

    k = M.Observable(1)
    fr = M.lift(i -> frames[i], k)
    fig = M.Figure(size = (1900, 980))
    M.Label(fig[1, 1:2], M.lift(f -> @sprintf("%s — κ = %.3f", title_stub, f.kappa), fr), fontsize = 22)
    M.Label(fig[2, 1:2], M.lift(f -> @sprintf(
        "ΔP (GDP deflator) %+.2f%%    ΔL %+.2f%%    ΔW (consumption) %+.2f%%    Π/Y %.2f%%    μ ∈ [%.3f, %.3f]",
        100 * (f.defl - 1), 100 * (f.employment - 1), 100 * (f.consumption - 1),
        100 * f.profit_share, f.mu_min, f.mu_max), fr), fontsize = 16)
    axl = M.Axis(fig[3, 1]; xlabel = "Sector", ylabel = "%", ytickformat = "{:.2f}%",
        xticks = (1:length(idx), labels[idx]), xticklabelrotation = π / 4,
        xticklabelsize = 11, yticklabelsize = 11)
    M.hlines!(axl, [0.0]; color = :gray, linestyle = :dash)
    M.barplot!(axl, xs, M.lift(heights, fr); width = width, color = bar_colors)
    M.ylims!(axl, left_lim...)
    M.axislegend(axl, [M.PolyElement(polycolor = colors[si]) for si in 1:ns],
        ["Programme demand (% of baseline gross output)", "Change in output",
            "Change in consumption", "Price deviation from numeraire"];
        position = :rt, labelsize = 11)
    axr = M.Axis(fig[3, 2]; xlabel = "Change in quantities (%)", ylabel = "Change in prices (%)",
        xtickformat = "{:.2f}%", ytickformat = "{:.2f}%", xticklabelsize = 11, yticklabelsize = 11)
    M.hlines!(axr, [0.0]; color = :gray, linestyle = :dash)
    M.vlines!(axr, [0.0]; color = :gray, linestyle = :dash)
    # Ghost of the kappa = 0 frame: where the sectors start.
    f0 = frames[1].cell
    M.scatter!(axr, pct(f0.quantities), pct(f0.prices); color = (:gray, 0.25), markersize = 8)
    pts(v, js) = [M.Point2f(pct(v.cell.quantities[j]), pct(v.cell.prices[j])) for j in js]
    M.scatter!(axr, M.lift(f -> pts(f, others), fr); color = colors[1], markersize = 10)
    M.scatter!(axr, M.lift(f -> pts(f, idx), fr); color = colors[2], markersize = 16)
    # Programme movers sit apart; the non-programme movers crowd the origin,
    # so their labels get separate vertical lanes.
    nprog = min(3, length(idx))
    lanes = [(8, -16), (8, -32), (8, 14)]
    offsets = [k <= nprog ? (6, 4) : lanes[mod1(k - nprog, 3)] for k in eachindex(tagged)]
    M.text!(axr, M.lift(f -> pts(f, tagged), fr); text = labels[tagged], fontsize = 11,
        align = (:left, :bottom), offset = offsets)
    M.xlims!(axr, xlim...)
    M.ylims!(axr, ylim...)
    M.axislegend(axr, [M.MarkerElement(marker = :circle, color = colors[2]),
            M.MarkerElement(marker = :circle, color = colors[1]),
            M.MarkerElement(marker = :circle, color = (:gray, 0.4))],
        ["Programme sectors", "Other sectors", "κ = 0 (start)"]; position = :rb, labelsize = 12)
    anchor_txt = join([@sprintf("κ = %.2f ↔ %s", κ, anchors[κ]) for κ in sort(collect(keys(anchors)))], ";  ")
    M.Label(fig[4, 1:2], M.lift(f -> haskey(anchors, f.kappa) ?
        "This frame reproduces the executed run $(anchors[f.kappa]) (validated)." :
        "Interpolating re-solve between executed runs — not a run, not citable (ADR-0025).", fr),
        fontsize = 12, color = :gray25)
    M.Label(fig[5, 1:2], _anim_note(anchor_txt), fontsize = 10)
    M.rowsize!(fig.layout, 4, M.Fixed(18))
    M.rowsize!(fig.layout, 5, M.Fixed(18))
    # Hold the first and last frames, and briefly each validated anchor.
    seq = Int[]
    for (i, f) in enumerate(frames)
        reps = i == 1 ? fps : i == length(frames) ? 2 * fps :
               haskey(anchors, f.kappa) ? fps ÷ 2 : 1
        append!(seq, fill(i, max(reps, 1)))
    end
    out = String[]
    for ext in formats
        path = paths_base * "." * ext
        M.record(fig, path, seq; framerate = fps) do i
            k[] = i
        end
        push!(out, path)
    end
    return out
end

_anim_note(anchors) = "Shock effect relative to the design reference (no-programme baseline). " *
    "Frames: re-solves by continuation in κ, fixed axes. Validated against: " * anchors * "."

"""Driver: returns the written file paths."""
function anim_main(args::Vector{String} = ARGS; root::AbstractString = default_root())
    o = parse_anim_args(args)
    M = try
        Base.require(Main, :CairoMakie)
    catch e
        println(stderr, "CairoMakie is not available: add it in a stacked environment " *
            "(see the header of experiments/animate_pricing_panels.jl).")
        rethrow()
    end
    ok, reason = preregistration_status(o.design; root = root)
    ok || error("refusing to animate: $reason")
    design_d = load_design(o.design; root = root)
    order = cell_order(o.design; root = root)
    ids = if o.cells === nothing
        # One cell per labour x financing: the one with the largest kappa.
        best = Dict{Tuple{String,String,Float64},Tuple{Float64,String}}()
        for id in order
            c = design_cell(design_d, id)
            pc, _ = cell_pricing(design_d, c)
            pc isa MarkupPricing || continue
            key = (string(c["labor"]), string(c["financing"]), pc.s_pi)
            (!haskey(best, key) || pc.κ > best[key][1]) && (best[key] = (pc.κ, id))
        end
        sort([v[2] for v in values(best)])
    else
        o.cells
    end
    for id in ids
        id in order || throw(ArgumentError("unknown cell \"$id\" for design \"$(o.design)\""))
    end
    ref = build_reference(design_d; root = root)
    n = length(ref.data.factor_share)
    ψ, g = programme_vectors(design_d, ref.data, collect(1:n); root = root)
    labels = plot_sector_labels(design_d, n; root = root)
    kappas = [o.kappa_max * s / o.steps for s in 0:o.steps]
    mkpath(o.outdir)
    written = String[]
    for id in ids
        cell = design_cell(design_d, id)
        pc, _ = cell_pricing(design_d, cell)
        println("$id: solving $(length(kappas)) frames")
        frames, anchors = solve_frames(cell, o.design, design_d, ref, ψ, g, labels, kappas; root = root)
        stub = @sprintf("Markup door (P2, %s) — %s-%s",
            pc.s_pi == 1 ? "rent retained" : @sprintf("rent distributed, s_π = %.4g", pc.s_pi),
            cell["labor"], cell["financing"])
        base = joinpath(o.outdir, @sprintf("%s_%s_%s%s_kappa", o.design, cell["labor"],
            cell["financing"], pc.s_pi == 1 ? "" : @sprintf("_spi%.4g", pc.s_pi)))
        paths = Base.invokelatest(record_cell, M, frames, anchors, g, ref.sol.quantities, stub, base;
            fps = o.fps, formats = o.formats)
        foreach(p -> println("  wrote $p"), paths)
        append!(written, paths)
    end
    return written
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    anim_main()
end
