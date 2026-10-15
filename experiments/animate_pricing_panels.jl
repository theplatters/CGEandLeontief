# experiments/animate_pricing_panels.jl — animate the manuscript panels across
# the pricing-door parameter: the markup slope kappa (P2) or the capacity
# elasticity delta (P1) (ROADMAP Phase 7, ADR-0029/ADR-0030; see
# `experiments/README.md`, "Animating the pricing panels").
#
# Usage (from the repo root; needs a Makie backend, see below):
#
#   julia --project=. experiments/animate_pricing_panels.jl
#       [--design pricing_markup_v1] [--cells a,b,c] [--kappa-max 0.5]
#       [--delta-max 0.5] [--steps 40] [--fps 8] [--formats mp4,gif]
#       [--outdir output/animations]
#
# The swept parameter follows the design's pricing kind: kappa on the grid
# `kappa_max * (0:steps)/steps` for a markup design, delta on the symmetric
# grid `delta_max * (2s - steps)/steps`, s = 0..steps (steps even), for a
# capacity design. One animation per selected cell (default: one per labour x
# financing [x s_pi]). Each frame is the manuscript panel layout of
# `experiments/plot_matrix.jl` — grouped bars over the programme sectors
# (left) and the quantity/price scatter over all sectors (right) — at one grid
# value, with axis limits FIXED across frames so the motion is the model's,
# not the autoscaler's; a grey ghost marks the door-closed (parameter 0)
# positions.
#
# Validation contract (ADR-0025: figures are not runs). Frames are re-solves by
# continuation in the parameter (each warm-started from its neighbour nearer
# 0) from the NoPricing twin; they are never written to `runs/` or the
# registry. Every grid value that coincides with an executed run of the same
# labour, financing (and s_pi) in the design is HARD-VALIDATED against that
# run's `solution.csv`
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
    o = Dict{String,String}("design" => "pricing_markup_v1", "kappa-max" => "0.5", "delta-max" => "0.5",
        "steps" => "40", "fps" => "8", "formats" => "mp4,gif",
        "outdir" => joinpath(default_root(), "output", "animations"))
    i = 1
    while i <= length(args)
        key = startswith(args[i], "--") ? args[i][3:end] : ""
        key in ("design", "cells", "kappa-max", "delta-max", "steps", "fps", "formats", "outdir") ||
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
        kappa_max = parse(Float64, o["kappa-max"]), delta_max = parse(Float64, o["delta-max"]), steps = steps, fps = parse(Int, o["fps"]),
        formats = formats, outdir = o["outdir"])
end

"""One frame's data: relative sectoral vectors and the aggregate readouts.
`band` is the range of the door's sectoral wedge: the markup μ (P2) or the
technology factor A_eff/A = (y/λ)^(-δ) (P1). `cons_h` / `cons_g` split the
change in final consumption into its household and government parts (the
government part includes the programme where it is an additive purchase,
F2/F3; under F1 the programme is a household-demand tilt), each in % of the
sector's baseline total final consumption c_h0 + G0."""
struct AnimFrame
    param::Float64
    cell::MatrixCellData
    cons_h::Vector{Float64}
    cons_g::Vector{Float64}
    defl::Float64
    employment::Float64
    consumption::Float64
    profit_share::Float64
    band_min::Float64
    band_max::Float64
end

"""The swept door: its symbol, its value in a pricing, and a rebuild at a value."""
door_param(pr::MarkupPricing) = pr.κ
door_param(pr::CapacityPricing) = pr.δ
door_at(pr::MarkupPricing, v) = MarkupPricing(v; s_pi = pr.s_pi)
door_at(pr::CapacityPricing, v) = CapacityPricing(v)
door_symbol(::MarkupPricing) = "κ"
door_symbol(::CapacityPricing) = "δ"
door_name(::MarkupPricing) = "kappa"
door_name(::CapacityPricing) = "delta"
door_s_pi(pr::MarkupPricing) = pr.s_pi
door_s_pi(::CapacityPricing) = 1.0
door_band(pr::MarkupPricing, y, λ) = pricing_markup(pr, y, λ)
door_band(pr::CapacityPricing, y, λ) = (y ./ λ) .^ (-pr.δ)
door_band_label(::MarkupPricing) = "μ"
door_band_label(::CapacityPricing) = "A_eff/A"

"""
Solve the frames of one cell along the parameter grid `vals` (ascending,
containing 0) by continuation outward from 0 and validate them against the
design's executed runs. Returns `(frames, anchors)`, frames in grid order,
with `anchors` the validated `value => run_id` pairs.
"""
function solve_frames(cell::Dict{String,Any}, design::AbstractString, design_d::Dict{String,Any},
        ref, ψ::Vector{Float64}, g::Vector{Float64}, labels::Vector{String},
        vals::Vector{Float64}; root::AbstractString = default_root())
    pr0, _ = cell_pricing(design_d, cell)
    pr0 isa Union{MarkupPricing,CapacityPricing} ||
        throw(ArgumentError("the animation needs a markup (P2) or capacity (P1) cell"))
    name = door_name(pr0)
    # Anchors: executed runs of the same door, labour, financing (and s_pi) at a grid value.
    anchors = Dict{Float64,String}()
    for id in cell_order(design; root = root)
        c = design_cell(design_d, id)
        (c["labor"] == cell["labor"] && c["financing"] == cell["financing"]) || continue
        pc, _ = cell_pricing(design_d, c)
        (typeof(pc) == typeof(pr0) && door_s_pi(pc) == door_s_pi(pr0)) || continue
        k = findfirst(x -> isapprox(x, door_param(pc); atol = 1e-12), vals)
        k === nothing || (anchors[vals[k]] = id)
    end
    i0 = findfirst(iszero, vals)
    i0 === nothing && throw(ArgumentError("the parameter grid must contain 0"))
    twin = build_cell_model(cell, design_d, ref.data, ψ, g; pricing = NoPricing())
    x0 = OpenCGE._solution_vector(twin, solve(twin; init = ref.init_warm))
    frames = Vector{AnimFrame}(undef, length(vals))
    # Continue upward from 0, then downward from 0 (the δ = 0 point equals the
    # NoPricing twin), each frame warm-started from its neighbour.
    for branch in (i0:length(vals), i0-1:-1:1)
        x = x0
        for i in branch
            v = vals[i]
            m = build_cell_model(cell, design_d, ref.data, ψ, g; pricing = door_at(pr0, v))
            sol = solve(m; init = x)
            x = OpenCGE._solution_vector(m, sol)
            if haskey(anchors, v)
                id = anchors[v]
                man = TOML.parsefile(joinpath(root, "runs", id, "manifest.toml"))
                man["status"] == "executed" || error("anchor $id is not executed")
                p, q = read_stored_solution(joinpath(root, "runs", id, "solution.csv"))
                dev = max(maximum(abs.(sol.prices_raw .- p)), maximum(abs.(sol.quantities .- q)))
                dev < 1e-6 || error("frame $name = $v does not reproduce $id (max dev $dev)")
                println("  validated $name = $v against $id (max dev $(@sprintf("%.1e", dev)))")
            end
            mcd = matrix_cell_data(cell_id_label(cell, name, v), string(cell["labor"]),
                string(cell["financing"]), Float64(cell["eta"]), Float64(get(cell, "eta_s", 0.0)),
                "re-solve", ref.data, sol, ref.sol; labels = labels)
            band = door_band(m.pricing, sol.quantities, ref.data.λ)
            L = sum(sectoral_labor_demand(sol.prices_raw, sol.quantities, sol.wages_raw[1], m))
            ch, cg = consumption_split(m, sol, ref)
            frames[i] = AnimFrame(v, mcd, ch, cg, gdp_deflator(sol, ref.sol), L,
                real_consumption(sol) / real_consumption(ref.sol), profit_share(m, sol),
                minimum(band), maximum(band))
        end
    end
    return frames, anchors
end

"""
Change in final consumption by type, in % of baseline total final
consumption (gross household `c_h0` plus government `G0`, the reference
levels): household `100 (c_h - c_h0) / (c_h0 + G0)` and government
`100 (G + g_add - G0) / (c_h0 + G0)`, with `g_add` the financing's additive
programme purchase (zero under F1). Sectors with no baseline consumption are NaN.
"""
function consumption_split(m, sol, ref)
    n = length(ref.data.factor_share)
    ch0, G0 = Float64.(ref.sol.consumption), Float64.(ref.data.gov_demand)
    base = ch0 .+ G0
    gov = m.data.gov_demand .+ OpenCGE.additive_demand(m.financing, n)
    share(d) = [base[i] > 0 ? 100 * d[i] / base[i] : NaN for i in 1:n]
    return share(sol.consumption .- ch0), share(gov .- G0)
end

cell_id_label(cell, name, v) = @sprintf("%s-%s-%s%.4f", cell["labor"], cell["financing"], name, v)

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
        paths_base::String; fps::Int, formats::Vector{String}, sym::String = "κ",
        band_label::String = "μ", legend_pos::Symbol = :rb)
    M = M.Makie   # the backend re-exports Makie; utilities like wong_colors live there
    n = length(frames[1].cell.quantities)
    idx = [j for j in sortperm(programme; rev = true) if programme[j] > 0]
    labels = frames[1].cell.labels
    colors = M.wong_colors()
    pct(v) = 100 .* (v .- 1)
    # Left panel geometry: four groups over the programme sectors; the third
    # (final consumption) is a stacked bar, household solid and government
    # hatched, each segment stacked from 0 on its own side of the axis.
    ns, width = 4, 0.7 / 4
    slot(si) = [k + (si - 2.5) * width for k in eachindex(idx)]
    xs = vcat(slot(1), slot(2), slot(4))
    bar_colors = vcat(fill(colors[1], length(idx)), fill(colors[2], length(idx)),
        fill(colors[4], length(idx)))
    prog_pct = [100 * programme[j] / baseline_y[j] for j in idx]
    heights(fr) = vcat(prog_pct, pct(fr.cell.quantities[idx]), pct(fr.cell.prices[idx]))
    # Government segment sits on top of the household one when both share a sign.
    gov_from(fr) = [sign(fr.cons_g[j]) == sign(fr.cons_h[j]) ? fr.cons_h[j] : 0.0 for j in idx]
    gov_to(fr) = gov_from(fr) .+ fr.cons_g[idx]
    left_lim = fixed_limits(vcat(0.0, (vcat(heights(fr), fr.cons_h[idx], gov_to(fr)) for fr in frames)...))
    # Headroom above the tallest bar for the five-row legend.
    left_lim = (left_lim[1], left_lim[2] + 0.4 * (left_lim[2] - left_lim[1]))
    # Right panel: fixed limits over all frames and sectors.
    allq = vcat((pct(fr.cell.quantities) for fr in frames)...)
    allp = vcat((pct(fr.cell.prices) for fr in frames)...)
    xlim, ylim = fixed_limits(vcat(0.0, allq)), fixed_limits(vcat(0.0, allp))
    # Room on the right for the labels of the far-right movers.
    xlim = (xlim[1], xlim[2] + 0.15 * (xlim[2] - xlim[1]))
    inprog = Set(idx)
    others = [j for j in 1:n if !(j in inprog)]
    # Labelled movers, fixed by their largest motion over the frames.
    move(j) = maximum(max(abs(pct(f.cell.quantities[j])), abs(pct(f.cell.prices[j]))) for f in frames)
    tagged = vcat(sort(idx; by = move, rev = true)[1:min(3, length(idx))],
        sort(others; by = move, rev = true)[1:min(3, length(others))])

    k = M.Observable(1)
    fr = M.lift(i -> frames[i], k)
    fig = M.Figure(size = (1900, 980))
    M.Label(fig[1, 1:2], M.lift(f -> @sprintf("%s — %s = %+.3f", title_stub, sym, f.param), fr), fontsize = 22)
    M.Label(fig[2, 1:2], M.lift(f -> @sprintf(
        "ΔP (GDP deflator) %+.2f%%    ΔL %+.2f%%    ΔW (household consumption) %+.2f%%    Π/Y %.2f%%    %s ∈ [%.3f, %.3f]",
        100 * (f.defl - 1), 100 * (f.employment - 1), 100 * (f.consumption - 1),
        100 * f.profit_share, band_label, f.band_min, f.band_max), fr), fontsize = 16)
    axl = M.Axis(fig[3, 1]; xlabel = "Sector", ylabel = "%", ytickformat = "{:.2f}%",
        xticks = (1:length(idx), labels[idx]), xticklabelrotation = π / 4,
        xticklabelsize = 11, yticklabelsize = 11)
    M.hlines!(axl, [0.0]; color = :gray, linestyle = :dash)
    M.barplot!(axl, xs, M.lift(heights, fr); width = width, color = bar_colors)
    gov_fill = M.LinePattern(; direction = M.Vec2f(1, 1), width = 2.5, tilesize = (8, 8),
        linecolor = colors[3], backgroundcolor = :white)
    M.barplot!(axl, slot(3), M.lift(f -> f.cons_h[idx], fr); width = width, color = colors[3])
    M.barplot!(axl, slot(3), M.lift(gov_to, fr); fillto = M.lift(gov_from, fr), width = width,
        color = gov_fill, strokecolor = colors[3], strokewidth = 1)
    M.ylims!(axl, left_lim...)
    M.axislegend(axl, [M.PolyElement(polycolor = colors[1]), M.PolyElement(polycolor = colors[2]),
            M.PolyElement(polycolor = colors[3]),
            M.PolyElement(polycolor = gov_fill, strokecolor = colors[3], strokewidth = 1),
            M.PolyElement(polycolor = colors[4])],
        ["Programme demand (% of baseline gross output)", "Change in output",
            "Final consumption: household", "Final consumption: government (incl. programme, F2/F3)",
            "Price deviation from numeraire"];
        position = :rt, labelsize = 11)
    axr = M.Axis(fig[3, 2]; xlabel = "Change in quantities (%)", ylabel = "Change in prices (%)",
        xtickformat = "{:.2f}%", ytickformat = "{:.2f}%", xticklabelsize = 11, yticklabelsize = 11)
    M.hlines!(axr, [0.0]; color = :gray, linestyle = :dash)
    M.vlines!(axr, [0.0]; color = :gray, linestyle = :dash)
    # Ghost of the door-closed (parameter 0) frame.
    f0 = frames[findfirst(f -> iszero(f.param), frames)].cell
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
        ["Programme sectors", "Other sectors", "$sym = 0 (door closed)"]; position = legend_pos, labelsize = 12)
    anchor_txt = join([@sprintf("%s = %+.2f ↔ %s", sym, v, anchors[v]) for v in sort(collect(keys(anchors)))], ";  ")
    M.Label(fig[4, 1:2], M.lift(f -> haskey(anchors, f.param) ?
        "This frame reproduces the executed run $(anchors[f.param]) (validated)." :
        "Interpolating re-solve between executed runs — not a run, not citable (ADR-0025).", fr),
        fontsize = 12, color = :gray25)
    M.Label(fig[5, 1:2], _anim_note(sym, anchor_txt), fontsize = 10)
    M.rowsize!(fig.layout, 4, M.Fixed(18))
    M.rowsize!(fig.layout, 5, M.Fixed(18))
    # Hold the first and last frames, and briefly each validated anchor.
    seq = Int[]
    for (i, f) in enumerate(frames)
        reps = i == 1 ? fps : i == length(frames) ? 2 * fps :
               haskey(anchors, f.param) ? fps ÷ 2 : 1
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

_anim_note(sym, anchors) = "Shock effect relative to the design reference (no-programme baseline). " *
    "Frames: re-solves by continuation in $sym from $sym = 0, fixed axes. Validated against: " * anchors * "."

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
        # One cell per labour x financing [x s_pi]: the one with the largest parameter.
        best = Dict{Tuple{String,String,Float64},Tuple{Float64,String}}()
        for id in order
            c = design_cell(design_d, id)
            pc, _ = cell_pricing(design_d, c)
            pc isa Union{MarkupPricing,CapacityPricing} || continue
            key = (string(c["labor"]), string(c["financing"]), door_s_pi(pc))
            v = door_param(pc)
            (!haskey(best, key) || v > best[key][1]) && (best[key] = (v, id))
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
    kind = get(get(design_d, "pricing", Dict{String,Any}()), "kind", "none")
    vals = if kind == "markup"
        [o.kappa_max * s / o.steps for s in 0:o.steps]
    elseif kind == "capacity"
        iseven(o.steps) || throw(ArgumentError("--steps must be even for the symmetric delta grid"))
        [o.delta_max * (2s - o.steps) / o.steps for s in 0:o.steps]
    else
        throw(ArgumentError("design \"$(o.design)\" has no pricing door to animate"))
    end
    mkpath(o.outdir)
    written = String[]
    for id in ids
        cell = design_cell(design_d, id)
        pc, _ = cell_pricing(design_d, cell)
        println("$id: solving $(length(vals)) frames")
        frames, anchors = solve_frames(cell, o.design, design_d, ref, ψ, g, labels, vals; root = root)
        stub = pc isa MarkupPricing ?
            @sprintf("Markup door (P2, %s) — %s-%s",
                pc.s_pi == 1 ? "rent retained" : @sprintf("rent distributed, s_π = %.4g", pc.s_pi),
                cell["labor"], cell["financing"]) :
            @sprintf("Capacity door (P1, δ > 0 capacity pressure, δ < 0 Kaldor–Verdoorn) — %s-%s",
                cell["labor"], cell["financing"])
        base = joinpath(o.outdir, @sprintf("%s_%s_%s%s_%s", o.design, cell["labor"],
            cell["financing"], door_s_pi(pc) == 1 ? "" : @sprintf("_spi%.4g", door_s_pi(pc)),
            door_name(pc)))
        paths = Base.invokelatest(record_cell, M, frames, anchors, g, ref.sol.quantities, stub, base;
            fps = o.fps, formats = o.formats, sym = door_symbol(pc), band_label = door_band_label(pc),
            # The capacity sweep runs through the lower right (δ < 0) too.
            legend_pos = pc isa CapacityPricing ? :lt : :rb)
        foreach(p -> println("  wrote $p"), paths)
        append!(written, paths)
    end
    return written
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    anim_main()
end
