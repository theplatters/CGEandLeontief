# =====================================================================
# Table 4 (ALT) --- the WELFARE-COMPARISON table (GENERATED)
# ---------------------------------------------------------------------
# Do NOT edit by hand. Regenerate with:
#   julia --project=<OPEN_CGE_CHECKOUT> experiments/make_table4_alt.jl.jl
#
#   NOT `--project=.`: this script calls kernel functions added in OpenCGE
#   09e2b83 (diag(WELFARE)), which the paper repo's Manifest does not yet pin.
#   Point --project at a checkout that contains them, or the load fails.
#
# Job: compare the welfare reads of one and the same allocation. The macro
# outcomes and the indices already shown in Table 4 are dropped; only W_C is
# carried over, so the table reads as a welfare comparison.
#
#   W_C      household-consumption welfare (the Table 4 column, from the manifest)
#   W_D      quantity index of demand: the own-weights aggregate over c + g
#   W_CES    CES evaluation of the bundle at the BASELINE preference weights
#   W_CES*   CES evaluation of the bundle at the ADOPTED preference weights
#   W_CES*^L the same, first-order leisure-augmented, under both saving conventions
#
# Preferences are ADOPTED throughout: the weights of a state are read off that
# state's own bundle, alpha_i ∝ f_i p_i^sigma, with g SUBSUMED into the private
# part. F2/F3 therefore adopt preferences too --- the household is treated as
# having taken the programme into its own demand --- which is why W_CES* differs
# from W_CES in every financing rather than only under F1.
#
# The ALT layout keeps the matrix grammar with the elastic-supply row
# bracketed over the executed elasticities, as in Table 4.
# =====================================================================

using OpenCGE, TOML, Printf

if !isdefined(OpenCGE, :AbstractPricing)   # the pricing-door line renames the type;
    @eval const AbstractPricing = OpenCGE.AbstractPricingClosure   # the release exports it
end
isdefined(Main, :load_design) || include(joinpath(@__DIR__, "run.jl"))

const SIGMA = 0.9
const LAMBDA = 1.0
const GATE_TOL = 1e-6

const FIN = ["F1", "F2", "F3"]
const LABOUR = ["BF", "ALPHA", "BETA", "GAMMA", "DELTA"]
const ETAS = ["-etas05", "-etas1", "-etas2"]          # the bracket, low end first
const ETAS_LO, ETAS_HI = "0.5", "2"

const FINLABEL = Dict("F1" => "repurposing of consumption",
                      "F2" => "refinancing via taxes",
                      "F3" => "external / debt finance")
const LABEL = Dict("BF" => "No-reallocation", "ALPHA" => "Standard CGE",
                   "BETA" => "Elastic supply", "GAMMA" => "Fixed real wage",
                   "DELTA" => "Leontief")

cell_id(lab, fin, suf = "") = "matrix_5x3-v10-$lab-$fin" * (lab == "BETA" ? (suf == "" ? "-etas05" : suf) : "")

"""The bundle block: c, g, L_i and wage_bill at one solution."""
function bundle_block(sol::Solution; lam::Real = LAMBDA)
    final_demand_bundle(sol.model, sol.prices_raw, sol.quantities,
        OpenCGE._bundle_wage(sol); lambda = lam, external_transfer = sol.external_transfer)
end

"""The bundle vector f = stretch*c + lambda*g."""
fvec(b, stretch::Real; lam::Real = LAMBDA) = stretch .* b.c .+ lam .* b.g

"""
CES aggregate of q against q0 under weights w, in the kernel's index
parameterisation (`ces_quantity_index`): weights enter linearly, not as
w^(1/sigma). Support-restricted, since c_i = 0 in many of the 62 sectors.
"""
function U_idx(q::AbstractVector, q0::AbstractVector, w::AbstractVector)
    k = (q .> 0) .& (q0 .> 0)
    qk, q0k, wk = q[k], q0[k], w[k]
    wk = wk ./ sum(wk)
    r = qk ./ q0k
    return (sum(wk .* r .^ ((SIGMA - 1) / SIGMA)))^(SIGMA / (SIGMA - 1))
end

"""The five computed reads for one cell, under both saving conventions."""
function reads(sol::Solution, ref)
    b1, b0 = bundle_block(sol), bundle_block(ref.sol)
    ls, ls0 = sum(b1.L_i), sum(b0.L_i)
    w = b1.wage_bill / ls
    out = Dict{Symbol,Float64}()
    for (key, stretch) in ((:leak, 1.0), (:def, 1.0 / (1 - sol.model.data.saving_rate)))
        f, f0 = fvec(b1, stretch), fvec(b0, stretch)
        alpha = f .* sol.prices .^ SIGMA        # the state's own, ADOPTED weights
        r_ces, r_star = U_idx(f, f0, f0), U_idx(f, f0, alpha)
        v0 = sum(ref.sol.prices .* f0)
        out[Symbol(:ces, key)] = 100 * (r_ces - 1)
        out[Symbol(:star, key)] = 100 * (r_star - 1)
        out[Symbol(:starl, key)] = 100 * ((r_star - 1) - w * (ls - ls0) / v0)
    end
    out[:d] = 100 * (sum(fvec(b1, 1.0)) / sum(fvec(b0, 1.0)) - 1)     # W_D
    return out
end

function main()
    root = default_root()
    dd = load_design("matrix_5x3_v10"; root = root)
    @printf("building reference continuation for matrix_5x3_v10 ...\n")
    ref = build_reference(dd; root = root)
    drops = Vector{Int}(dd["data"]["drops"])
    kept = sort(setdiff(1:(length(ref.data.factor_share) + length(drops)), drops))
    ψ, g = programme_vectors(dd, ref.data, kept; root = root)

    rows = NamedTuple[]
    for fin in FIN, lab in LABOUR
        sufs = lab == "BETA" ? ETAS : [""]
        got, okall, wc = NamedTuple[], true, NaN
        for suf in sufs
            id = cell_id(lab, fin, suf)
            manpath = joinpath(root, "runs", id, "manifest.toml")
            man = isfile(manpath) ? TOML.parsefile(manpath) : nothing
            sol = solve_cell(design_cell(dd, id), dd, ref.data, ψ, g, ref.init_warm)
            r = reads(sol, ref)
            wc = man === nothing ? NaN : 100 * man["metrics"]["consumption_rel"]
            # canary: the baseline-weight CES must reproduce the kernel's index
            ker = real_final_demand_ces(sol, ref.sol; lambda = LAMBDA, sigma = SIGMA) - 1
            abs(r[:cesleak] - 100 * ker) < 1e-8 || @warn "W_CES differs from the kernel" id r[:cesleak] ker
            ok = man !== nothing && abs(wc - 100 * (real_consumption(sol) / real_consumption(ref.sol) - 1)) < GATE_TOL
            push!(got, (; wc = wc, ok = ok, d = r[:d], cesleak = r[:cesleak],
                starleak = r[:starleak], starlleak = r[:starlleak], starldef = r[:starldef]))
        end
        okall = all(x -> x.ok, got)
        lo, hi = x -> minimum(x), x -> maximum(x)
        agg(f) = (lo([f(x) for x in got]), hi([f(x) for x in got]))
        push!(rows, (lab = lab, fin = fin, ok = okall,
            wc = agg(x -> x.wc), wd = agg(x -> x[:d]),
            ces = agg(x -> x[:cesleak]), star = agg(x -> x[:starleak]),
            sl = agg(x -> x[:starlleak]), sd = agg(x -> x[:starldef])))
        @printf("  %-26s gate %s  Wc %+7.3f  Wd %+7.3f  CES %+7.3f  CES* %+7.3f  L %+7.3f/%+7.3f\n",
            cell_id(lab, fin), okall ? "ok" : "FAIL", got[1].wc, got[1][:d],
            got[1][:cesleak], got[1][:starleak], got[1][:starlleak], got[1][:starldef])
    end
    path = joinpath(root, "revision", "tables", "table4_results-5x3-ALT.tex")
    open(path, "w") do io
        render(io, rows)
    end
    println("wrote ", path)
    return all(r -> r.ok, rows) ? 0 : 1
end

# same convention as table4_results-5x3.tex: values in math mode, so the
# unary minus is a proper minus sign and "--" stays a text en-dash between them
num(x) = (!isfinite(x)) ? "n/a" : (abs(x) < 5e-5 ? "\$0.00\$" :
         (x > 0 ? @sprintf("\$+%.2f\$", x) : @sprintf("\$%.2f\$", x)))

"""min--max interval from a (lo, hi) pair; a point where the bracket is flat."""
function rngv(p)
    lo, hi = p
    abs(hi - lo) < 5e-5 && return num(lo)
    return num(lo) * "--" * num(hi)
end

# the variables shown, in the order they appear inside each cell
const MEAS = [(:wc,   "\$W_C\$"),
              (:wd,   "\$W_D\$"),
              (:ces,  "\$W_{CES}\$"),
              (:star, "\$W_{CES^*}\$"),
              (:sl,   "\$W^{L}\$")]

CELL_SPEC = raw"@{}l@{\hspace{3pt}}r@{}"
# elastic-supply cell: label column + one bracketed column
ETAS_CELL_SPEC = raw"@{}>{\raggedright\arraybackslash}p{1.05cm}@{\hspace{3pt}}>{\raggedleft\arraybackslash}p{2.2cm}@{}"

"""One plain cell: the six welfare variables stacked, as in Table 4."""
function plain_cell(io, r)
    println(io, "\\begin{tabular}{" * CELL_SPEC * "}")
    for (k, lab) in MEAS
        println(io, lab * " & " * num(getfield(r, k)[1]) * " \\\\")
    end
    print(io, "\\end{tabular}")
end

"""The elastic-supply cell: same variables, each as a min--max interval."""
function etas_cell(io, r)
    println(io, "\\begin{tabular}{" * ETAS_CELL_SPEC * "}")
    for (k, lab) in MEAS
        println(io, lab * " & " * rngv(getfield(r, k)) * " \\\\")
    end
    print(io, "\\end{tabular}")
end

"""The elasticity bracket, printed in the elastic-supply column header."""
function etas_head()
    return ("{\\footnotesize\\begin{tabular}{" * ETAS_CELL_SPEC * "}\n"
            * "& \$\\eta_{s,i} \\in [" * ETAS_LO * "," * ETAS_HI * "]\$\\end{tabular}}")
end

function render(io, rows)
    by = Dict((r.fin, r.lab) => r for r in rows)
    println(io, "% =====================================================================")
    println(io, "% Table 4 (ALT) --- welfare comparison across the closure matrix (GENERATED)")
    println(io, "% ---------------------------------------------------------------------")
    println(io, "% Do NOT edit by hand. Regenerate with:")
    println(io, "%   julia --project=<OPEN_CGE_CHECKOUT> experiments/make_table4_alt.jl")
    println(io, "% Layout mirrors table4_results-5x3.tex exactly --- same axes, same matrix")
    println(io, "%   grammar --- with the welfare variables in place of the macro outcomes.")
    println(io, "% W_C is the manifest field; W_D, W_CES and W_CES* are computed here and")
    println(io, "%   gated against runs/<id>/manifest.toml (consumption_rel, tol 1e-6). No run")
    println(io, "%   directory is created or modified, and no number is hand-typed.")
    println(io, "% =====================================================================")
    println(io, "\\begin{table}[htbp]")
    println(io, "\\centering")
    println(io, "\\footnotesize")
    println(io, "\\setlength{\\tabcolsep}{3pt}")
    println(io, "\\setlength{\\parindent}{0pt}")
    println(io, "\\setlength{\\parskip}{2pt}")
    println(io, "\\renewcommand{\\arraystretch}{1.25}")
    width = Dict("BETA" => "3.5cm")
    print(io, "\\begin{tabular}{@{}>{\\raggedright\\arraybackslash}p{1.5cm}")
    for lab in LABOUR
        print(io, " >{\\raggedright\\arraybackslash}p{" * get(width, lab, "1.95cm") * "}")
    end
    println(io, "@{}}")
    println(io, "\\toprule")
    println(io, "\\textbf{Financing closure}")
    println(io, "& \\textbf{No-reallocation}\\newline \\scriptsize immobile labour, sectoral wages")
    println(io, "& \\textbf{Standard CGE}\\newline \\scriptsize one flexible wage, mobile labour")
    println(io, "& \\textbf{Elastic labour supply}\\newline \\scriptsize sectoral supply on the real wage\\newline " * etas_head())
    println(io, "& \\textbf{Fixed real wage}\\newline \\scriptsize real wage pegged, employment uncapped")
    println(io, "& \\textbf{Leontief closure}\\newline \\scriptsize real wage pegged, quantities absorb")
    println(io, "\\\\")
    println(io, "\\midrule")
    for (i, fin) in enumerate(FIN)
        flag = all(by[(fin, lab)].ok for lab in LABOUR) ? "" : "\\textsuperscript{\\dag}"
        println(io, "\\textbf{" * fin * "}\\newline \\scriptsize " * FINLABEL[fin] * flag)
        for lab in LABOUR
            r = by[(fin, lab)]
            print(io, "& ")
            lab == "BETA" ? etas_cell(io, r) : plain_cell(io, r)
        end
        println(io, "\\\\")
        i < length(FIN) && println(io, "\\midrule")
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
    println(io, "\\caption{Welfare comparison across the closure matrix. The layout mirrors")
    println(io, "Table~\\ref{tab:results-matrix} exactly --- the same three financing closures")
    println(io, "(rows) crossed with the same five labour-market closures (columns), with the")
    println(io, "elastic-supply column bracketed the same way --- but the variables reported in")
    println(io, "each cell are welfare reads of the identical allocation rather than macro")
    println(io, "outcomes, all in percent against the reference continuation. \$W_C\$ is the")
    println(io, "household-consumption index of Table~\\ref{tab:results-matrix}; \$W_D\$ the")
    println(io, "quantity index of demand, i.e.\\ the own-weights aggregate over \$c+g\$;")
    println(io, "\$W_{CES}\$ the CES evaluation of that bundle at the \\emph{baseline}")
    println(io, "preference weights, which is what prices the composition shift the programme")
    println(io, "imposes; \$W_{CES^*}\$ the same bundle evaluated at the \\emph{adopted} weights,")
    println(io, "read off the counterfactual bundle itself with \$g\$ subsumed into the private")
    println(io, "part; and \$W^{L}\$ the first-order leisure-augmented form, which also")
    println(io, "charges the household for the additional work a programme requires.")
    println(io, "It is expressed on the household's resources rather than on a consumption basket,")
    println(io, "so the treatment of saving has to be stated; the conservative reading, in which")
    println(io, "saving is a leak that never returns, is the one tabulated here, and the deferred")
    println(io, "alternative is discussed in the text.")
    println(io, "\$\\eta_{s,i}\\in\\{" * ETAS_LO * ",1," * ETAS_HI * "\\}\$, as in")
    println(io, "Table~\\ref{tab:results-matrix}.")
    println(io, "}")
    println(io, "\\label{tab:welfare-comparison}")
    println(io, "\\end{table}")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    exit(main())
end
