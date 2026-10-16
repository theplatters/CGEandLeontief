# experiments/make_table4_alt.jl — Table 4 (ALT): the three welfare measures
#
# Re-solves the fifteen matrix cells of `matrix_5x3_v10` IN MEMORY through the
# harness and reports, per cell, the four macro outcomes plus the three welfare
# measures of Section 4.3:
#
#   W_a    Tornqvist quantity index of consumption      (`real_consumption`)
#   W_a'   exact money-metric of consumption            (`real_consumption_ev`)
#   W_b    CES quantity index of the c+g bundle         (`real_final_demand_ces`)
#   W_c    first-order leisure-augmented measure        (`welfare_leisure`),
#          under both saving conventions
#
# Nothing is hand-typed. The driver NEVER creates run directories and never
# writes to `runs/`: it only reads the design and the recorded manifests, which
# it also uses as a validity gate — each re-solved cell must reproduce the
# recorded `consumption_rel` to `#GATE_TOL`, otherwise the row is marked and the
# run reports it.
#
#   julia --project=. experiments/make_table4_alt.jl
#   julia --project=. experiments/make_table4_alt.jl --design matrix_5x3_v10 --out <path>

using OpenCGE, TOML, Printf

const GATE_TOL = 1e-6            # re-solve vs recorded consumption_rel
const SIGMA = 0.9                # decision D3
const LAMBDA = 1.0               # decision D2

# Compatibility shim (2026-10-16): the working-copy kernel renamed the abstract
# pricing type to `AbstractPricingClosure`, while `experiments/run.jl` (and
# `animate_pricing_panels.jl`) still refer to it as `AbstractPricing`, so the
# harness does not load against the checkout. This alias restores the old name
# for THIS script only; drop it once the harness is updated. Nothing else in
# `run.jl` is touched, and the reproduction gate below would expose any deeper
# drift between the working copy and the kernel that minted the v10 runs.
const AbstractPricing = OpenCGE.AbstractPricingClosure

if !isdefined(Main, :load_design)
    include(joinpath(@__DIR__, "run.jl"))
end

function _parse(args)
    design, out = "matrix_5x3_v10", nothing
    i = 1
    while i <= length(args)
        if args[i] == "--design"; design = args[i + 1]; i += 2
        elseif args[i] == "--out"; out = args[i + 1]; i += 2
        else; i += 1
        end
    end
    return design, out
end

const LABOUR = ["BF", "ALPHA", "BETA", "GAMMA", "DELTA"]
const FIN = ["F1", "F2", "F3"]
const LABEL = Dict("BF" => "No-reallocation", "ALPHA" => "Standard CGE",
    "BETA" => "Elastic supply", "GAMMA" => "Fixed real wage",
    "DELTA" => "Leontief")
const FINLABEL = Dict("F1" => "repurposing", "F2" => "taxes", "F3" => "external")

cell_id(lab, fin) = "matrix_5x3-v10-$lab-$fin" * (lab == "BETA" ? "-etas05" : "")

"""
	own_bundle(sol, base) -> Float64

Own-weights bundle aggregate (spec v6): each state's `c + lambda*g` bundle valued
by its own composition weights, which reduces identically to the ratio of bundle
quantity sums. Equals the bundle's value ratio where prices are at base level.
"""
function own_bundle(sol::Solution, base::Solution)::Float64
	f(s) = final_demand_bundle(s.model, s.prices_raw, s.quantities,
		OpenCGE._bundle_wage(s); lambda = LAMBDA, external_transfer = s.external_transfer).f
	return sum(f(sol)) / sum(f(base))
end

function table4_alt_main(args = ARGS)
    design, out = _parse(args)
    root = default_root()
    dd = load_design(design; root = root)
    @printf("building reference continuation for %s ...\n", design)
    ref = build_reference(dd; root = root)
    drops = Vector{Int}(dd["data"]["drops"])
    n_full = length(ref.data.factor_share) + length(drops)
    kept = sort(setdiff(1:n_full, drops))
    ψ, g = programme_vectors(dd, ref.data, kept; root = root)

    rows = NamedTuple[]
    for fin in FIN, lab in LABOUR
        id = cell_id(lab, fin)
        cell = design_cell(dd, id)
        manpath = joinpath(root, "runs", id, "manifest.toml")
        man = isfile(manpath) ? TOML.parsefile(manpath) : nothing
        sol = solve_cell(cell, dd, ref.data, ψ, g, ref.init_warm)
        w_a = 100 * (real_consumption(sol) / real_consumption(ref.sol) - 1)
        w_ap = 100 * (real_consumption_ev(sol, ref.sol) - 1)
        w_b = 100 * (real_final_demand_ces(sol, ref.sol; lambda = LAMBDA, sigma = SIGMA) - 1)
        w_bo = 100 * (own_bundle(sol, ref.sol) - 1)
        # W_c hangs off the SAME bundle as the bundle column: own weights (spec v6),
        # so the leisure term rides on a regime-consistent aggregate.
        lk = welfare_leisure(sol, ref.sol; lambda = LAMBDA, sigma = SIGMA,
            saving = :leak, aggregator = :ownweights)
        df = welfare_leisure(sol, ref.sol; lambda = LAMBDA, sigma = SIGMA,
            saving = :deferred, aggregator = :ownweights)
        rec = man === nothing ? NaN : 100 * man["metrics"]["consumption_rel"]
        ok = isfinite(rec) ? abs(w_a - rec) < GATE_TOL : false
        push!(rows, (lab = lab, fin = fin, id = id, ok = ok,
            dY = man === nothing ? NaN : 100 * man["metrics"]["gdp_rel"],
            dL = man === nothing ? NaN : 100 * (man["metrics"]["employment"] - 1),
            dP = man === nothing ? NaN : 100 * (man["metrics"]["gdp_deflator"] - 1),
            dD = man === nothing ? NaN : 100 * man["metrics"]["external_position"],
            wa = w_a, wap = w_ap, wbo = w_bo, wb = w_b,
            wcl = 100 * lk.welfare, wcd = 100 * df.welfare,
            q = 100 * (lk.value_ratio - 1)))
        @printf("  %-22s gate %s  Wa %+7.3f  Wap %+7.3f  Bown %+7.3f  Wb %+7.3f  Wc %+7.3f/%+7.3f\n",
            id, ok ? "ok" : "FAIL", w_a, w_ap, w_bo, w_b, 100 * lk.welfare, 100 * df.welfare)
    end
    path = out === nothing ? joinpath(root, "revision", "tables", "table4_results-5x3-ALT.tex") : out
    open(path, "w") do io
        render(io, rows, design)
    end
    println("wrote ", path)
    return all(r -> r.ok, rows) ? 0 : 1
end

num(x) = (!isfinite(x)) ? "n/a" : (abs(x) < 5e-5 ? "0.00" : @sprintf("%+.2f", x))
num4(x) = (!isfinite(x)) ? "n/a" : (abs(x) < 5e-5 ? "0.000" : @sprintf("%+.3f", x))

function render(io, rows, design)
    println(io, "% =====================================================================")
    println(io, "% Table 4 (ALT) --- the closure matrix with the three welfare measures")
    println(io, "% ---------------------------------------------------------------------")
    println(io, "% Do NOT edit by hand. Regenerate with:")
    println(io, "%   julia --project=. experiments/make_table4_alt.jl")
    println(io, "% Source: the ", design, " design, re-solved in memory through the harness and")
    println(io, "%   gated against runs/<id>/manifest.toml (consumption_rel, tol 1e-6). The")
    println(io, "%   macro columns are the recorded manifest fields; the welfare columns are")
    println(io, "%   computed in the kernel by real_consumption / real_consumption_ev /")
    println(io, "%   real_final_demand_ces / welfare_leisure. No number is hand-typed and no")
    println(io, "%   run directory is created or modified.")
    println(io, "% All in percent against the reference continuation. B^own is the OWN-WEIGHTS")
    println(io, "%   bundle aggregate over c+lambda*g (lambda=1): each state under its own")
    println(io, "%   composition weights, which reduces to the bundle value ratio --- the")
    println(io, "%   regime-consistent measure (spec v6). W_b is the fixed-weight CES at")
    println(io, "%   sigma=0.9 --- NOT the value ratio --- and is unidentified under F1. W_c is")
    println(io, "%   the first-order leisure-augmented measure on the OWN-WEIGHTS bundle under")
    println(io, "%   the two saving conventions; the saving convention, not the leisure")
    println(io, "%   preference, sets the sign at the wage-pinned pole.")
    println(io, "% =====================================================================")
    println(io, "\\begin{table}[htbp]")
    println(io, "\\centering")
    println(io, "\\scriptsize")
    println(io, "\\setlength{\\tabcolsep}{4pt}")
    println(io, "\\renewcommand{\\arraystretch}{1.15}")
    print(io, "\\begin{tabular}{@{}>{\\raggedright\\arraybackslash}p{1.9cm}")
    for _ in 1:10
        print(io, ">{\\raggedleft\\arraybackslash}p{1.05cm}")
    end
    println(io, "@{}}")
    println(io, "\\toprule")
    println(io, "& \\multicolumn{4}{c}{\\textbf{Macro outcomes}} & \\multicolumn{6}{c}{\\textbf{Welfare measures}} \\\\")
    println(io, "\\cmidrule(lr){2-5}\\cmidrule(lr){6-11}")
    println(io, "\\textbf{Cell} & \$\\Delta Y\$ & \$\\Delta L\$ & \$\\Delta P\$ & \$\\Delta D\$ " *
                "& \$W_a\$ & \$W_a'\$ & \$B^{\\mathrm{own}}\$ & \$W_b\$ & \$W_c^{\\mathrm{leak}}\$ & \$W_c^{\\mathrm{def}}\$ \\\\")
    println(io, "\\midrule")
    lastfin = ""
    for r in rows
        if r.fin != lastfin
            r.fin == "F1" || println(io, "\\midrule")
            mark = r.fin == "F1" ? "\\textsuperscript{\\ddag}" : ""
            println(io, "\\multicolumn{11}{@{}l}{\\textit{Financing \\textbf{$(r.fin)}: $(FINLABEL[r.fin])}$mark}\\\\")
            lastfin = r.fin
        end
        flag = r.ok ? "" : "\\textsuperscript{\\dag}"
        # The bundle indices are NOT identified under F1: the preference tilt
        # replaces the aggregator's own weights, so the CES aggregate loses its
        # utility foundation and reads the intended reallocation as a loss
        # (measured: ~ -0.47 pp against the exact money-metric on the SAME
        # consumption basket, in every F1 cell; see welfare-f1-decomp.jl).
        wb = r.fin == "F1" ? "--" : num4(r.wb)
        println(io, "\\quad \\textbf{$(r.lab)} \\scriptsize $(LABEL[r.lab])$flag" *
                    " & $(num(r.dY)) & $(num(r.dL)) & $(num(r.dP)) & $(num(r.dD))" *
                    " & $(num4(r.wa)) & $(num4(r.wap)) & $(num4(r.wbo)) & $wb" *
                    " & $(num4(r.wcl)) & $(num4(r.wcd)) \\\\")
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
    println(io, "\\caption{The closure matrix with the three welfare measures. \$W_a\$ is")
    println(io, "the T\\\"ornqvist index of household consumption (the measure currently")
    println(io, "reported); \$W_a'\$ the exact money-metric of consumption, on the same base;")
    println(io, "\$W_b\$ the CES quantity aggregate over \$c+g\$ at \$\\lambda=1\$, \$\\sigma=0.9\$ ---")
    println(io, "which is \\emph{not} the value ratio, and differs from it at first order")
    println(io, "wherever the basket does not grow proportionally, as the programme cells")
    println(io, "do; \$B^{\\mathrm{own}}\$ the own-weights bundle aggregate --- each state's")
    println(io, "\$c+g\$ valued by its own composition weights, which reduces identically to")
    println(io, "\$\\sum_i f_i\$, the bundle value ratio at base prices; it is the")
    println(io, "preference-consistent bundle measure and carries the F1 rows, where the")
    println(io, "fixed-weight aggregates are not defined (see the dagger);")
    println(io, "\$W_c\$ that same own-weights bundle less the value of the")
    println(io, "additional work \$w\\,\\Delta N\$, to first order, under the two saving")
    println(io, "conventions; it therefore inherits the bundle column's regime-consistency,")
    println(io, "to first order, under the two saving conventions. \$W_a\$ and \$W_a'\$ coincide")
    println(io, "except under F1, whose preference shifter moves the basket off proportional")
    println(io, "growth. Macro columns: the recorded manifest fields, all in percent against")
    println(io, "the reference continuation. \\textsuperscript{\\ddag}~F1 rows report the")
    println(io, "consumption measures, \$B^{\\mathrm{own}}\$ and \$W_c\$, but not the fixed-weight")
    println(io, "aggregates. Under F1 the preference weights are themselves shifted, so the")
    println(io, "fixed-weight CES has no utility foundation: on the same consumption basket it")
    println(io, "sits about \$0.47\$ pp below the exact money-metric in every F1 cell, and")
    println(io, "would report the intended reallocation as a loss; the own-weights aggregate,")
    println(io, "being regime-consistent, is unaffected. \\texttt{experiments/welfare-f1-decomp.jl}")
    println(io, "and \\texttt{experiments/welfare-ownweights.jl} reproduce the two.")
    if !all(r -> r.ok, rows)
        println(io, "\\textsuperscript{\\dag}~cell did not reproduce through the in-memory solve within \$10^{-6}\$; regenerate and inspect.")
    end
    println(io, "}")
    println(io, "\\label{tab:results-matrix-alt}")
    println(io, "\\end{table}")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    exit(table4_alt_main())
end
