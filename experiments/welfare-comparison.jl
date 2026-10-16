# experiments/welfare-comparison.jl — SUPERSEDED (kept as the exploratory record)
#
# This script predates the operator's ADOPTED-PREFERENCES decision (2026-10-16):
# it assigns the public part its baseline share (w_new ∝ c_0*(alpha/omega) + g_0)
# and therefore collapses W_CES* onto W_CES under F2/F3. The table generator
# experiments/make_table4_alt.jl now SUBSUMES g into the private part, so W_CES*
# differs from W_CES in every financing. Numbers here differ from the table by
# construction; do not cite them.
#
# Original header follows.
#
# the Table 5 comparison columns, all 15 cells
#
#   W_C     household-consumption welfare (the Table 4 column; from the manifest)
#   W_D     quantity index of demand: the own-weights bundle aggregate over c+g
#   W_CES   CES evaluation of the bundle at BASE preferences  -> prices the
#           composition shift the programme imposes (sigma = 0.9)
#   W_CES*  CES evaluation at NEW preferences against the old bundle at old
#           preferences, i.e. U_new(f_1) / U_old(f_0)
#   plus the first-order leisure-augmented variants.
#
# Preference weights are RECOVERED from the solutions rather than from the tilt
# formula: at any price vector the demand is c_i ∝ omega_i p_i^(-sigma), so
# omega_i ∝ c_i p_i^sigma. The bundle weights extend them to the public part at
# its baseline share, w_new ∝ c_0*(alpha/omega) + g_0 and w_old ∝ f_0. This
# collapses to W_CES whenever preferences do not move (F2/F3: alpha = omega),
# which is the reduction the specification requires.
#
#   julia --project=. experiments/welfare-comparison.jl

using OpenCGE, Printf

const AbstractPricing = OpenCGE.AbstractPricingClosure
isdefined(Main, :load_design) || include(joinpath(@__DIR__, "run.jl"))

const SIGMA = 0.9

"""CES aggregate of q under weights w (both q0-relative), w normalised inside."""
function U(q, q0, w)
    # common support: c_i = 0 in many sectors and the CES needs a negative
    # exponent to be finite (same guard the kernel's ces_quantity_index uses)
    k = (q .> 0) .& (q0 .> 0)
    qk, q0k, wk = q[k], q0[k], w[k]
    wk = wk ./ sum(wk)
    r = qk ./ q0k
    # the kernel's index parameterisation (ces_quantity_index): weights enter
    # linearly, NOT raised to 1/sigma -- the two forms differ by w <-> w^sigma
    # and only this one reproduces W_CES = W_b.
    return (sum(wk .* r .^ ((SIGMA - 1) / SIGMA)))^(SIGMA / (SIGMA - 1))
end

function main()
    root = default_root()
    dd = load_design("matrix_5x3_v10"; root = root)
    ref = build_reference(dd; root = root)
    drops = Vector{Int}(dd["data"]["drops"])
    kept = sort(setdiff(1:(length(ref.data.factor_share) + length(drops)), drops))
    ψ, g = programme_vectors(dd, ref.data, kept; root = root)

    bundle(sol; lam = 1.0) = final_demand_bundle(sol.model, sol.prices_raw,
        sol.quantities, OpenCGE._bundle_wage(sol); lambda = lam,
        external_transfer = sol.external_transfer)
    b0 = bundle(ref.sol)
    # baseline preference weights and the private/public split of the bundle
    ω = b0.c .* ref.sol.prices .^ SIGMA; ω = ω ./ sum(ω)
    f0 = b0.c .+ b0.g
    v0 = sum(ref.sol.prices .* f0)
    w_old = f0

    @printf("%-4s %-6s %8s %8s %8s %8s | %8s %8s %8s %8s\n", "fin", "labour",
        "W_C", "W_D", "W_CES", "W_CES*", "CES*L leak", "CES*L def", "CESL leak", "CESL def")
    for fin in ["F1", "F2", "F3"], lab in ["BF", "ALPHA", "BETA", "GAMMA", "DELTA"]
        id = "matrix_5x3-v10-$lab-$fin" * (lab == "BETA" ? "-etas05" : "")
        sol = solve_cell(design_cell(dd, id), dd, ref.data, ψ, g, ref.init_warm)
        b = bundle(sol)
        α = b.c .* sol.prices .^ SIGMA; α = α ./ sum(α)
        w_new = b0.c .* (α ./ ω) .+ b0.g
        f, f1 = b.c .+ b.g, sol.prices
        w_c = 100 * TOML.parsefile(joinpath(root, "runs", id, "manifest.toml"))["metrics"]["consumption_rel"]
        w_d = 100 * (sum(f) / sum(f0) - 1)
        w_ces = 100 * (U(f, f0, w_old) - 1)
        w_star = 100 * (U(f, f0, w_new) / U(f0, f0, w_old) - 1)
        # leisure: (ratio - 1) - w*dN/V0, private part stretched under :deferred
        ls, ls0 = sum(b.L_i), sum(b0.L_i)
        wage = b.wage_bill / ls
        s = sol.model.data.saving_rate
        lc(rat, stretch) = (rat - 1) - wage * (ls - ls0) / sum(ref.sol.prices .* (stretch .* f0))
        r_ces = U(f, f0, w_old); r_star = U(f, f0, w_new) / U(f0, f0, w_old)
        @printf("%-4s %-6s %7.3f%% %7.3f%% %7.3f%% %7.3f%% | %7.3f%% %7.3f%% %7.3f%% %7.3f%%\n",
            fin, lab, w_c, w_d, w_ces, w_star,
            100 * lc(r_star, 1.0), 100 * lc(r_star, 1 / (1 - s)),
            100 * lc(r_ces, 1.0), 100 * lc(r_ces, 1 / (1 - s)))
    end
    return 0
end

exit(main())
