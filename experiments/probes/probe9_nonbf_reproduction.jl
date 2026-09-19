# experiments/probes/probe9_nonbf_reproduction.jl
#
# Focused reproduction check: the twelve non-BF v5 cells vs their committed
# manifests. Used to separate "the eta = 0 promotion moved GAMMA/DELTA" from
# "the near-singular fixed-wage system drifts with a non-reproducible warm
# start" by running it with and without the ADR-0020 kernel changes.

include(joinpath(@__DIR__, "..", "run.jl"))
using Printf, TOML

const NONBF = ("ALPHA", "BETA", "GAMMA", "DELTA")

function nonbf_main()
    root = default_root()
    design_d = load_design("matrix_5x3_v5"; root = root)
    ref = build_reference(design_d; root = root)
    data = ref.data
    N = length(data.factor_share)
    ψ, g = programme_vectors(design_d, data, collect(1:N); root = root)
    @printf("reference real_consumption = %.17g\n", real_consumption(ref.sol))
    keys = ("gdp_rel", "consumption_rel", "employment", "wage", "external_transfer",
        "programme_financing", "external_position", "gdp_wedge")
    dkeys = ("canary_s", "canary_ixm", "canary_diff")
    for labor in NONBF, fin in ("F1", "F2", "F3")
        id = "matrix_5x3-v5-$labor-$fin"
        cell = design_d["cells"][id]
        model = build_cell_model(cell, design_d, data, ψ, g)
        sol = solve_cell(cell, design_d, data, ψ, g, ref.init_warm)
        ev = evaluate_gates(cell, design_d, model, sol, ref.sol)
        man = TOML.parsefile(joinpath(root, "runs", id, "manifest.toml"))
        d, dk = 0.0, ""
        for k in keys
            v = abs(ev.metrics[k] - man["metrics"][k]); v > d && (d = v; dk = k)
        end
        for k in dkeys
            v = abs(ev.diagnostics[k] - man["diagnostics"][k]); v > d && (d = v; dk = k)
        end
        @printf("  %-12s %.3e (%s)   resid %.2e vs %.2e\n", "$labor-$fin", d, dk,
            ev.gates["residual"]["value"], man["gates"]["residual"]["value"])
    end
    return nothing
end

nonbf_main()
