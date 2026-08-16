# BFReplication.jl — top-level module for the B&F (2019) replication package.
#
# Phases:
#   DataLoader          R1  — load US 1980 IO tables (BFdata.csv)
#   BFModel             R2  — pedagogical equilibrium solver
#   InflationAnalysis        — no-inflation decomposition
#   core_solver.jl      R2–R5  — production solver (fixed + mobile labor)
#   shocks.jl           R3  — TFP shock generation
#   monte_carlo.jl      R3+R5 — fixed-labor + reallocation MC
#   elasticity_gradient.jl  R4 — elasticity sweep
#   second_order.jl     R6  — second-order approximation

module BFReplication

include("data_loader.jl")
include("model.jl")
include("inflation_analysis.jl")
include("core_solver.jl")
include("shocks.jl")
include("monte_carlo.jl")
include("elasticity_gradient.jl")
include("second_order.jl")
include("oil_shock.jl")

export DataLoader, BFModel, InflationAnalysis,
       gdp_nominal, gdp_welfare, mean_price, solve_bf, good_init,
       solve_bf_realloc, BFReallocSolution,
       real_gdp_mc,
       empirical_covariances, draw_logA, mvnrnd,
       run_monte_carlo, run_monte_carlo_realloc, moments_loggdp,
       run_elasticity_gradient,
       second_order_hessian_norealloc, second_order_hessian_realloc,
       second_order_mc, run_second_order_norealloc, run_second_order_realloc,
       run_oil_shock

end # module