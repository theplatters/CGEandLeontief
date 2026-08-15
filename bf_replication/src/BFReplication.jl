
# BFReplication.jl -- top-level module for the B&F (2019) replication package.
#
# Structure:
#   DataLoader          R1  -- load US 1980 IO tables (BFdata.csv)
#   BFModel             R2  -- pedagogical equilibrium solver (numerical Jacobian)
#   InflationAnalysis        -- the "no inflation" decomposition module
#   BFCore                   -- PRODUCTION solver (analytical Jacobian, port of
#                               Simulation.m + Simulation_Derivs.m) for R3/R4
#   BFShocks                 -- TFP shock generation (R3)
#   BFMonteCarlo             -- 50k-draw robustness Monte Carlo (R3)
#   BFGradient               -- elasticity-gradient sweep (R4)

module BFReplication

include("data_loader.jl")
include("model.jl")
include("inflation_analysis.jl")
include("core_solver.jl")
include("shocks.jl")
include("monte_carlo.jl")
include("elasticity_gradient.jl")

export DataLoader, BFModel, InflationAnalysis,
       gdp_nominal, gdp_welfare, mean_price, solve_bf, good_init,
       real_gdp_mc,
       empirical_covariances, draw_logA, mvnrnd,
       run_monte_carlo, moments_loggdp,
       run_elasticity_gradient

end # module
