# tests/runtests.jl — paper-repo suite (post kernel split).
#
# The kernel package `OpenCGE` (separate repo) owns the synthetic-only tests;
# this suite covers the paper-side tooling and the real-data contracts that
# need data/, cbase2/data_raw/impulses.csv and runs/. Run it with
#   julia --project=. tests/runtests.jl
# (this repo is an environment, not a package: Pkg.test() does not apply).
using Test

include("test_helpers.jl")
include("test_calibration.jl")
include("test_gdp_measurement.jl")
include("test_external_closure.jl")
include("test_sectoral_labour.jl")
include("test_supply_arm_recovery.jl")
include("test_run_manifest.jl")
include("test_f1_tilt.jl")
include("test_check_repo.jl")
include("test_matrix_plots.jl")
