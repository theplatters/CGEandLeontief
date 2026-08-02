"""
    BFReplication.jl — Main module for the B&F (2019) replication

Baqaee & Farhi (2019), "The Macroeconomic Impact of Microeconomic Shocks:
Beyond Hulten's Theorem" (Econometrica, 87(4), 1155–1203)
"""
module BFReplication

include("data_loader.jl")
include("model.jl")
include("inflation_analysis.jl")

export DataLoader, BFModel, InflationAnalysis

end # module
