"""
    data_loader.jl — Data loading from CSV (no MAT.jl dependency)

Loads B&F (2019) data from the CSV exports:
  - BFdata.csv: raw Jorgenson 88-sector data (1960–2005)
  - stfp.csv: sectoral TFP growth rates (Carvalho-Gabaix method)

The processing pipeline replicates getData.m and the data setup in
GDP_Simulation_88sectorKLEMS.m.
"""
module DataLoader

using DelimitedFiles
using LinearAlgebra
using Statistics
using Printf

export BFData, load_bf_data, load_tfp_data, describe_data

"""
    BFData

Container for the processed B&F (2019) replication data.
"""
struct BFData
    Ω::Matrix{Float64}        # cost share matrix (N×N)
    α::Vector{Float64}        # factor (value-added) shares (N)
    β::Vector{Float64}        # consumption shares (N)
    λ::Vector{Float64}        # Domar weights (N)
    L::Vector{Float64}        # baseline labor allocation (N)
    N::Int                    # number of sectors
    year::Int                 # base year for IO matrix
    grossy::Matrix{Float64}   # nominal gross output (N × 46)
    vadd::Matrix{Float64}     # nominal value added (N × 46)
    capital::Matrix{Float64}  # nominal capital (N × 46)
    labor::Matrix{Float64}    # nominal labor (N × 46)
end

"""
    load_bf_data(csvpath; year=1980) -> BFData

Load and process the Jorgenson 88-sector data from BFdata.csv.

BFdata.csv format (matching us80dbasedata.mat):
  - 4048 rows: 88 sectors × 46 years each
  - Each row: [year, sector_id, grossy, capital, labor, intermediates...]

Processing steps (matching MATLAB getData.m / GDP_Simulation_88sectorKLEMS.m):
  1. Reshape into 88×46 sector-year matrices
  2. Remove government sectors (60, 80:88) and owner-occupied housing
  3. Remove zero-sales sectors (8=uranium ores, 62=renting machinery)
  4. Extract IO matrix for base year
  5. Compute Ω (row-normalized cost shares)
  6. Compute α = vadd / grossy (factor shares)
  7. Compute β = grossy' · (I - diag(1-α)·Ω), normalize
  8. Compute λ = (I - diag(1-α)·Ω)⁻¹' · β (Domar weights)
  9. Compute L = λ ⊙ α (baseline labor allocation)
"""
function load_bf_data(csvpath::String; year::Int=1980)
    raw = readdlm(csvpath, ',')

    n_sectors_raw = 88
    n_years = 46

    # Reshape: data stored as 88 blocks of 46 rows
    # Each block: [year, sector_id, grossy, capital, labor, intermediates...]
    grossy = zeros(n_sectors_raw, n_years)
    capital = zeros(n_sectors_raw, n_years)
    labor = zeros(n_sectors_raw, n_years)

    # BFdata.csv: column 1=year, column 2=sector_id, column 3=grossy,
    # column 4=capital, column 5=labor, columns 6:93=intermediate inputs
    for row in 1:size(raw, 1)
        yr = Int(raw[row, 1]) - 1959  # year index (1960→1)
        sec = Int(raw[row, 2])        # sector index (1:88)
        if 1 <= sec <= n_sectors_raw && 1 <= yr <= n_years
            grossy[sec, yr] = Float64(raw[row, 3])
            capital[sec, yr] = Float64(raw[row, 4])
            labor[sec, yr] = Float64(raw[row, 5])
        end
    end

    vadd = capital .+ labor

    # Remove government sectors and owner-occupied housing (60, 80:88)
    remove_sectors = [60; 80:88]
    keep = setdiff(1:n_sectors_raw, remove_sectors)
    grossy = grossy[keep, :]
    vadd = vadd[keep, :]
    capital = capital[keep, :]
    labor = labor[keep, :]

    # Remove zero-sales sectors
    zero_sales = findall(==(0.0), vec(sum(grossy, dims=2)))
    if !isempty(zero_sales)
        keep2 = setdiff(1:size(grossy, 1), zero_sales)
        grossy = grossy[keep2, :]
        vadd = vadd[keep2, :]
        capital = capital[keep2, :]
        labor = labor[keep2, :]
    end

    N = size(grossy, 1)
    yr_idx = year - 1959  # 1980 → 21

    # Extract IO matrix for base year from raw data
    # Columns 6:93 in BFdata.csv are intermediate inputs (88 sectors)
    # But after removing sectors, we need to build the IO table carefully

    # Build IO: find all rows for the base year
    year_rows = findall(==(Float64(year)), raw[:, 1])

    # Sort by sector id
    sort!(year_rows, by = i -> raw[i, 2])

    # Get sector IDs for the base year
    sec_ids = Int.(raw[year_rows, 2])

    # Build full IO table (88×88): row=buying sector, col=supplying sector
    io_full = zeros(n_sectors_raw, n_sectors_raw)
    for (idx, row) in enumerate(year_rows)
        sec = Int(raw[row, 2])
        # Columns 6:93 = intermediate inputs from sectors 1:88
        for j in 1:min(n_sectors_raw, size(raw, 2) - 5)
            io_full[sec, j] = Float64(raw[row, 5 + j])
        end
    end

    # Remove government and zero-sales sectors from IO
    io_reduced = io_full[keep, keep]
    # Also remove zero-sales from io_reduced
    if !isempty(zero_sales)
        keep2 = setdiff(1:size(io_reduced, 1), zero_sales)
        io_reduced = io_reduced[keep2, keep2]
    end

    # Normalize to cost shares: Ω[i,j] = IO[i,j] / sum_j(IO[i,j])
    row_sums = sum(io_reduced, dims=2)
    Ω = io_reduced ./ row_sums
    # Handle any NaN from zero rows
    Ω[isnan.(Ω)] .= 0.0

    # Factor shares for the base year
    α = vadd[:, yr_idx] ./ grossy[:, yr_idx]

    # Consumption shares: β = grossy' · (I - diag(1-α)·Ω)
    # In MATLAB: beta = grossy(:,year-1959)' * (eye(N) - diag(1-alpha)*Omega)
    β = grossy[:, yr_idx]' * (I - diagm(0 => 1 .- α) * Ω)
    β = vec(β)
    β[β .< 0] .= 0.0  # remove negative implied final sales
    β = β ./ sum(β)   # normalize

    # Domar weights: λ = (I - diag(1-α)·Ω)⁻¹' · β
    λ = (I - diagm(0 => 1 .- α) * Ω)' \ β

    # Baseline labor allocation: L = λ ⊙ α
    L = λ .* α

    return BFData(Ω, α, β, λ, L, N, year, grossy, vadd, capital, labor)
end

"""
    load_tfp_data(csvpath) -> (stfp::Matrix{Float64}, Sigma::Matrix{Float64}, mu::Vector{Float64})

Load sectoral TFP growth rates from stfp.csv.
stfp.csv: 76 rows (sectors), each row is a time series of TFP growth rates.
"""
function load_tfp_data(csvpath::String)
    raw = readdlm(csvpath, ',')
    stfp = Float64.(raw)

    # MATLAB: stfp(:,1) = [] (remove first year, which is blank)
    stfp = stfp[:, 2:end]
    # MATLAB: stfp(end,:) = [] (remove private household sector)
    stfp = stfp[1:end-1, :]

    Sigma = cov(stfp, dims=2)
    mu = vec(mean(stfp, dims=2))

    return stfp, Sigma, mu
end

"""
    describe_data(data::BFData)

Print summary statistics of the loaded data.
"""
function describe_data(data::BFData)
    @printf("B&F Replication Data Summary\n")
    @printf("============================\n")
    @printf("Sectors: %d\n", data.N)
    @printf("Base year: %d\n", data.year)
    @printf("\n")
    @printf("Ω matrix: %dx%d (cost shares)\n", size(data.Ω)...)
    @printf("α (factor shares): mean=%.4f, min=%.4f, max=%.4f\n",
        mean(data.α), minimum(data.α), maximum(data.α))
    @printf("β (consumption shares): sum=%.4f, max sector=%.4f\n",
        sum(data.β), maximum(data.β))
    @printf("λ (Domar weights): sum=%.4f, max=%.4f\n",
        sum(data.λ), maximum(data.λ))
    @printf("L (labor allocation): sum=%.4f\n", sum(data.L))
    @printf("\n")
    @printf("Key checks:\n")
    @printf("  λ ⊙ α ≈ L? max diff = %.2e\n", maximum(abs.(data.λ .* data.α .- data.L)))
    @printf("  (I - diag(1-α)·Ω)⁻¹' · β ≈ λ? max diff = %.2e\n",
        maximum(abs.((I - diagm(0 => 1 .- data.α) * data.Ω)' \ data.β .- data.λ)))
    @printf("  β sums to 1? sum = %.10f\n", sum(data.β))
end

end # module
