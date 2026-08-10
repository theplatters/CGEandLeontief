# src/io_table.jl — Load and calibrate the BEA 66-sector IO table
#
# Ports the calibration block of Master_file_3.m (lines 59–142).
# Reads IO_data_2018.mat and returns the calibrated matrices.
# Works with year-specific 2D slices for clarity and correctness.

using MAT
using LinearAlgebra
using Printf

export load_io_table, IOData

"""
    IOData

Holds all calibrated IO table quantities from the BEA 66-sector data.

# Fields
- `N::Int`                   — number of sectors (66)
- `Omega::Matrix{Float64}`   — row-normalized input-output coefficients (N×N)
- `alphaL::Vector{Float64}`  — labor share of value added (N)
- `alphaK::Vector{Float64}`  — capital share of value added (N)
- `beta::Vector{Float64}`    — final consumption shares (N, normalized)
- `int_share::Vector{Float64}`  — intermediate share of gross output (N)
- `va_share::Vector{Float64}`   — value-added share of gross output (N)
- `grossout::Vector{Float64}`   — gross output by sector (N)
- `labor::Vector{Float64}`      — labor compensation by sector (N)
- `gos::Vector{Float64}`        — gross operating surplus / capital (N)
- `va_sales::Vector{Float64}`   — value-added = labor + capital (N)
- `indname::Vector{String}`     — sector names (N or full)
- `year::Int`                   — base year used (2015)
"""
struct IOData
    N::Int
    Omega::Matrix{Float64}
    alphaL::Vector{Float64}
    alphaK::Vector{Float64}
    beta::Vector{Float64}
    int_share::Vector{Float64}
    va_share::Vector{Float64}
    grossout::Vector{Float64}
    labor::Vector{Float64}
    gos::Vector{Float64}
    va_sales::Vector{Float64}
    indname::Vector{String}
    year::Int
end

"""
    load_io_table(data_path::AbstractString; N=66, year=2015) -> IOData

Load and calibrate the BEA IO table from `IO_data_2018.mat`.

The .mat file contains `Data_raw` (Nb×Ncols×T) where T = 22 years (1997–2018).
Column map (per Master_file_3.m):
- Tot_use = col 94
- Tot_int = col 74
- Pce     = col 75

With year = 2015, the labor/capital rows are 77/79 (pre-2016 layout).
"""
function load_io_table(data_path::AbstractString; N::Int=66, year::Int=2015)
    @printf("Loading IO table from %s ... ", data_path)
    mat = matread(data_path)
    Data_raw = mat["Data_raw"]  # 3-D array: [row, col, year_index]
    raw_indname = mat["indname"]

    # year_index = year − 1996  (years 1997..2018 → indices 1..22)
    yi = year - 1996

    Ti = size(Data_raw, 3)
    if yi < 1 || yi > Ti
        error("year index $yi out of range (1..$Ti) for year=$year")
    end

    # Convert indname to Vector{String}
    if eltype(raw_indname) <: AbstractString
        indname = collect(raw_indname)
    else
        indname = [string(raw_indname[i]) for i in 1:length(raw_indname)]
    end

    # ---- Calibration (port of MATLAB lines 59–99, year < 2016 branch) ----
    # MATLAB: IO = Data_raw;  (3D)
    # Then: IO(:,1:2,:) = []; — removes columns 1 and 2 from all layers
    # After this, the 3D array has (rows, cols-2, years)
    # grossout = sum(IO(1:N,1:N,:)) + IO(77,1:N,:) + IO(79,1:N,:)
    # This is a 3D operation. We replicate it column-wise, then extract the year.

    if year < 2016
        labor_row = 77
        gos_row   = 79
    else
        labor_row = 75
        gos_row   = 76
    end

    # Extract the full 3D slice to replicate the column-removal step
    # Data_raw with columns 1:2 removed (same as MATLAB IO(:,1:2,:) = [])
    Clone_red = Data_raw[:, 3:end, :]  # remove first two columns entirely

    # Gross output: sum(intermediate block) + labor row + GOS row → per year
    # MATLAB: sum(IO(1:N,1:N,:)) over dim 1 → 1×N×T → squeeze → N×T
    # First extract the core IO block and clean it (MATLAB: IO(isnan(IO))=0; IO(IO<0)=0)
    core_block = Clone_red[1:N, 1:N, :]  # N×N×T
    for t in 1:size(core_block, 3)
        slice = @view core_block[:, :, t]
        replace!(x -> isnan(x) ? 0.0 : x, slice)
        replace!(x -> x < 0 ? 0.0 : x, slice)
    end
    int_sum = dropdims(sum(core_block, dims=1); dims=1)  # N×T

    # labor row and GOS row
    labor_3d = Clone_red[labor_row, 1:N, :]
    gos_3d = Clone_red[gos_row, 1:N, :]

    # grossout = intermediate + labor + GOS
    grossout_all = int_sum .+ labor_3d .+ gos_3d  # N×T

    # Specific year's IO block (already cleaned via core_block above)
    temp = core_block[:, :, yi]'  # N×N, transposed per MATLAB's convention
    row_sums = sum(temp, dims=2)
    Omega = temp ./ row_sums

    # Labor and capital for specific year
    L = vec(labor_3d[:, yi])      # N-vector
    K = vec(gos_3d[:, yi])        # N-vector
    va_sales = L .+ K
    alphaL = L ./ va_sales
    alphaK = K ./ va_sales

    # Gross output
    grossout_vec = vec(grossout_all[:, yi])

    # Final consumption (branch-dependent per Master_file_3.m)
    if year < 2016
        # Pre-2016: Final is directly stored in Data_raw column 99
        # Final = IO(:,99,:) before column removal
        Final = vec(Data_raw[1:N, 99, yi])
    else
        # Year ≥ 2016: Final = Tot_use - Tot_int
        # Tot_use = Data_raw(1:N, 94, yi), Tot_int = Data_raw(1:N, 74, yi)
        Tot_use = Data_raw[1:N, 94, yi]
        Tot_int = Data_raw[1:N, 74, yi]
        Final = Tot_use .- Tot_int
    end

    # Beta: clean and normalize
    beta = copy(Final)
    replace!(x -> isnan(x) ? 0.0 : x, beta)
    replace!(x -> x < 0 ? 0.0 : x, beta)
    beta_sum = sum(beta)
    if beta_sum > 0
        beta ./= beta_sum
    end

    # Intermediate and value-added shares (use cleaned core_block)
    int_sales = vec(sum(core_block[:, :, yi], dims=1))  # row sums of Ω block
    int_share_vec = int_sales ./ (int_sales .+ va_sales)
    va_share_vec = va_sales ./ (int_sales .+ va_sales)

    @printf("done. N=%d, year=%d, beta sum=%.6f\n", N, year, sum(beta))

    return IOData(
        N,
        Omega,
        alphaL, alphaK, beta,
        vec(int_share_vec), vec(va_share_vec),
        grossout_vec, L, K, va_sales,
        indname, year
    )
end