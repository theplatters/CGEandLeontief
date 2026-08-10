# src/shocks.jl — Load COVID-19 shock data from .xlsx files
#
# Ports the shock-loading section of Master_file_3.m (lines 22–55).

using XLSX
using LinearAlgebra
using Printf

export load_shocks, ShockData

"""
    ShockData

Holds the four shock/outcome vectors read from pre-built .xlsx files.

# Fields
- `N::Int`                       — number of sectors
- `BLS_shock::Vector{Float64}`   — supply shock (hours change, col 4 diff_2005)
- `PCE_shock::Vector{Float64}`   — demand shock (PCE change, col 4 diff_2005_pce)
- `wages::Vector{Float64}`       — wage change (col 6 w_adj_20Q1_20Q2)
- `PPI::Vector{Float64}`         — PPI Feb→May change (col 5 av_p_change_Feb_May)
"""
struct ShockData
    N::Int
    BLS_shock::Vector{Float64}
    PCE_shock::Vector{Float64}
    wages::Vector{Float64}
    PPI::Vector{Float64}
end

"""
    read_xlsx_column(filepath::String, col::Int, nrows::Int; start_row::Int=1) -> Vector{Float64}

Read a single column of numeric data from an .xlsx file using cell-access API.
"""
function read_xlsx_column(filepath::String, col::Int, nrows::Int; start_row::Int=1)
    data = Vector{Float64}(undef, nrows)
    XLSX.openxlsx(filepath) do xf
        sheet = xf[1]  # first sheet
        for i in 1:nrows
            r = start_row + i - 1
            cell_val = sheet[r, col]
            if cell_val isa Number
                data[i] = Float64(cell_val)
            elseif cell_val isa AbstractString
                parsed = tryparse(Float64, cell_val)
                data[i] = parsed === nothing ? 0.0 : parsed
            else
                data[i] = 0.0
            end
        end
    end
    return data
end

"""
    load_shocks(data_dir::AbstractString; N::Int=66) -> ShockData

Load the four shock/outcome .xlsx files and return structured data.

File mapping (from Master_file_3.m):
- BLS_labor_shock_202108.xlsx : col 4 = diff_2005  (supply shock)
- PCE_shock_202107.xlsx       : col 4 = diff_2005_pce (demand shock)
- wage_change_final.xlsx      : col 6 = w_adj_20Q1_20Q2
- ppi_data.xlsx               : col 5 = av_p_change_Feb_May
"""
function load_shocks(data_dir::AbstractString; N::Int=66)
    @printf("Loading shocks from %s ...\n", data_dir)

    bls_path = joinpath(data_dir, "BLS_labor_shock_202108.xlsx")
    pce_path = joinpath(data_dir, "PCE_shock_202107.xlsx")
    wage_path = joinpath(data_dir, "wage_change_final.xlsx")
    ppi_path  = joinpath(data_dir, "ppi_data.xlsx")

    BLS_shock_raw = read_xlsx_column(bls_path, 4, N; start_row=2)
    PCE_shock_raw = read_xlsx_column(pce_path, 4, N; start_row=2)
    wages_raw     = read_xlsx_column(wage_path, 6, N; start_row=2)
    PPI_data_raw  = read_xlsx_column(ppi_path,  5, N; start_row=2)

    # HS/ORE merge: row 49 (HS, Housing) ← row 48 (ORE, Other real estate)
    # Both map to NAICS 531 (Real estate). MATLAB: PPI_data_raw(49,3:end) = PPI_data_raw(48,3:end)
    if length(PPI_data_raw) >= 49
        PPI_data_raw[49] = PPI_data_raw[48]
    end

    @printf("  BLS shock: min=%.4f, max=%.4f\n", minimum(BLS_shock_raw), maximum(BLS_shock_raw))
    @printf("  PCE shock: min=%.4f, max=%.4f\n", minimum(PCE_shock_raw), maximum(PCE_shock_raw))

    return ShockData(N, BLS_shock_raw, PCE_shock_raw, wages_raw, PPI_data_raw)
end