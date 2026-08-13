# src/network.jl — Build the standard-form production network Ω_re
#
# Ports Master_file_3.m lines 144–211: relabels the IO table into a
# "standard form" with N goods, N value-added, N intermediates, N labor,
# N capital, and three auxiliary agents (HtM, Ricardian, tomorrow's consumption).

using LinearAlgebra
using Printf

export build_standard_form, StandardForm

"""
    StandardForm

Holds the standard-form production network matrices and categorical vectors.

# Fields
- `N::Int`              — number of sectors
- `D::Int`              — total dimension = 5N + 4
- `Omega_re::Matrix`    — standard-form input-output matrix (D×D)
- `Psi_re::Matrix`      — Leontief inverse (I − Ω_re)⁻¹
- `factor::Vector{Int}` — type codes: 1=good, 0=factor, 2=Ricardian, 3=HtM
- `keynes::Vector{Int}` — price rigidity: 1=normal, 0=flexible, -1=sticky
- `chi::Vector{Float64}` — factor ownership indicator
- `phi::Vector{Float64}`  — factor supply elasticity (all 0)
- `Domar::Vector{Float64}` — Domar weights (Psi_re[1, 2:N+1])
"""
struct StandardForm
    N::Int
    D::Int
    Omega_re::Matrix{Float64}
    Psi_re::Matrix{Float64}
    factor::Vector{Int}
    keynes::Vector{Int}
    chi::Vector{Float64}
    phi::Vector{Float64}
    Domar::Vector{Float64}
end

"""
    build_standard_form(io::IOData) -> StandardForm

Build the standard-form Ω_re from calibrated IO data.

Block layout (1-indexed):
   1           : consumption today
   2..N+1      : N goods today
   N+2..2N+1   : N value-added
   2N+2..3N+1  : N intermediates
   3N+2..4N+1  : N labor (sticky, keynes=-1)
   4N+2..5N+1  : N capital (flexible, keynes=0)
   5N+2        : HtM consumer (factor=3)
   5N+3        : Ricardian consumer (factor=2)
   5N+4        : consumption good tomorrow (factor=0, numeraire)

D = 5N + 4
"""
function build_standard_form(io::IOData)
    N = io.N
    D = 5 * N + 4

    @printf("Building standard form Ω_re: N=%d, D=%d\n", N, D)

    Omega_re = zeros(Float64, D, D)

    # Row 1: consumption buys N goods with shares beta
    Omega_re[1, 2:(N+1)] = io.beta

    # Rows 2:N+1 (goods): use VA (N+2:2N+1) and intermediates (2N+2:3N+1)
    mu = 1.0  # markup per Master_file_3.m
    for i in 1:N
        Omega_re[1+i, N+1+i] = io.va_share[i] / mu     # goods → VA (diagonal)
        Omega_re[1+i, 2*N+1+i] = io.int_share[i] / mu  # goods → intermediates (diagonal)
    end

    # Rows N+2:2N+1 (value-added): use labor (3N+2:4N+1) and capital (4N+2:5N+1)
    for i in 1:N
        total_factor = io.alphaL[i] + io.alphaK[i]
        Omega_re[N+1+i, 3*N+1+i] = io.alphaL[i] / total_factor  # VA → labor
        Omega_re[N+1+i, 4*N+1+i] = io.alphaK[i] / total_factor  # VA → capital
    end

    # Rows 2N+2:3N+1 (intermediates): use goods (2:N+1) via Omega matrix
    Omega_re[(2*N+2):(3*N+1), 2:(N+1)] = io.Omega

    # Row 5N+2 (HtM consumer, index 332): consumes tomorrow's good (col 334)
    Omega_re[5*N+2, D] = 1.0

    # Row 5N+3 (Ricardian consumer, index 333): consumes today (col 1) and tomorrow (col 334)
    Omega_re[5*N+3, 1]  = 0.5
    Omega_re[5*N+3, D]  = 0.5

    # Leontief inverse: Psi_re = (I - Ω_re)⁻¹
    Psi_re = inv(I - Omega_re)

    # Factor type codes (categories for the AMPL model)
    factor = ones(Int, D)             # 1 = good
    factor[(3*N+2):(5*N+1)] .= 0     # 0 = factor (labor and capital)
    factor[D] = 0                     # consumption tomorrow = factor
    factor[D-2] = 3                   # HtM consumer
    factor[D-1] = 2                   # Ricardian consumer

    # Keynes rigidity codes
    keynes = ones(Int, D)             # 1 = normal good
    keynes[(4*N+2):(5*N+1)] .= 0     # capital = flexible
    keynes[D] = 0                     # tomorrow's consumption = flexible
    keynes[(3*N+2):(4*N+1)] .= -1    # labor = sticky

    # Chi: ownership of factors.
    # MATLAB: chi = (factor==0).*rand(D,1), then labor/capital/tomorrow rows are
    # explicitly zeroed → chi ≡ 0 deterministically. We set it directly.
    chi = zeros(Float64, D)
    # (All factor rows are zeroed in MATLAB, so chi is identically zero.)

    # Phi: factor supply elasticity (all 0 in this model)
    phi = zeros(Float64, D)

    # Domar weights = Psi_re[1, 2:N+1] (GDP share of each sector)
    Domar = Psi_re[1, 2:(N+1)]

    @printf("  Psi_re[1,2:N+1] sum = %.6f\n", sum(Domar))

    return StandardForm(N, D, Omega_re, Psi_re, factor, keynes, chi, phi, Domar)
end