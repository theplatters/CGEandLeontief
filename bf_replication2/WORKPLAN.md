---
title: "Replication Workplan — Baqaee & Farhi (2022) COVID-19 Application in Julia"
project: BFRep / RepAEA2022 → bf_replication2
date: 2026-08-09
revised: 2026-08-09 (corrected — all data inputs confirmed present; Route B solver confirmed)
tags: [workplan, replication, julia, keynesian, production-network, covid, baqaee-farhi, bfrep]
---

# Replication Workplan — Baqaee & Farhi (2022) COVID-19 Application in Julia

## Revision Note

An earlier draft of this plan wrongly concluded that `IO_data_2018.mat`,
`BLS_labor_shock_202108.xls`, and `PCE_shock_202107.xls` were missing. That was
an artifact of a malfunctioning directory search — a direct `find` of
`RepAEA2022/` shows **every input the MATLAB driver needs is present** in
`Replication code_ver2/`, including all the raw microdata. There is no
missing-data blocker. The three `.xls` files have been manually converted to
`.xlsx` and confirmed present in the same folder. This revision also switches
the solver strategy from the initially-proposed custom Newton+Fischer–Burmeister
(Route A) to **JuMP + PATHSolver.jl (Route B)**, per user direction.

## Overview

The `RepAEA2022` folder contains a multi-language replication of Baqaee & Farhi
(2022), *Supply and Demand in Disaggregated Keynesian Economies with an
Application to the COVID-19 Crisis* (AEA Papers and Proceedings). The
"application" calibrates a **static, disaggregated Keynesian production-network
model** to the US economy (66 BEA sectors, 2018) and feeds it the 2020 COVID
shocks to decompose the recession into supply- vs demand-driven components.

The model is solved in **AMPL + KNITRO** as a system of nonlinear equations with
a **Calvo-style sticky-price complementarity** (the Keynesian core) and
**Hand-to-Mouth (HtM) vs Ricardian** consumers. This workplan ports the whole
pipeline to Julia, following the clean-room approach already used for the 2019
paper in `bf_replication/` (see the `beyond-hulten-julia-modeling` skill).

The four source code files:

- `Master_file_1.R` — builds the sectoral **labor shock** from BLS employment tables.
- `Master_file_2.do` — Stata: builds **BLS labor**, **PCE demand**, **wage**, and **PPI** shock files from raw BLS/BEA/CPS microdata.
- `Master_file_3.m` — MATLAB driver: loads IO table + shocks, builds the **standard-form Ω**, calls AMPL/KNITRO across an elasticity × shock-type × HtM grid, and produces all figures.
- `Standard_form_Covid_HtM3_logCES.mod` — the AMPL equilibrium model (prices, inverse demands λ, factor clearing, numeraire).

## Source Material Inventory

All files confirmed present in `RepAEA2022/Replication code_ver2/`.

Code files:

| File | Role |
|------|------|
| `Master_file_1.R` | Labor-shock construction (BLS → EPI 70-code) |
| `Master_file_2.do` | Stata: BLS / PCE / wage / PPI shocks |
| `Master_file_3.m` | MATLAB driver + AMPL bridge + figures |
| `Standard_form_Covid_HtM3_logCES.mod` | AMPL equilibrium system |
| `Readme.pdf` | data README (read on host — no PDF tool in container) |

Data files consumed directly by `Master_file_3.m` (the model's essential inputs):

| File | Feeds | Format |
|------|-------|--------|
| `IO_data_2018.mat` | BEA 66-sector IO table (`Data_raw` 3-D: sectors × variables × years 1997–2018; `indname`) | MATLAB 5.0 `.mat` |
| `BLS_labor_shock_202108.xlsx` | sectoral supply shock `A` (col 4 = `diff_2005`) | OOXML `.xlsx` |
| `PCE_shock_202107.xlsx` | sectoral demand shock `B` (col 4 = `diff_2005_pce`) | OOXML `.xlsx` |
| `wage_change_final.xlsx` | wage changes (col 6 = `w_adj_20Q1_20Q2`) | OOXML `.xlsx` |
| `ppi_data.xlsx` | PPI % change Feb→May (col 5 = `av_p_change_Feb_May`) | OOXML `.xlsx` |

All files are now in a Julia-readable format: `XLSX.jl` reads the `.xlsx` files
directly; the `.mat` is MATLAB 5.0 format (readable by `MAT.jl`).

Raw inputs used by the R/Stata prep scripts (all present — only needed if full
end-to-end provenance is desired):

| File | Used by |
|------|---------|
| `Employment_TableB1.xlsx`, `Xwalk.xlsx` | `Master_file_1.R` |
| `Expenditure_202107.xlsx`, `bridge_2018.dta`, `IO_to_NAICS(1).xlsx`, `2017-industry-code-list-with-crosswalk.xlsx`, `Missing-industries.dta`, `industry_sales_census.dta/_new.xlsx/_org.xlsx`, `series_id.dta/.xlsx/_ppi.txt`, `ppi_bls.xlsx` | `Master_file_2.do` |
| `cps_update.dta` (193 MB), `QCEW_2019Q4/2020Q1/2020Q2.csv`, `iocode_skeleton.dta`, `LaborShock_EPImodel.csv` | `Master_file_2.do` |

## Model Specification (reverse-engineered)

The `.mod` file is a **standard-form relabeling** of the production network.
With `N` sectors and `F = 2N` factors, the dimension is `D = 1 + 3N + F + 3 = 5N + 4`.
For `N = 66` this gives `D = 334` (confirmed: the `.mod` references `member(334, Ind)`).

```julia
# Block layout of the standard-form matrix Omega_re (D x D)
# 1              : consumption today
# 2   .. N+1     : N goods today
# N+2 .. 2N+1    : N value-added today
# 2N+2 .. 3N+1   : N intermediates today
# 3N+2 .. 4N+1   : N labor (sticky factor, keynes = -1)
# 4N+2 .. 5N+1   : N capital (flexible factor, keynes = 0)
# 5N+2           : HtM consumer       (factor code 3)
# 5N+3           : Ricardian consumer (factor code 2)
# 5N+4           : consumption good tomorrow (numeraire)
```

The equilibrium system (port of `Standard_form_Covid_HtM3_logCES.mod`):

- **Price equations** — CES (`cobb_douglas=0`): `p_k = (μ/in_μ)/A_k · (Σ_j B_kj^θ Ω_kj p_j^(1-θ))^(1/(1-θ))`; or Cobb-Douglas (`cobb_douglas=1`): `log p_k = −log A_v + Σ_j B_kj Ω_kj log p_j`.
- **Consumption share** — for goods: `cons_share_k = B_1k^θ Ω_1k (p_k/p_1)^(1-θ) A_1^(θ−1)`.
- **Inverse-demand λ equations** — for goods/factors (`factor<2`): cost-push propagation of λ; for the two consumers (`factor=3,2`): HtM vs Ricardian budget constraints with the `φ_htm` social-insurance parameter.
- **Factor clearing** — capital & tomorrow's consumption (`keynes=0`): equality `λ = p·(p/p_1)^φ·λ̄·A`. Labor (`keynes=-1`): a **complementarity** `(p − ((λ/p)/(A λ̄))^φ)·(λ/p − A λ̄) = 0` with `p ≥ ((λ/p)/(A λ̄))^φ` and `λ/p ≤ A λ̄`. This splits sectors into **demand-constrained** (price stuck above flexible level → unemployment) vs **supply-constrained** (price falls to clear the market).
- **Numeraire** — `p[end] = 1` (consumption tomorrow is the price anchor).

Shocks are applied through `A` (labor supply rows `3N+2..4N+1`) and `B` (row 1
sectoral demand, row `5N+4` aggregate demand), scaled by a grid `t ∈ [0,1]`.

Calibration (from `Master_file_3.m`):

| Regime | θ (intermediate) | ε (VA↔intermediate) | σ (consumption) | η (across VA) | note |
|--------|------------------|---------------------|-----------------|---------------|------|
| Benchmark (complementarity) | 0.2 | 0.6 | 1.0 | 0.6 | Figure 2 |
| Cobb-Douglas | 1.0 | 1.0 | 1.0 | 1.0 | Figure 3 |

Five shock types: (1) baseline = supply + sectoral demand + aggregate demand;
(2) supply only; (3) sectoral + aggregate demand; (4) aggregate demand only;
(5) supply + sectoral demand. HtM share `φ_htm` swept `0, 0.2, …, 1.0` (Figure 4).

## Replication Strategy

Replicate `Master_file_3.m`'s outputs (Figures 2–4, Appendix A1–A3 / Tables A1–A2,
and the out-of-sample fit in Parts B–D) in Julia.

Solver route — **Route B (JuMP + PATHSolver.jl)**:

- Primary: express the entire `.mod` equilibrium as a **Mixed Complementarity
  Problem (MCP)** in JuMP and solve with PATHSolver.jl. This is the most faithful
  approach — PATH is the same solver family as KNITRO (both handle MCP/complementarity
  natively) and maps directly from the `.mod`'s complementarity structure.
- Fallback: if PATHSolver.jl precompilation causes container issues, either (a)
  switch to JuMP + Ipopt (NLP formulation with complementarity constraints), or
  (b) drop back to Route A (custom Newton + Fischer–Burmeister, stdlib-only).

Data route — all inputs are present in clean Julia-readable formats.

- **Fast path (recommended):** consume the five essential model inputs directly
  — `IO_data_2018.mat`, `BLS_labor_shock_202108.xlsx`, `PCE_shock_202107.xlsx`,
  `wage_change_final.xlsx`, `ppi_data.xlsx`. No need to re-run the R/Stata prep
  unless full end-to-end provenance is required.
- `.mat` reading: load with `MAT.jl` (344 KB, non-HDF5 format — straightforward).
  If `MAT.jl` precompilation hangs in the container, pre-export
  `Data_raw` + `indname` to Arrow/CSV on the host Mac.
- `.xlsx` reading: `XLSX.jl` reads all five files directly.
- Reuse the existing `bf_replication/` (2019) code patterns for IO-construction
  logic where applicable, and apply the `matlab-comparison-protocol`
  verification method: equation-by-equation against the `.mod` source, then
  against the paper's reported figures.

## Scaffold and Julia Environment

- Create `bf_replication2/` with `src/` (data, network, model, solver, driver),
  `data/` (symlink or copy the present inputs), `figures/`, and a `Project.toml`.
- Minimal `Project.toml`:

```toml
name = "BFAEA2022"
uuid = "5a3a2a4e-1b0e-4a1a-9b2d-2e1c3a4b5c6d"

[deps]
JuMP = "4076af6c-e467-56ae-b986-b466b2749572"
PATHSolver = "024fb59e-2b6b-5e2b-9e4d-8f0a6f7e18b0"
XLSX = "fd0916cc-8b04-5870-babc-2e136c64c1e6"
MAT = "23992714-dd62-5051-b70f-ba57cb901cac"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
```

- Run `Pkg.instantiate()` and monitor for OOM. If the container OOMs during
  precompile, kill from the host, switch to a lighter set of deps (drop
  `MAT.jl`, pre-export data on host; drop `PATHSolver`, fall back to Route A).
- Symlink/verify Julia 1.12.6 at `/usr/local/bin/julia`. Confirm
  `julia -e 'println("OK")'` returns instantly.
- Add a `README.md` documenting how to drop the data files into `data/`.

## Data Layer

- **IO table** (`src/io_table.jl`): load `IO_data_2018.mat` with `MAT.jl`.
  `Data_raw` is indexed `Data_raw[row, column, year_index]` with
  `year_index = year − 1996` (years 1997–2018) and the column map used by
  `Master_file_3.m`: `Final = col 99`, `Tot_use = col 94`, `Tot_int = col 74`,
  `Pce = col 75`; labor/capital rows are 75/76 (year ≥ 2016 branch). Reproduce
  the calibration exactly: `grossout`, `labor`, `gos`, `Omega` (row-normalized),
  `alphaL/alphaK`, `beta` (final, normalized, negatives zeroed), `int_share`,
  `va_share`. Keep `indname` for sector labels.
- **Shocks** (`src/shocks.jl`): read `BLS_labor_shock_202108.xlsx` (col 4) → `A`
  supply shock; `PCE_shock_202107.xlsx` (col 4) → `B` sectoral demand;
  `wage_change_final.xlsx` (col 6) → wages; `ppi_data.xlsx` (col 5) → PPI.
  Apply the HS/ORE merge (`PPI_data_raw[49,:] ← PPI_data_raw[48,:]`, both NAICS 531).
- **Standard-form Ω** (`src/network.jl`): build `Omega_re` (`D × D`) and
  `Psi_re = (I − Ω_re)⁻¹`, plus the `factor`, `keynes`, `chi`, `phi`,
  `phi_htm` categorical vectors, exactly per `Master_file_3.m` lines 144–211.
- Unit tests: `Psi_re[1, 2:N+1]` (Domar weights) sum to 1; `Omega_re` block
  stochasticity; `D = 334` for `N = 66`.

## Equilibrium System and Solver

- `src/model.jl`: express every `.mod` constraint in JuMP:

  - **Variables**: `p[1:D] ≥ 0`, `λ[1:D] ≥ 0`
  - **Price equations** — CES or CD (branch on `cobb_douglas` vector):
    `@NLconstraint(model, p[k] == ...)` for CES; `@constraint(model, log(p[k]) == ...)` for CD.
  - **Consumption share** — computed algebraically from `p`.
  - **Lambda equations** — linear propagation of λ through Ω; the HtM/Ricardian
    lambda blocks become linear equations.
  - **Factor clearing (flexible)** — equality for `keynes == 0` (capital).
  - **Factor clearing (sticky)** — for `keynes == -1` (labor): use JuMP's
    complementarity support (`@complementarity`) or reformulate as a
    complementarity constraint via `@constraint` with a slack. PATHSolver
    handles MCP natively, so the three conditions (one equality with product = 0
    plus two inequalities) become a mixed complementarity problem.
  - **Numeraire** — `p[end] == 1`.

- Expose `solve_equilibrium(Omega_re, shocks, theta, phi_htm,
  cobb_douglas) → (p, λ)` that builds the JuMP model, attaches `PATHSolver`,
  and returns the solution.
- Validation gate: on the real `N=66` Ω with `t=0`, confirm `p[end] == 1`,
  all equations hold, and the complementarity residual is near zero.

## Driver Loop and Calibration Grid

Port the four nested loops of `Master_file_3.m` (loop 1/2 × elasticity 1/2 ×
shock_type × HtM share `s` × `t` grid). For each `(elasticity, shock_type, s,
loop)` store: real GDP `λ[1]/p[1]`, nominal GDP `λ[1]`, Hulten benchmark
`exp(Psi_re[1,:]·log A)`, prices, consumption shares, unemployment measures
(`unemployment1/2/3`), inflation `p[1]`, and the `Trunc_A` constraint tracker.
Reproduce the `B` row renormalization (`B[1,:] / sum(Ω_re[1,:]·B[1,:]')`, same
for row `5N+3`).

## Outputs, Figures, and Out-of-Sample Fit

Reproduce the `disp()` blocks and figures:

- **Part A** — calibration: targeted nominal-GDP reduction.
- **Part B** — out-of-sample fit: model-implied vs data PPI inflation in
  demand-constrained (sticky) vs supply-constrained (flex) sectors; wage
  inflation and worked-hours split by constraint type.
- **Part C / Figure 2** — benchmark bars (Real GDP, Inflation, Keynesian
  Unemployment) for Baseline / Demand / Supply; numeric reductions.
- **Part D** — tightness/slackness: count of demand vs supply sectors;
  model-implied unemployment for `supply_con=[19,65,7,41]`,
  `demand_con=[32,34,33,24,3]`, `healthcare=[59,58]`.
- **Part E / Figure 3** — Cobb-Douglas analogs.
- **Part F / Figure 4** — policy: aggregate-demand effect, HtM with/without,
  and the HtM-share sweep (RGDP, Inflation, Unemployment vs `φ_htm ∈ 0..1`).
- **Appendix** — Fig A1 (PCE & employment histograms), A2 (model vs data PPI
  scatter, 45° line), A3 (sectoral hours bar); Table A1/A2 (price growth,
  model vs data, by constraint type).

Use `Plots.jl`/`GLMakie` only on the host (they OOM the container); for the
container, emit CSVs of every series and generate figures on the Mac.

## Validation

Primary target is the MATLAB driver's own `disp()` outputs (run the original in
MATLAB on the host if available, else compare structure + the paper's reported
numbers). Secondary target is the published paper's Figures 2–4 and Appendix
tables. Apply the two-tier check from the `matlab-comparison-protocol`:
equation-by-equation against the `.mod` source first, then against the paper's
claims. Record a `VALIDATION.md` with side-by-side MATLAB-vs-Julia numbers.

## Documentation

- `WORKPLAN.md` (this file), `README.md` (setup + `.mat`/`.xlsx` sources),
  `VALIDATION.md` (replication diffs), and inline docstrings in `src/`.
- Note deviations explicitly (solver choice, any data approximation, sector
  ordering) so the replication is auditable.

## Risk Register

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|-----------|
| JuMP + PATHSolver.jl precompile OOMs container | Medium | High | Drop to JuMP + Ipopt (binary libraries); if still OOM, fall back to Route A (stdlib Newton+FB) |
| `MAT.jl` precompilation hangs | Low | Low | Pre-export `.mat` contents to Arrow/CSV on host (344 KB trivial) |
| Complementarity (sticky labor) mishandled by JuMP → PATH | Low | High | Test on a synthetic small Ω first; PATH natively solves MCP, so the `.mod` maps 1-to-1 |
| Container OOM on heavy precompile | High | Medium | Minimal `Project.toml`; run heavy operations on host |
| Sector/code crosswalk subtleties (HS/ORE, iocode ordering) | Medium | Medium | Mirror `Master_file_3.m` row ops exactly; unit-test `D=334` |
| N ambiguity (66 vs other sector aggregations) | Low | Medium | Fix `N=66`, `year=2018` per the driver; document |
| PDF not machine-readable here | Certain | Low | Read `Readme.pdf` / paper on host for exact magnitudes |

## Execution Order

| Order | Phase | Blocks | Effort |
|-------|-------|--------|--------|
| 1 | Scaffold & environment (`Project.toml`, instantiate) | — | 0.5 day |
| 2 | Data layer (load `.mat` + `.xlsx`; build Ω) | — | 2–3 days |
| 3 | Standard-form Ω construction + tests | Phase 2 | 1–2 days |
| 4 | Equilibrium system in JuMP + PATHSolver | Phase 3 | 3–5 days |
| 5 | Driver loop & calibration grid | Phase 4 | 2–3 days |
| 6 | Outputs, figures, out-of-sample fit | Phase 5 | 2–4 days |
| 7 | Validation vs MATLAB/paper | Phase 6 | 2–3 days |
| 8 | Documentation | — | 1 day |

The critical path is now Phase 4 (the JuMP + MCP model) — the data is fully
present and in Julia-readable format.

## Open Questions for Jakob

- **Data-prep scope:** consume the four already-built model inputs directly
  (fast path, recommended), or also re-implement `Master_file_1.R` /
  `Master_file_2.do` for full end-to-end provenance? The built `.xlsx` shock
  files are present and sufficient for the model replication.
- **First milestone scope:** full Figures 2–4 + Appendix, or benchmark
  (Figure 2) + out-of-sample fit first?
- **OOM contingency:** if JuMP/PATHSolver.jl precompilation kills the
  container, do you prefer (a) JuMP + Ipopt fallback, or (b) straight to Route
  A (stdlib Newton+FB, no packages at all)?