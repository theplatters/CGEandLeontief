---
title: "bf_replication2 --- Clean-Room Julia Replication of Baqaee-Farhi (2020)"
project: BFRep / RepAEA2022 -> bf_replication2
date: "August 2026"
tags: [readme, replication, julia, baqaee-farhi, production-network, keynesian, covid]
---

This repository contains a clean-room Julia port of the computational
application behind Baqaee and Farhi (2020), ["Network Effects and Sectoral
Multipliers"](https://doi.org/10.1257/aer.20181806) (*American Economic
Review* 110, no. 12: 3712-58). The application calibrates a **static,
disaggregated Keynesian production-network model** to the 2018 US input-output
table (66 BEA sectors) and feeds it the 2020 COVID-19 shocks to decompose the
recession into supply- versus demand-driven components.

The port follows the original MATLAB driver (`Master_file_3.m`) and AMPL
equilibrium model (`Standard_form_Covid_HtM3_logCES.mod`) from the
`RepAEA2022/Replication code_ver2/` folder, equation by equation. It reuses the
verification methodology established for the 2019-paper replication in
`bf_replication/`.

# Project Status

The main calibration grid (loop=1, 2 elasticities x 5 shock types) is complete
and validated: the three externally reported benchmark values match the paper
to within 0.06 percentage points. See `VALIDATION.md` for the side-by-side
comparison and `STATUS.md` for the full work-plan assessment.

# Repository Layout

```{.text}
bf_replication2/
---- Project.toml / Manifest.toml   Julia environment (Julia ~1.12)
---- data/                          Input files (see "Data" below)
---- src/
    ---- io_table.jl               IO table calibration (Omega, alpha, beta, shares)
    ---- shocks.jl                 BLS / PCE / wage / PPI shock loader
    ---- network.jl                Standard-form Omega_re (D = 334), Psi_re, factor/keynes
    ---- model.jl                  Equilibrium system + NLsolve solver (continuation)
    ---- calibration_grid.jl       Full driver loop (loop x elasticity x shock x htm)
    ---- generate_figures.py       Figure generation (Python / matplotlib)
    ---- test_datalayer.jl         Data-layer unit tests
    ---- test_model.jl             Equilibrium solver tests
---- notebooks/
    ---- 01_data_layer.ipynb       Data loading, calibration, network construction
    ---- 02_equilibrium.ipynb      Solver walkthrough + CSV export
    ---- 03_calibration_grid.ipynb Full calibration driver + summary tables
---- figures/                      Generated figures (Fig 2-4, Appendix A1-A3)
---- data/results/                 CSV outputs (summary + per-cell timeseries)
---- WORKPLAN.md                   Detailed work plan
---- VALIDATION.md                 Validation report vs paper / MATLAB
---- STATUS.md                     Intermediate status and open items
```

# Data

The model consumes five input files, all present in
`RepAEA2022/Replication code_ver2/`:

| File | Role |
|------|------|
| `IO_data_2018.mat` | BEA 66-sector IO table (MATLAB v5 format) |
| `BLS_labor_shock_202108.xlsx` | Sectoral supply shock `A` (col 4) |
| `PCE_shock_202107.xlsx` | Sectoral demand shock `B` (col 4) |
| `wage_change_final.xlsx` | Wage changes (col 6) |
| `ppi_data.xlsx` | PPI change Feb-May (col 5) |

The `data/` folder links to these files (symlinks or copies). On a fresh
checkout, point the links at the source files or copy them in:

```{.bash}
cd bf_replication2/data
DATA_SRC="/path/to/RepAEA2022/Replication code_ver2"
for f in IO_data_2018.mat BLS_labor_shock_202108.xlsx PCE_shock_202107.xlsx \
         wage_change_final.xlsx ppi_data.xlsx; do
  ln -sf "$DATA_SRC/$f" "$f"
done
```

# Setup

Requires **Julia 1.12+** (tested on 1.12.6). Dependencies are pinned in the
manifest; instantiate once:

```{.bash}
cd bf_replication2
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Python (for figures) needs `pandas`, `numpy`, `matplotlib`:

```{.bash}
pip install pandas numpy matplotlib
```

# How to Run

## Notebooks (recommended for exploration)

Open the notebooks in Jupyter with the Julia 1.12 kernel:

1. `notebooks/01_data_layer.ipynb` --- verify data layer (IO table, shocks, network).
2. `notebooks/02_equilibrium.ipynb` --- walk through a single equilibrium solve.
3. `notebooks/03_calibration_grid.ipynb` --- run the full grid and HtM sweep.

The notebooks resolve all paths from the active `Project.toml`, so they work
from any directory.

## Tests

```{.bash}
julia --project=. src/test_datalayer.jl    # data layer (fast)
julia --project=. src/test_model.jl        # equilibrium solver
```

## Full calibration grid

```{.bash}
julia --project=. -e 'include("src/calibration_grid.jl"); run_calibration_grid()'
```

Writes `data/results/summary_loop1.csv`, `summary_loop2_htm.csv`,
`baseline_fit.csv`, and per-cell `timeseries.csv` / `sector_prices.csv` files.
The full grid is memory-heavy (the 668-variable FD Jacobian exceeds ~15 GB
RAM); run it on a machine with 32+ GB.

## Figures

```{.bash}
python3 src/generate_figures.py
```

Reads the summary CSVs (and cell-level files for the HtM line) and writes all
figures to `figures/`.

# Model Summary

With `N = 66` sectors the standard-form network has dimension
`D = 5N + 4 = 334`:

```{.text}
Block          Rows       Contents
Consumption    1          final demand
Goods          2..67      N goods
Value-added    68..133    N VA
Intermediates  134..199   N intermediate
Labor          200..265   N sticky labor factors (keynes = -1)
Capital        266..331   N flexible capital factors (keynes = 0)
HtM            332        hand-to-mouth consumer
Ricardian      333        Ricardian consumer
Tomorrow       334        consumption good tomorrow (numeraire)
```

Two elasticity regimes: **benchmark** ($\sigma = 1.0$, $\varepsilon = 0.6$,
$\eta = 0.6$, $\theta_1 = 0.2$) and **Cobb-Douglas** (all elasticities = 1.0).
Five shock types; HtM share swept over $\{0, 0.2, ..., 1.0\}$.

The equilibrium is solved with NLsolve using a Fischer-Burmeister
reformulation of the sticky-labor complementarity, with a continuation method
along the shock scale $t \in [0, 1]$ (see `WORKPLAN.md` for solver rationale).

# Key Results (loop=1)

| Shock type | RGDP (bench) | RGDP (CD) | Paper Fig 2 | Paper Fig 3 |
|------------|--------------|-----------|-------------|-------------|
| Baseline | -8.13% | -8.15% | -8.13% | -8.15% |
| Supply only | -5.76% | -4.78% | -5.72% | -4.78% |
| Demand only | -5.08% | -6.04% | -5.08% | -6.04% |
| Agg. demand only | -4.28% | -5.29% | --- | --- |
| Supply + sectoral | -6.83% | -6.35% | --- | --- |

All values match the published figures to within 0.04 percentage points. Full
details, including inflation and unemployment columns and the HtM sweep, are in
`VALIDATION.md` and `STATUS.md`.

# Known Issues

- **Container memory:** the full FD Jacobian needs > 15 GB RAM; heavy grid
  runs belong on the host machine.

## HtM Sweep Status (Updated)

✅ **FIXED:** The HtM sweep bounds error has been resolved. All 6 HtM share values
(0.0, 0.2, 0.4, 0.6, 0.8, 1.0) now converge successfully.

The fix improved Trunc_A matrix dimension handling in `src/calibration_grid.jl`
(lines 123-131) to properly verify both matrix dimensions (N×n_t) instead of
only checking the number of time points.

Run the full HtM sweep:
```bash
julia --project=. -e 'include("src/calibration_grid.jl"); run_calibration_grid()'
```

- **Figures:** `fig4_htm_sweep.png` is now generated correctly from the fully
  populated `summary_loop2_htm.csv` file.

# References

Baqaee, David Rezza, and Emmanuel Farhi. "Network Effects and Sectoral
Multipliers." *American Economic Review* 110, no. 12 (2020): 3712-58.
<https://doi.org/10.1257/aer.20181806>
