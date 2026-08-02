# Repository Guidelines

## Project Structure & Module Organization

`src/BeyondHulten.jl` defines the Julia module and includes the implementation files in `src/`. Keep model-specific logic in files such as `ces.jl`, `cobbdouglas.jl`, or `leontief.jl`; shared data structures and helpers belong in `interface.jl`, `solution.jl`, and `util.jl`. Tests live in `tests/`, with `tests/runtests.jl` as the package test entry point. Exploratory and replication work is stored in root-level and `Notebooks/` Jupyter notebooks. Input data belongs in the gitignored `data/` directory. Treat `Replication Files/` and `Dokumente/` as source/archive material, while `plots/` and `paper/pictures/` contain generated figures. Paper and presentation sources are under `paper/`.

## Build, Test, and Development Commands

- `julia --project=. -e 'using Pkg; Pkg.instantiate()'` installs the dependencies pinned in `Manifest.toml`.
- `julia --project=. -e 'using BeyondHulten'` checks that the package loads and source includes resolve.
- `julia --project=. -e 'using Pkg; Pkg.test()'` runs the test suite through `tests/runtests.jl`.
- `julia --project=.` starts an interactive Julia session in the repository environment; use IJulia/Jupyter for notebooks.

Julia 1.9 or newer is required. Some analyses also require an input-output CSV placed in `data/`; do not commit proprietary or machine-local datasets.

## Coding Style & Naming Conventions

Follow the existing Julia conventions: `CamelCase` for types (`Solution`, `CESElasticities`), lowercase `snake_case` for functions and variables (`read_data`, `real_gdp`), and a trailing `!` for mutating functions. Indent nested blocks consistently with surrounding code and avoid whitespace-only rewrites; tabs are common in current sources. Add docstrings to exported APIs and use explicit module imports when they clarify dependencies. No formatter or linter is currently enforced.

## Testing Guidelines

Use Julia's `Test` standard library. Add focused `@testset` blocks and name new files `test_<feature>.jl`; include them from `tests/runtests.jl`. Cover numerical results with tolerances (`isapprox`) rather than exact floating-point equality. Run both package loading and `Pkg.test()` before submitting changes.

## Commit & Pull Request Guidelines

Recent commits use short, imperative summaries such as `Add plotting functions...`, `Refactor...`, and `Fix...`; occasional `refactor:` prefixes are accepted. Keep commits focused and avoid mixing generated figures with model changes. Pull requests should explain the economic or computational change, list validation commands, identify required data, and link related issues. Include before/after figures when plotting output changes, but do not commit notebook checkpoint files.
