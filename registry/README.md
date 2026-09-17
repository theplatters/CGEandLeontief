# `registry/` — machine-readable tracking (Phase 0)

This directory is the **single source of truth** for what exists in the closure
exploration and what state it is in. Other documents link here; they must not
restate status. See `docs/decisions/ADR-0003-registry-single-source-of-truth.md`.

| File | Contents | Maintained by |
| --- | --- | --- |
| `closures.toml` | Labour and financing closures: formulation, status, code symbols, references, tests, evidence, open gates | Humans/agents by hand; changes require an ADR or a linked run |
| `scenarios.csv` | Experiment runs: the planned 5×3 matrix plus the recorded cbase2-v3 verification runs (provisional/failed) | Humans/agents; from Phase 3, `experiments/run.jl` appends rows |
| `freeze.toml` | Frozen snapshots and read-only zones with recorded commits and git tree hashes | Only via ADR (Phase 0: ADR-0001) |

## Status vocabulary — closures (`closures.toml`)

```
idea -> spec -> implemented -> tested -> validated
                                   \-> rejected     (must link docs/dead-ends/DE-*)
                                   \-> superseded   (must name the successor id)
```

- **idea** — recorded in a plan/ADR, no equations agreed.
- **spec** — formulation and interpretation agreed and written down, no code.
- **implemented** — code exists and runs the benchmark/shock.
- **tested** — repo test suite covers the closure contract and passes.
- **validated** — the Phase 4 gates of `ROADMAP.md` pass, with a recorded run
  manifest.
- **rejected** — tried and abandoned; the entry must link a `DE-` record.
- **superseded** — replaced; the entry must name the replacing closure id.

## Status vocabulary — scenarios (`scenarios.csv`)

```
planned | running | executed | provisional | failed | superseded | cancelled
```

- **executed** — ran and passed the gates recorded with it.
- **provisional** — ran but the result is not gate-clean; the caveat is in `notes`.
- **failed** — ran and did not pass; evidence is kept (never silently dropped).
- **superseded** — replaced by a later `run_id`.

## Schema — `closures.toml`

```toml
schema_version = 1
generated_at  = "YYYY-MM-DD"
generated_by  = "..."

[labor.BF]                 # axis table: [labor.*] or [financing.*]
id             = "BF"      # stable id used in scenarios.csv and designs
name           = "..."
formulation    = "..."     # plain text, Greek letters allowed
interpretation = "..."
status         = "implemented"
symbols        = ["MobileLaborCES"]        # Julia types/functions, grep-able
files          = ["src/mobile_labor.jl"]   # repo-relative paths
references     = ["docs/DOCS_ASSESSMENT.md"]  # where the economics is argued
tests          = ["tests/test_mobile_labor.jl"]
evidence       = ["cbase2/process_comments.md"]  # runs/notebooks/scripts
dead_ends      = ["DE-0001"]   # optional
adrs           = ["ADR-0002"]  # optional
open_gates     = ["..."]       # what blocks tested/validated
notes          = ""
```

## Schema — `scenarios.csv`

Header:

```
run_id,design,status,labor,financing,eta,eta_s,theta,epsilon,sigma,shock,magnitude,data_vintage,evidence,commit,actor,notes
```

- `run_id` — stable, unique, `<design>-<labor>-<financing>[-<variant>]`.
- `labor`/`financing` — ids from `closures.toml`.
- `eta` — BF reallocation parameter; `eta_s` — BETA supply elasticity.
- Unknown or not-yet-pinned cells: `TBD`. Empty cells are allowed.
- Fields containing commas must be double-quoted.

## Generate the status board

```bash
julia --project=. scripts/status.jl          # writes docs/status.md
julia --project=. scripts/status.jl --check  # fails if the board is stale
```

`docs/status.md` is generated from this directory and must never be edited by
hand. Broken references (missing files, unknown ids, stale run directories) are
reported in the board's warnings section.
