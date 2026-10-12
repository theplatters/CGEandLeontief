# `registry/` — machine-readable tracking (Phase 0)

This directory is the **single source of truth** for what exists in the closure
exploration and what state it is in. Other documents link here; they must not
restate status. See `docs/decisions/ADR-0003-registry-single-source-of-truth.md`.

| File | Contents | Maintained by |
| --- | --- | --- |
| `closures.toml` | Labour and financing closures: formulation, status, code symbols, references, tests, evidence, open gates | Humans/agents by hand; changes require an ADR or a linked run |
| `scenarios.csv` | Experiment runs: the planned 5×3 matrix plus the recorded cbase2-v3 verification runs (provisional/failed) | Humans/agents; `experiments/run.jl --design` updates rows on status transitions |
| `freeze.toml` | Frozen snapshots and read-only zones with recorded commits and git tree hashes | Only via ADR (Phase 0: ADR-0001) |
| `preregistration.toml` | Pinned design hashes (`[designs.<design>]`: SHA-256, timestamp, commit, actor, pinned scalars) | Only by `experiments/run.jl --preregister` (ADR-0006) |

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
schema_version = 2
generated_at  = "YYYY-MM-DD"
generated_by  = "..."

[labor.BF]                 # axis table: [labor.*] or [financing.*]
id             = "BF"      # stable id used in scenarios.csv and designs
name           = "..."
formulation    = "..."     # plain text, Greek letters allowed
interpretation = "..."
status         = "implemented"
symbols        = ["MobileLaborCES"]        # Julia types/functions, grep-able
files          = ["OpenCGE@src/core/equilibrium.jl"]   # see "Cross-repo references"
references     = ["docs/DOCS_ASSESSMENT.md"]  # where the economics is argued
tests          = ["OpenCGE@tests/test_mobile_labor.jl", "tests/test_f1_tilt.jl"]
evidence       = ["cbase2/process_comments.md"]  # runs/notebooks/scripts
dead_ends      = ["OpenCGE@DE-0001"]   # optional
adrs           = ["ADR-0002"]          # optional
open_gates     = ["..."]       # what blocks tested/validated
notes          = ""
```

## Cross-repo references (schema v2)

Since the kernel split (2026-10-12) the model kernel is the separate `OpenCGE`
package (pinned in the committed root `Manifest.toml`) and this repository is
the paper repo. Registry entries therefore use a repo-qualified convention for
everything that lives in the kernel:

- `Name@rest` in `files`, `tests`, `references` and `evidence` — `rest` is a
  path relative to the root of the package `Name`; `OpenCGE@src/core/equilibrium.jl`
  resolves against the pinned OpenCGE checkout
  (`Base.find_package("OpenCGE")` → `<root>/src/OpenCGE.jl`, root =
  `dirname(dirname(...))`), and the existence checks and the symbol grep of
  `scripts/status.jl` run against that pinned kernel. Unqualified paths
  resolve in this repo as before.
- `Name@id` in `dead_ends` and `adrs` — cross-repo document references
  (`OpenCGE@DE-…`, `OpenCGE@ADR-…` live in the kernel's DECISIONS.md). They
  skip the local `docs/dead-ends/` / `docs/decisions/` existence check and
  render on the status board as cross-repo references. Unqualified ids stay
  local to `docs/`.

Concretely: kernel formulation paths and kernel test files, and the ADRs/DE
records that moved to the kernel's DECISIONS.md (ADR-0001/0005/0010/0011/
0014/0015/0016/0017/0019/0020/0021/0022 and the formulation-level DE records),
are qualified `OpenCGE@…`; paper-side test paths, `evidence` (`runs/…`) and
`references` (`docs/…`) stay unqualified.

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

## Run manifests and the run workflow (ADR-0004, ADR-0006)

Every quantitative claim is traceable: claim → `run_id` → manifest →
commit → closure/design → data vintage. Runs go through the single entry
point `experiments/run.jl` (schemas: `experiments/README.md`):

```bash
julia --project=. experiments/run.jl --list <design>
julia --project=. experiments/run.jl --preregister <design> [--actor NAME]
julia --project=. experiments/run.jl --design <design> [--cell <run_id>] [--cells a,b,c] [--runs-dir DIR] [--budget-seconds N] [--actor NAME]
```

Design files live in `experiments/designs/<design>.toml`. `--design`
refuses to start unless the design file's SHA-256 matches
`registry/preregistration.toml` (written only by `--preregister`).
`runs/<run_id>/` holds `manifest.toml` (`running` → `executed`/`failed`),
`log.txt`, and `solution.csv` on success; existing run dirs are never
overwritten (reruns register a `-v2` variant). The committed
`runs/index.csv` holds one row per run (`run_id`, date, design, closures,
status, gate summary, headline metrics, commit); only the index and the
manifests are tracked (`.gitignore`, checked with `git check-ignore`).
`--design` also updates the cell's `registry/scenarios.csv` row (status,
pinned params, `evidence`, recorded HEAD).

## Generate the status board

```bash
julia --project=. scripts/status.jl          # writes docs/status.md
julia --project=. scripts/status.jl --check  # fails if the board is stale
```

`docs/status.md` is generated from this directory and must never be edited by
hand. Broken references (missing files, unknown ids, stale run directories) are
reported in the board's warnings section.

## Repository gate and WIP limit (ADR-0007)

```bash
julia --project=. scripts/check_repo.jl         # eight checks, exit 1 on violations
julia --project=. scripts/check_repo.jl --quiet # violations only
```

Run the gate before and after every batch. It covers the board/warnings check,
closure bidirectionality against the pinned OpenCGE kernel,
manifest/index/scenario consistency, failed/provisional visibility, tracked
artifacts (the root `Manifest.toml` is a required committed artifact and
exempt), living kernel copies, preregistration integrity, and the WIP limit:
at most one `scenarios.csv` row may be `running` at a time. The pre-split
closure-subtype and kernel-reachability checks moved with the kernel sources
and are owned by the OpenCGE kernel's CI.
