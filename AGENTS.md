# AGENTS.md — operating contract

## What this repository is

The **paper repo** for the production-network model of green public
investment (BeyondHulten, Baqaee–Farhi lineage) used to compare labour-market
and financing closures: the manuscript, the run/registry provenance, and the
reproduction materials. The model kernel is the separate `OpenCGE` package
(ADR-0028), pinned here via the committed root `Manifest.toml` (tag
`v0.2.1`); there is exactly one kernel, and it is not in this tree
(ADR-0001). Research plan and validation gates live in `ROADMAP.md`; the
documentation audit, closure catalogue and the executed results are in
`docs/DOCS_ASSESSMENT.md` (currently Version 8).

## Layout

| Path | What it is |
| --- | --- |
| `OpenCGE` (external, <https://github.com/theplatters/OpenCGE>; local `../opencge`) | The model kernel — model, closures, diagnostics. Package/module `OpenCGE`, one kernel only (ADR-0001, summarized in the kernel's `DECISIONS.md`); pinned in `Manifest.toml` at rev `v0.2.1`. Kernel implementation changes happen there and land here only via a Manifest pin bump (ADR-0028). |
| `tests/` | The paper/reproduction test suite (entry `tests/runtests.jl`): real-data and tooling contracts. The kernel's own tests live in `OpenCGE` (`OpenCGE@tests/…`). |
| `baseline/` | The pre-split reproduction baseline: everything the paper pipeline produced before the kernel moved out, reproducible within 1e-8 by `scripts/reproduce.jl --check`. |
| `Manifest.toml` | The committed pin of the `OpenCGE` dependency (`repo-rev = "v0.2.1"`, `git-tree-sha1`) — the only tracked "generated" artifact, and a generation-level object (ADR-0028). |
| `scripts/` | Repository tooling — `scripts/status.jl` generates the status board, `scripts/check_repo.jl` is the pre/post-batch gate (ADR-0007), `scripts/reproduce.jl` is the reproduction baseline tool. |
| `registry/` | Machine-readable single source of truth: `closures.toml`, `scenarios.csv`, `freeze.toml`, `preregistration.toml`. Schema: `registry/README.md`; kernel paths and kernel doc ids are repo-qualified `OpenCGE@…` (schema v2). |
| `docs/` | `status.md` (generated status board), `decisions/` (ADRs), `dead-ends/` (DE register), `ideas/` (unexplored-idea register, `IDEA-NNNN`; see its README for the status vocabulary), `log/` (lab session log), `archive/` (closed historical docs). |
| `experiments/` | Single run entry point `run.jl` plus pinned designs in `experiments/designs/` (schemas: `experiments/README.md`; ADR-0006). |
| `runs/` | One `runs/<run_id>/manifest.toml` + `log.txt` (`solution.csv` on success) per run; committed `runs/index.csv` holds one row per run (ADR-0004). Only the index and manifests are tracked. |
| `paper/tables/` | Paper-facing accounting and flow tables generated from the run manifests (cite `run_id`s; never hand-typed numbers). |
| `paper/`, `revision/` | Manuscript sources (`revision/main.tex` is the revision build, `latexmk -pdf`). |
| `data/` | Input data; mostly gitignored, a few small reference files are tracked. |
| `archive/src-orphans/` | Read-only archive of abandoned sources that were never wired into the module. |
| `.opencode/skills/` | Agent workflow skills `tracking`, `closures`, `experiments` (byte-identical copies in `.claude/skills/`). |
| `cbase2/`, `bf_replication/`, `bf_replication2/`, `Replication Files/`, `Dokumente/`, `Notebooks/`, the archived legacy notebooks (`docs/archive/notebooks/`), and `archive/` | Frozen/read-only snapshots and archives (ADR-0001, ADR-0009); never edited in place. |

## Commands

From the repository root:

| Purpose | Command |
| --- | --- |
| Install dependencies | `julia --project=. -e 'using Pkg; Pkg.instantiate()'` |
| Kernel load check | `julia --project=. -e 'using OpenCGE'` |
| Test suite | `julia --project=. tests/runtests.jl` |
| Reproduction baseline | `julia --project=. scripts/reproduce.jl --check` (1e-8 gate against `baseline/`) |
| Regenerate status board | `julia --project=. scripts/status.jl` |
| CI gate (stale board) | `julia --project=. scripts/status.jl --check` |
| Repository gate | `julia --project=. scripts/check_repo.jl` (0 violations required, before and after every batch) |
| List design cells | `julia --project=. experiments/run.jl --list <design>` |
| Preregister a design | `julia --project=. experiments/run.jl --preregister <design> [--actor NAME]` (commit the record first) |
| Run a design | `julia --project=. experiments/run.jl --design <design> [--cell <run_id>] [--cells a,b,c] [--runs-dir DIR] [--budget-seconds N] [--actor NAME]` |
| Build the revision | `latexmk -pdf` in `revision/` (entry `revision/main.tex`) |
| Optional plotting | `julia --project=. -e 'using Pkg; Pkg.add("GLMakie")'`, then `using GLMakie` (the plotting extension lives in `OpenCGE` and loads lazily). |

Julia ≥ 1.9 is required. The root `Manifest.toml` **is committed on
purpose**: it pins `OpenCGE` (rev `v0.2.1`). Change it only as a deliberate
generation-level event (see the generation rule below).

`tests/runtests.jl` includes `test_check_repo.jl` ("real repo is clean"), so
the board and every preregistration record must be current *before* the suite
is run: regenerate the board (`scripts/status.jl`) and re-preregister
(`run.jl --preregister <design>`) first, or a green-code suite reports a red
gate.

## Model invariants pinned by ADR (do not re-derive, do not drift)

The implementations live in the `OpenCGE` kernel; the paper cites these
invariants and must not drift from them. Record ids follow the `OpenCGE@`
convention (`registry/README.md`): `OpenCGE@ADR-…`/`OpenCGE@DE-…` are
summarized in the kernel's `DECISIONS.md`, unqualified ids are paper-repo
records in `docs/`.

- **Mobile labour closure** (OpenCGE@ADR-0019, superseding
  OpenCGE@ADR-0010's N−1 decision): every regime enforces **all N
  goods-market clearings**. The mobile η = 1 system carries the endogenous
  net external transfer `F` (canonical vector `[p; y; w; F]`,
  `E = (1−τ)·w·ΣL + F`); the η = 0 endpoint pins `F = 0`; the fixed-wage
  system clears all N with `F ≡ 0`. Employment allocation uses the endpoints
  η ∈ {0, 1} only; the interpolated η* and the allocation wedge are retired.
- **A-bill intermediate-bill accounting** (ADR-0012, ADR-0013):
  intermediate demand is charged with the **domestic** bill `A_bill`
  (raw-table row 73). The other two components of the purchaser-price bill
  are explicit leaks booked in `external_balance_canary`: imported
  intermediates `M_int` (row 74) and product taxes on intermediate use
  `T_int` (row 75). Row 76 is their sum and
  `ΣA + ΣM_int + ΣT_int = Σλ − 1` holds to 2e-16. Leaving row 75 unbooked
  makes the canary short by exactly `ΣT_int` (2.5745e-2 of GDP).
- **The external-account identity is the acceptance test**
  (OpenCGE@ADR-0019): with zero legacy manna, every η = 1 solution satisfies
  `S + T_int + M − (I+X) = F + B_gov` (`B_gov = Σp·g` under F3, else 0) and
  `market_clearing_residuals ≈ 0`; the fixed-wage and η = 0 endpoints clear
  all N markets too, with the η = 0 canary carrying the documented
  factor-market gap. The legacy-manna path adds `p·(A+G)` to the gap and `F`
  absorbs it. Never threshold-fit it: if it fails, a term is missing from the
  accounting.
- **Labour supply is on the real wage** (OpenCGE@ADR-0014, OpenCGE@DE-0004):
  `L^s = L̄·[(w/P)/(w₀/P₀)]^{η_s}`, deflated by the CPI.
- **Scale determinacy is verified** (OpenCGE@ADR-0014,
  OpenCGE@ADR-0017): the fixed-wage η = 1 system is admitted only when
  `max(A_bill/λ + (1−m)(1−s)·fs) < 1` (the round-gain criterion). Closed
  fixtures sit at exactly 1 and are rejected — including when manna is
  present, since manna is a constant and cannot remove a unit root.
- **Cell metrics must not depend on the solver's stopping point**
  (OpenCGE@ADR-0015): the polish is monotone and targets ~1e-10; the
  acceptance gates stay at 1e-6 (fixed) / 1e-5 (mobile).
- **The cbase2 solver ladder is retired** (OpenCGE@DE-0010): the standing
  `solve()` with its residual-gated LM polish covers the pipeline.

## Operating contract

### Before changing anything

1. Read `docs/status.md` (human entry point) and the relevant entries in
   `registry/` (`closures.toml`, `scenarios.csv`, `freeze.toml`).
2. Read the ADRs in `docs/decisions/` and the dead ends in `docs/dead-ends/`
   that touch the area; kernel-side decisions are summarized in the
   `OpenCGE` kernel's `DECISIONS.md`.
3. Load the matching skill from `.opencode/skills/` when the task matches:
   - `tracking` — status, decisions, dead ends, log;
   - `closures` — closure formulations, promotion, registry entries;
   - `experiments` — runs, scenarios, manifests.
4. Frozen zones are read-only (ADR-0001). Never edit a frozen file in place,
   not even to "fix" it: a change inside one requires a superseding ADR.
   Frozen hashes are recorded in `registry/freeze.toml`.
5. Kernel implementation changes do **not** happen here: they are made in
   `../opencge` (the `OpenCGE` repository) and land in this repo only via a
   Manifest pin bump (`Pkg.add(url=…, rev=…)` + committed `Manifest.toml`).
   This repo owns `registry/closures.toml`, `registry/scenarios.csv` and the
   paper-side tests.

### Status transitions require evidence

| Status | Evidence required |
| --- | --- |
| `implemented` | Code exists in the listed file and runs. |
| `tested` | The listed repo tests pass. |
| `validated` | `ROADMAP.md` Phase 4 gates pass, with a recorded run manifest. |
| `rejected` | A linked `docs/dead-ends/DE-*` record exists. |
| `superseded` | The replacing id is named. |

Never upgrade a status without the evidence. Failed and provisional runs stay
visible (ADR-0004).

### After changing something

- Update `registry/*`: `closures.toml` status/symbols/tests/open gates,
  `scenarios.csv` rows; `freeze.toml` only via ADR.
- Regenerate `docs/status.md` (`julia --project=. scripts/status.jl`); never
  hand-edit it, and it must show **0 warnings**.
- Add an ADR for a decision, a DE record for an abandoned approach, an
  `IDEA-NNNN` note under `docs/ideas/` for a route worth remembering but not
  tested (an idea is not a decision: it carries no registry entry and no
  weight in the paper, and promoting it follows the ADR path), and one entry
  in `docs/log/YYYY-MM.md` for the session.
- Revising `docs/DOCS_ASSESSMENT.md`: increment the version, colour every new
  or changed word with `\textcolor{revisionV<N-1>}{...}` (Version 5 = blue),
  add the line to the top version block and an entry to the Revision Log —
  see the `md-style-corrections` scheme. Inside `\textcolor{...}` write plain
  ASCII math (no `$...$`) and `\texttt{}` with escaped underscores instead of
  backticks; never wrap a markdown table in `\textcolor`.
- Never restate status in prose (ADR-0003); other documents link to
  `registry/` or `docs/status.md`.

### Runs

- Every experiment is a scenario in `registry/scenarios.csv`; ids follow
  `<design>-<labor>-<financing>[-<variant>]` (ADR-0002).
- Workflow per batch: register the scenario row (`planned`) → preregister the
  design (`run.jl --preregister`, commit `registry/preregistration.toml`) →
  run via `experiments/run.jl --design` (the only entry point; it refuses on a
  preregistration mismatch before creating any run dir) → `run.jl` updates the
  scenario row and `runs/index.csv` on every status transition.
- Run `scripts/check_repo.jl` before and after every batch; it must report
  0 violations. It enforces the board/warnings gate, registry bidirectionality
  (against the pinned `OpenCGE` kernel registry), manifest/index/scenario
  consistency, preregistration integrity, and the WIP limit.
- WIP limit: at most one scenario row may be `running` at a time (ADR-0007).
- Failed and provisional runs stay visible with their caveats; never drop
  them.
- Paper tables and figures cite `run_id`s (ADR-0004).
- **Generation rule:** an `OpenCGE` rev bump invalidates the provenance of
  existing runs (as a `src/` change used to): batch such changes and re-run as
  a new generation (`-v2`, `-v3`, …) rather than mixing generations in one
  table. Every run manifest records the OpenCGE pin (version, repo-url,
  repo-rev, git-tree-sha1) in its `[provenance]` block.
- Current state: the matrix is executed on the full-71 A-bill calibration in
  ten generations — `matrix_5x3` (v1, two cells blocked by the retired
  heuristic guard), `matrix_5x3_v2` (ADR-0014), `matrix_5x3_v3` (ADR-0015),
  `matrix_5x3_v4` (ADR-0018 measurement), `matrix_5x3_v5` (ADR-0019,
  external-account closure), `matrix_5x3_v6` (ADR-0020 option C, the η = 0
  sectoral-wage endpoint), `matrix_5x3_v7` / `matrix_5x3_v8` (the record of
  the two ADR-0022 promotion aborts), `matrix_5x3_v9` (first complete
  sectoral generation) and **`matrix_5x3_v10`** (C1-corrected measurement,
  33/33 executed, all gates pass — the generation to cite), plus the
  supply-side identification arm `supply_etas_s1` / `_prog` / `_unif`
  (63 cells each, 189 total, ADR-0023), and the Phase 7 pricing generation
  on OpenCGE `v0.2.0` (ADR-0029, ADR-0030; pricing closures `P1`/`P2`):
  `pricing_capacity_v1`, `pricing_markup_v1`, `pricing_profit_v1` and the
  mobile robustness row `pricing_mobile_v1` (42 cells, all executed). The capacity
  door is re-minted on OpenCGE `v0.2.1` (ADR-0031, value-added-scaled
  exponent, homogeneous and services-zero δ profiles): `pricing_capacity_v2`
  and `pricing_mobile_v2` (27 cells, all executed; they supersede the v1
  capacity cells for citation). Resolved since: the BETA row is
  identified under a supply shock and by the ADR-0022 sectoral form, and the
  η = 0 endpoint carries sectoral wages (ADR-0020 option C). Open items, in
  `registry/closures.toml`, `docs/DOCS_ASSESSMENT.md` §4.2 and the ADRs: the
  recombination of the two labour margins (allocation × supply elasticity) is
  recorded as a deliberate omission; the raw table's own
  production-vs-expenditure residual (5.387 %), the government-side recycling
  of `T_int`, and the ADR-0021 wage-structure dimension (proposed) are open
  modelling items; the `70s` variant is a documented robustness variant and
  `reduced` is deferred.

### Commits

- `closure(<id>): …` — closure formulation or implementation change
- `exp(<design>): run <id>` — executed run
- `adr(NNNN): …`, `deadend(NNNN): …`, `log(YYYY-MM): …` — tracking records
- Otherwise short imperative subjects. Keep commits focused; never mix
  generated figures with model changes.

## Code and tests

- Match the style of surrounding code.
- Tests live in `tests/`; new files are named `test_<feature>.jl` and included
  from `tests/runtests.jl`. Kernel tests live in the `OpenCGE` package
  (`OpenCGE@tests/test_<feature>.jl`, included from its own `tests/runtests.jl`)
  and are referenced from `registry/closures.toml` via the `OpenCGE@`
  convention.
- Use `@testset` and `isapprox` tolerances, not exact floating-point equality.
- A synthetic fixture cannot test the intermediate-bill identity: its
  calibration is consistent under any split of `(1−fs)·λ` into `A + M + T`.
  Assert that identity on the real table (`tests/test_calibration.jl`); on
  fixtures assert only the canary's linearity and solution invariance.
- Scale determinacy is a property of the fixture, not of the closure: closed
  fixtures (m = s = 0) have unit round-gain column sums and are rejected at
  η = 1; give a fixture `s > 0` when a test needs a determinate fixed-wage
  solve.
- Experimental code lives outside the kernel and outside the tracked tooling
  (scratch, notebooks, `experiments/` smoke designs). It enters `OpenCGE`
  only promoted: registry entry + tests + docs, or it is archived with a DE
  record. No zombie files here either — every tool in `scripts/` and
  `experiments/` is referenced from `AGENTS.md`, `experiments/README.md` or
  the registry; kernel-side reachability is enforced by the kernel's own CI.
- Run `julia --project=. tests/runtests.jl` before finishing work.

## House rules

- No plan files in the repository. Plans live in `ROADMAP.md`, ADRs, or the
  skill/registry workflow.
- Do not commit proprietary or machine-local data or generated artifacts;
  `data/` and generated figures are largely gitignored.
- Empty-vector broadcasting is a trap: once `OpenCGE` is loaded
  (NonlinearSolve/SciML), inference returns `Any` for `Int(::Any)`, so
  `Int.(Any[]) === Any[]`. Use `Vector{Int}(x)` for values that can be empty
  (`drops = []`).
- `final_demand_split` returns sector × category MATRICES: `fd.tot[k]` is a
  row, not a category. Use `vec(sum(fd.tot; dims = 1))`.
