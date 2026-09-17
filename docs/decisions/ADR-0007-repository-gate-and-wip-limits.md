# ADR-0007 — Repository gate and WIP limits

- **Status:** accepted
- **Date:** 2026-09-17
- **Supersedes:** —
- **Related:** ADR-0001, ADR-0003, ADR-0004, ADR-0006, `scripts/check_repo.jl`,
  `scripts/status.jl`, `tests/test_check_repo.jl`, `registry/preregistration.toml`

## Context

Phase 3 built the experiment entry point (`experiments/run.jl`), run manifests,
and preregistration (ADR-0006), but enforcement was split: `run.jl` gates only
its own invocation (preregistration match, no-overwrite), while `status.jl`
validates only registry references. Nothing checked that manifests, index,
scenarios, designs, and the kernel include graph agree with each other — and
nothing stopped two batches from interleaving `running` rows in
`registry/scenarios.csv`. The tooling contract tests (independent verification)
also needed a machine-checkable gate to test against.

## Decision

**`scripts/check_repo.jl` is the mandatory pre/post-batch gate.** It runs ten
checks over an explicit `root` and exits 1 on any violation; every check is a
testable function so `tests/test_check_repo.jl` can point it at fixtures:

1. `board` — the status board rebuilds with 0 warnings and `docs/status.md`
   matches (reuses `build_board`/`normalize` from `status.jl`, which is now
   include-safe via a `PROGRAM_FILE` guard).
2. `closures` — `registry/closures.toml` is bidirectional with the kernel:
   id sets agree with `closure_ids()`, axes agree with `closure_axis()`,
   non-`idea` closures have constructors, `ZETA` has none with status `idea`.
3. `closure-subtypes` — every `struct X <: AbstractLaborClosure` /
   `AbstractFinancing` in `src/` appears in some closure entry's `symbols`,
   except `ExogenousLaborClosure` (legacy CES callback wrapper) and
   `NoFinancing` (the unfinanced baseline, explicitly inadmissible per
   ROADMAP §4.3) — an explicit allowlist, not taxonomy closures.
4. `kernel-reachability` — every `src/**/*.jl` is reachable from
   `src/BeyondHulten.jl` via `include()` (no zombie files; ADR-0001).
5. `runs` — every `runs/<id>/` has a parseable `manifest.toml` with the
   ADR-0004/ADR-0006 keys and tables (`run_id` matches the directory, status
   in `running`/`executed`/`provisional`/`failed`, `log.txt` journal present),
   `runs/index.csv` has the exact header with exactly one row per manifest and
   vice versa, and each manifest has a `scenarios.csv` row with matching
   `labor`/`financing`.
6. `runs-visible` — every `executed`/`provisional`/`failed` scenario row is
   manifest-backed, except the six pre-manifest historical `cbase2-v3` rows
   (explicit allowlist; ADR-0004 never retrofits invented manifests).
7. `tracked-artifacts` — no tracked, non-frozen file is a generated artifact
   (`Manifest.toml`, `*.log`, `.ipynb_checkpoints/`, `plots/*.png`,
   `output/` contents, `runs/` contents beyond the index and manifests).
8. `kernel-copies` — no tracked, non-frozen file `include`s a `cbase2/src`
   path (ADR-0001).
9. `preregistration` — every `registry/preregistration.toml` record matches
   the SHA-256 of its `experiments/designs/<design>.toml`.
10. `wip-limit` — at most one scenario row is `running` at a time.

**WIP limit.** One `running` scenario at a time: batches are sequential, so a
crash can never leave two rows claiming the worker. The operator clears a
stale `running` row explicitly (to `failed` with a note, or back to `planned`).

**Sandbox/promotion path.** Experimental code lives outside `src/` (scratch
directories, notebooks, test-only smoke designs). It enters the kernel only
promoted — registry entry, tests, docs — or it is archived with a DE record.
`check_repo.jl` checks 3/4 enforce the boundary from the `src/` side.

## Consequences

- The batch workflow is: `check_repo.jl` (clean) → `--preregister` (commit) →
  `--design` → `check_repo.jl` (clean) → `status.jl` → commit. Both gates must
  pass; `status.jl --check` remains the board-staleness CI gate.
- The current `matrix_5x3` batch stays `planned`: the real-data reference
  continuation stalls at the first rung and `--design` refuses loudly rather
  than warm-starting from a bad root. That is the recorded solver open item,
  not a gate failure.
- Frozen zones are read off `registry/freeze.toml` by checks 7/8, so a future
  superseding ADR that moves a zone updates the gate automatically.

## Enforcement

- `scripts/check_repo.jl` (exit 1 on violations); `tests/test_check_repo.jl`
  pins the real repo clean plus negative fixtures for the subtype, zombie,
  manifest/index, visibility, preregistration, and WIP checks.
- `AGENTS.md` requires the gate before and after every batch.
