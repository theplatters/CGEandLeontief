# ADR-0001 — One kernel, no living copies

- **Status:** accepted
- **Date:** 2026-09-17
- **Related:** ADR-0003, ADR-0004, `registry/freeze.toml`, `docs/status.md`

## Context

The repository grew four model implementations side by side:

1. the root Julia package (`src/`), declared the reproducible baseline in
   `ROADMAP.md` §3 and described by `AGENTS.md`;
2. `cbase2/src/core/`, byte-identical copies of seven root kernel files plus
   cbase2-only additions, kept synchronized by `cbase2/scripts/diff_kernel.jl`;
3. `bf_replication/src/`, a port of Baqaee–Farhi (2019);
4. `bf_replication2/src/`, a clean-room port of Baqaee–Farhi (2022).

The copies existed for two reasons: `cbase2` had to be self-contained for an
assessor, and the root package pulls heavy dependencies (GLMakie, XLSX, JuMP,
Ipopt) that hinder headless runs. In practice, fixes were applied to one copy
and recorded as "pending backport" in `cbase2/process_comments.md`, and
status claims in the copies drifted from the parent. Every new closure
currently means touching at least two kernels plus notebooks.

## Decision

- The root package (`src/`, `Project.toml`) is the **single canonical
  kernel** for all new development: closures, financing, diagnostics,
  sensitivity.
- Closures and financing rules are **plug-in modules** registered centrally
  (ADR-0002); no new residual system may be built by copying the kernel.
- `cbase2`, `bf_replication`, and `bf_replication2` are **frozen snapshots**.
  Their freeze point (recorded commit, git tree hash, tracked-file count,
  open items) lives in `registry/freeze.toml`. Frozen zones are read-only.
- Remaining `cbase2` work (BETA, DELTA, F1–F3 and the matrix) is promoted
  into `src/` in Phase 2. If a self-contained submission artifact is needed
  later, a **new snapshot** is cut from the canonical kernel — it does not
  become a second living copy.
- Legacy source material (`Replication Files/`, `Dokumente/`, root and
  `Notebooks/` legacy notebooks) stays at its path but is marked read-only in
  `registry/freeze.toml`.

## Consequences

- There is exactly one place to fix a bug and one place to add a closure.
- `cbase2` as a *living* development area ends at the freeze commit; the
  v3 open items recorded there must be resolved in the root package (or
  explicitly abandoned with a dead-end record).
- The drift checker `cbase2/scripts/diff_kernel.jl` retires; freeze hashes
  in `registry/freeze.toml` are checked by the repo checker instead
  (Phase 3, `scripts/check_repo.jl`).
- Any change to a frozen zone requires a superseding ADR. New snapshots get
  a new freeze record.

## Enforcement

- `registry/freeze.toml` is authoritative; `docs/status.md` renders it.
- The repo checker (Phase 3) fails when a frozen zone's tree hash or working
  tree differs from its freeze record.
- `AGENTS.md` states the rule for agents: never edit a frozen zone, never
  copy the kernel.

## Amendment 2026-09-17 — The root notebooks leave their path (ADR-0009)

The Decision line that legacy source material "stays at its path" no longer
applies to the three root notebooks (`DemandShocks.ipynb`,
`CompareModels.ipynb`, `CobbDouglas.ipynb`): ADR-0009 moved them to
`docs/archive/notebooks/`, where they stay read-only. Everything else in
this ADR — in particular the one-kernel rule and the freeze of `cbase2`,
`bf_replication`, `bf_replication2` and the remaining legacy zones — is
unchanged. The original wording is kept per the append-only rule; this
amendment is authoritative where the two disagree.
