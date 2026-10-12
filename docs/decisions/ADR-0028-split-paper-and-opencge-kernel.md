# ADR-0028 — Split into paper repo and OpenCGE kernel repo

- **Status:** accepted
- **Date:** 2026-10-12
- **Supersedes:** —
- **Related:** ADR-0001 (one kernel — its "one kernel in root `src/`" layout is
  superseded by this split; the record itself is now summarized in the kernel's
  `DECISIONS.md`), ADR-0003, ADR-0004, ADR-0005 (layout superseded in part),
  ADR-0007

## Context

The monorepo held two products with different lifecycles: a CGE kernel
package (model, closures, diagnostics — the `src/`, `ext/` and kernel test
files) and a paper with its run/registry provenance (manuscript sources,
`runs/`, `registry/`, `experiments/`, `baseline/`, the gate and reproduction
tooling). The kernel is reusable software that other projects should be able
to depend on; the paper is a dated, gate-checked record of executed runs whose
claims must stay traceable to manifests. One repository forced one version
number, one commit stream and one set of operating rules onto artifacts that
change for entirely different reasons (a kernel refactor vs. a manuscript
pass), and it made the kernel unshareable: depending on "the model" meant
depending on a paper repository full of frozen snapshots and run records.

## Decision

- The kernel becomes its own repository at `../opencge`
  (https://github.com/theplatters/opencge), package/module name **`OpenCGE`**
  (new UUID `ce182077-ed3b-4f5d-9409-09ef275840d9`). The git history of the
  kernel files is preserved there via `git filter-repo` (199 commits carried
  over), and the extraction point is tagged **`v0.1.0`**.
- This repository becomes the **paper repo**: the manuscript and revision
  sources, the run/registry provenance, the reproduction baseline and the
  tooling around them. It depends on the kernel via
  `Pkg.add(url, rev="v0.1.0")`, pinned in the **committed root
  `Manifest.toml`** (`repo-rev = "v0.1.0"`, `git-tree-sha1`, `version`).
- Kernel files are removed from this repository **forward-only** — no history
  rewrite here; the pre-split history of this repo stays intact and remains
  the recovery source for the migrated records.
- The kernel carries **no registry/MDR system**. Its records are converted to
  plain `STATUS.md` + `DECISIONS.md` in the kernel repository: one-to-two
  sentence summaries of the kernel-relevant ADRs (0001, 0005, 0010, 0011,
  0014–0017, 0019–0022) and DE records (0001–0004, 0007–0011), with the full
  texts living on in this repo's git history.
- Cross-repository references use the **`OpenCGE@…`** convention in
  `registry/` (schema v2, documented in `registry/README.md`):
  `OpenCGE@<path>` for kernel files/tests, `OpenCGE@<id>` for kernel
  ADR/DE records. Unqualified paths and ids stay local to this repo.

## Consequences

- **An OpenCGE rev bump — not a `src/` change — starts a new run
  generation.** Kernel implementation changes now happen in `../opencge` and
  land here only via a Manifest pin bump (`Pkg.add(url=…, rev=…)`), which
  invalidates the provenance of existing runs exactly as a kernel code change
  used to: re-run as a new generation (`-v2`, `-v3`, …) rather than mixing
  generations in one table.
- **Run manifests record the OpenCGE pin** — `version`, `repo-url`,
  `repo-rev` and `git-tree-sha1` read from the committed `Manifest.toml`
  into every manifest's `[provenance]` block — so every run is traceable to
  the exact kernel tree that produced it.
- **`Pkg.test()` no longer applies here** (this repo is a plain environment,
  not a package). The paper/reproduction suite runs with
  `julia --project=. tests/runtests.jl`; the kernel's own test suite runs in
  the kernel repository.
- **Kernel-reachability and closure-subtype gate checks move to the kernel's
  CI** (the kernel include graph and the `AbstractLaborClosure` /
  `AbstractFinancingClosure` subtype audit are no longer observable from this
  repo); the paper-side gate keeps board, registry-bidirectionality (now
  against `OpenCGE.closure_ids` & friends), run/manifest/index consistency,
  failed-run visibility, tracked-artifact, preregistration and WIP checks.
- The status vocabulary, the ADR/DE/IDEA registers and the lab log stay in
  this repository for everything paper-side; the kernel-side records are
  summarized in OpenCGE's `DECISIONS.md` (numbering kept, never renumbered).

## Enforcement

- `scripts/check_repo.jl` — the repository gate (8 checks; the two
  kernel-structure checks are delegated to the kernel's CI).
- `scripts/status.jl --check` — the generated board must be current and show
  0 warnings (bare `ADR-`/`DE-` ids in `registry/closures.toml` must resolve
  to files kept here; `OpenCGE@…` ids are cross-repo and skip the local
  existence check).
- `scripts/reproduce.jl --check` — the pre-split reproduction baseline in
  `baseline/` must be reproduced within **1e-8** (abs + rel) against the
  pinned kernel.
- The committed root `Manifest.toml` — the OpenCGE pin itself (rev
  `v0.1.0`); it is the only tracked "generated" artifact (exempted in the
  tracked-artifacts check) and any change to it is a generation-level event.
