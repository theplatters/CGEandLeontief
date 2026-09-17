---
name: tracking
description: Use when changing or checking repository status in the BeyondHulten repo — adding or updating registry entries, closing or rejecting an approach, recording a decision (ADR), a dead end (DE), or a lab-log entry, or when unsure which document is authoritative. Covers registry/closures.toml, registry/scenarios.csv, registry/freeze.toml, the generated docs/status.md board, status vocabularies and evidence gates, freeze zones, and commit message conventions.
---

# Tracking: registry, status, decisions

## The one rule

`registry/` is the single source of truth for what exists and in which state
(ADR-0003). `docs/status.md` is generated from it by `scripts/status.jl` and
must never be hand-edited. Prose documents never restate status; they link to
`docs/status.md` or the registry.

Read first, always:

1. `docs/status.md` — current board (generated).
2. `registry/README.md` — schemas and status vocabularies.
3. `docs/decisions/README.md`, `docs/dead-ends/README.md` — ADR and DE rules.
4. `docs/log/2026-09.md` — lab-log format (newest `docs/log/YYYY-MM.md`).
5. `registry/freeze.toml` — what is frozen or read-only.

## Status vocabularies and evidence gates

Closures (`registry/closures.toml`): `idea -> spec -> implemented -> tested -> validated`, plus `rejected` and `superseded`.

| Status | Evidence required |
| --- | --- |
| `idea` | recorded in a plan/ADR; no equations agreed |
| `spec` | formulation and interpretation written down; no code |
| `implemented` | code exists in a listed `files` path and runs |
| `tested` | listed `tests` cover the contract and the suite passes |
| `validated` | the Phase 4 gates of `ROADMAP.md` pass **and** a recorded run manifest is linked |
| `rejected` | a linked `docs/dead-ends/DE-*.md` record |
| `superseded` | the replacing closure id is named |

Changes to `registry/closures.toml` require an ADR or a linked run. `registry/freeze.toml` changes only via a superseding ADR.

Scenarios (`registry/scenarios.csv`): `planned | running | executed | provisional | failed | superseded | cancelled`.
`executed` = ran and passed its recorded gates; `provisional` = ran but not gate-clean (caveat in `notes`); `failed` = ran and did not pass (evidence kept). Failed and provisional rows are never deleted.

## Update flow

```bash
# 1. edit registry/* (closures.toml, scenarios.csv; freeze.toml only via ADR)
# 2. regenerate the board and inspect the warnings count
julia --project=. scripts/status.jl              # must print "0 warnings"
# 3. verify the committed board is not stale
julia --project=. scripts/status.jl --check      # "docs/status.md is up to date."
# 4. commit registry changes and the regenerated docs/status.md together
```

Warnings mean a dangling reference: missing `files`/`tests`/`evidence` path, unknown closure id in a scenario, missing `DE-`/`ADR-` file, duplicate `run_id`, or a frozen zone modified since freeze. Fix the reference; never delete the row and never edit `docs/status.md` directly. `--check` fails only on staleness, so regenerate and commit the board in the same change.

## ADRs (`docs/decisions/`)

- One file per decision, `ADR-####-short-slug.md`, numbered consecutively, never reused; append-only.
- A decision is immutable once `accepted`. To change it, write a new ADR with `Supersedes: ADR-####` and change the old record's status to `superseded by ADR-####`; do not rewrite its content.
- Status values: `proposed | accepted | superseded | rejected`. Cite related `DE-` ids, registry entries, and code paths.
- Required when changing a closure formulation, a status rule, or the freeze policy; when unfreezing or moving a frozen zone; when adding a closure id. Template and index: `docs/decisions/README.md`.

## Dead ends (`docs/dead-ends/`)

- One file per abandoned approach, `DE-####-short-slug.md`, numbered consecutively, never reused; append-only.
- Required when an approach taken seriously is abandoned (not for typos or ordinary bug fixes; when in doubt, record it).
- Must state what was tried, evidence, what replaced it, and the conditions under which it could be revived. Link the DE from `registry/closures.toml` (`dead_ends = [...]`) and add it to the index table in `docs/dead-ends/README.md`.
- Never rewrite an old DE; a revival gets a new ADR plus a registry note.

## Lab log (`docs/log/YYYY-MM.md`)

- Append-only, chronological, one entry per working session or experiment batch; failed batches are logged too and never deleted.
- Format: intent → what ran (run ids and commit hashes) → observation → decision → next.
- Reference run ids/commits and link ADRs/DEs; do not restate status.

## Freeze rules

- ADR-0001 and `registry/freeze.toml` are authoritative: frozen snapshots and read-only zones are never edited.
- Frozen/read-only zones: `cbase2/`, `bf_replication/`, `bf_replication2/`, `Replication Files/`, `Dokumente/`, `Notebooks/`, the legacy root notebooks (`CobbDouglas.ipynb`, `CompareModels.ipynb`, `DemandShocks.ipynb`), and `archive/` (including the orphan-source archive `archive/src-orphans/`).
- Read frozen sources if needed, but write all changes in the root kernel. Never copy the kernel (ADR-0001).
- Unfreezing or moving a zone requires a superseding ADR and a new freeze record; `scripts/status.jl` warns when a zone differs from its recorded commit/tree hash.

## Commit conventions

- `closure(<id>): …` — closure registry or code change
- `exp(<design>): run <run_id>` — experiment/run record
- `adr(NNNN): …` — decision record
- `deadend(NNNN): …` — dead-end record
- `log(YYYY-MM): …` — lab-log entry
- Otherwise a short imperative subject (`Add …`, `Fix …`, `Refactor …`); keep generated figures out of model commits.

## Routing

| You are about to… | Do this |
| --- | --- |
| state or change status | edit `registry/*`, regenerate `docs/status.md` |
| settle a design question | add an ADR |
| abandon an approach | add a DE record and link it |
| run an experiment | register a scenario (see `experiments` skill) |
| implement or change a closure | root `src/` only (see `closures` skill) |
