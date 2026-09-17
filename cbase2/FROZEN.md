# FROZEN — cbase2 submission/analysis snapshot

| | |
| --- | --- |
| **Status** | frozen (read-only) |
| **Frozen** | 2026-09-17 |
| **Re-pinned** | 2026-09-17 by ADR-0011 (notebook 01 conformed) — see `registry/freeze.toml` for the current `recorded_commit` / `tree_hash` |
| **Original freeze commit** | `dfd60a21f3d816ae9a8d3c55eaa8a3c7a8f45fc7` (`dfd60a2`) |
| **Original tree hash** | `a8241e98ef61e572910c3a43d49de2fecc32b9e2` |
| **Recorded in** | `registry/freeze.toml` · `docs/status.md` |
| **Authority** | `docs/decisions/ADR-0001-one-kernel-no-living-copies.md` |

## Do not edit this folder

Active development continues in the root package (`src/`, `experiments/`).
Remaining cbase2 work — notebooks 04–08 equivalents, BETA/DELTA/F1–F3
promotion, and the v3 open items — moves to the canonical kernel in Phase 2
(see `docs/log/2026-09.md`). If a self-contained submission artifact is
needed later, cut a **new** snapshot from the canonical kernel; do not
revive this copy.

## Recorded kernel DIFFs at freeze

`src/core/interface.jl` and `src/core/mobile_labor.jl` diverged from the
parent `src/` (financing hooks). Recorded in
`cbase2/process_comments.md`, notebook 03 entry; backport is pending and is
resolved by the Phase 2 promotion.

## Open items at freeze

See `registry/freeze.toml` (`open_items`) for the authoritative list:
Törnqvist base still in v2 form; F3 mobile stall and BETA continuation
timeout; S = I + X − M canary mismatch and the omitted N-th clearing
equation (see `cbase2/review.md`); θ bounded at 1 with its own stall;
sector-71/high-self-loop instability and `drop_sectors` normalization;
clamp mass ≈ 8.7% of GDP and stale 71-sector headlines.

`cbase2/documentation.md` and `cbase2/process_comments.md` remain the
pipeline and history record; `cbase2/review.md` is the external review.

## Amendment 2026-09-17 — notebook 01 conformed to the canonical kernel (ADR-0011)

`01_data_wrangling.ipynb` was conformed **in place** under ADR-0011, the single
authorised exception to this freeze:

- the schema is resolved **by label** (`FD_NAMES`, the import and product-tax
  rows) instead of the hard-coded table positions 75:81 / 74 / 75, matching
  `src/core/accounting.jl` (`final_demand_columns`) and surviving a rebuilt
  retained-sector table;
- the `Plots.jl` figure — `Plots` is not a dependency of the `BeyondHulten`
  project — is replaced by an unconditional numeric summary plus an opt-in
  figure (`CBASE2_PLOT_FIGURE=1` with the GLMakie extension);
- an optional, non-fatal cross-check cell compares the notebook's own objects
  against `BeyondHulten.generate_data`.

**The artifact contract is unchanged.** All seven tracked
`data_processed/` artifacts regenerate byte-identically on Julia 1.12.7, the
notebook writes nothing else into this zone, and the cross-check never blocks
those writes. `02_accounting_consistency.ipynb`, `03_financing_closures.ipynb`,
`cbase2/src/` and `cbase2/plots/` are untouched. The freeze record in
`registry/freeze.toml` is re-pinned to the conforming commit; the original
freeze point stays visible above and in git history.

Any further edit to this zone — including conformance edits of notebooks 02 and
03 — requires a new ADR.