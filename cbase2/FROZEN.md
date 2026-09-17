# FROZEN — cbase2 submission/analysis snapshot

| | |
| --- | --- |
| **Status** | frozen (read-only) |
| **Frozen** | 2026-09-17 |
| **Recorded commit** | `dfd60a21f3d816ae9a8d3c55eaa8a3c7a8f45fc7` (`dfd60a2`) |
| **Tree hash** | `a8241e98ef61e572910c3a43d49de2fecc32b9e2` |
| **Tracked files** | 45 |
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
