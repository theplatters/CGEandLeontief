# DE-0008 — Open-economy calibration without a saving rate (cbase2 v2)

- **Status:** dead end — no interior equilibrium
- **Recorded:** 2026-09-17
- **Origin:** cbase2 open-absorption recalibration, decided 2026-09-15,
  reversed 2026-09-16
- **Related:** closures `F1`–`F3`, `DELTA`; cbase2 v3 structure
- **Evidence:** `cbase2/process_comments.md` 2026-09-16 ("STRUCTURAL FINDING:
  the v2 import margin without a saving/export block has NO interior
  equilibrium"; "The 2026-09-15 decision to skip the saving rate was wrong")

## What was tried

An open-economy demand block with import margins on household and programme
spending, government transfers, and no household saving rate: household
consumes all income and the import leak circulates nowhere.

## Why it seemed promising

Import margins answer the referee point about imports at the margin and at
baseline without adding a saving–investment block; the accounting looked
complete because domestic budget identities held.

## How it failed

The import leak `m̄E` has no offsetting injection, so `Y = C + I + G + X − M`
with a fully spending household and balanced government forces `p'imp = 0`:
income collapses until the leak vanishes. The only solution is the corner
`E = 0`, `L = Σ gG`. All-N market clearing is jointly inconsistent for
`E > 0`; numeric "successes" were spurious (an overwritten clearing equation
parked the leakage in one sector).

## What replaced it

The cbase2 v3 open-economy Keynesian structure: exogenous investment and
exports (injections), household saving rate `s` calibrating the accounting
partner of `I + X − M` (`S = I + X − M` with `T = G`), and import margins on
consumption/programme demand only. `s = 0.398` on the 71-sector data
(0.426 on the 70-sector variant). Budget identities hold exactly at solved
points.

## Revival conditions

None under this accounting. Any alternative demand block must satisfy
`S = I + X − M` (with `T = G`) and assert it as a canary at the solved
equilibrium, not at a reference point.
