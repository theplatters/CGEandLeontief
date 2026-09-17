# DE-0005 — Fixed-base quantity sum called real GDP

- **Status:** dead end — not a quantity index
- **Recorded:** 2026-09-17
- **Origin:** early model versions; mobile-labour system, 2026-09
- **Related:** ADR-0004, `ROADMAP.md` §3 and §4.2
- **Evidence:** `roadmaps/vertdict.md` item 6 ("the mobile model calculates a
  fixed-base sum instead of using the shared Törnqvist index"),
  `ROADMAP.md` §3 (replacement recorded as done for the baseline)

## What was tried

Aggregate real GDP as a fixed-base weighted sum of quantities
(`Σ_i q_i` with base-period weights), recomputed for each experiment.

## Why it seemed promising

It is simple, additive, and was already used in the legacy results and
figures, so it preserved continuity.

## How it failed

A fixed-base sum is base-year dependent, not a quantity index; it cannot be
compared across equilibria with different relative prices. Later, even the
Törnqvist replacement was found to be a **household-consumption** index
rather than GDP (`cbase2/review.md` §3.8), and the cbase2 v3 open items show
the base vector itself must be the actual baseline (`c0_gross`,
`cbase2/process_comments.md` 2026-09-16).

## What replaced it

The shared `tornqvist_quantity_index` helper (`src/util.jl`, exported) for
every model variant, plus a supplementary total-final-demand index in
`cbase2` when financing shifts expenditure between institutions
(`cbase2/process_comments.md` 2026-09-15). Real GDP is reported separately
from welfare (`ROADMAP.md` §4.2).

## Revival conditions

None as a GDP measure. A fixed-base sum may only appear as an explicitly
labelled diagnostic, never as "real GDP".
