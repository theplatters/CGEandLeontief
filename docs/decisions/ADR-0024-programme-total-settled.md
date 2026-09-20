# ADR-0024 — Programme total: G0 = 40,300 kept while the size question stays open

- **Status:** proposed (records the interim course; the underlying question
  is open)
- **Date:** 2026-09-20
- **Supersedes:** —
- **Related:** ADR-0003 (registry single source of truth), ADR-0004 (runs),
  ADR-0006 (preregistration — G0 is pinned by the design file's SHA-256);
  `experiments/designs/matrix_5x3_v6.toml` and `matrix_5x3_v9.toml`
  (`[programme] total_eur_m = 40300.0`);
  `cbase2/data_raw/impulses.csv` (the 27 x 71 annual impulse, 2024-2050);
  `cbase2/process_comments.md` (the factor-2 flag, 2026-09-15);
  `docs/log/2026-09.md` (the 2026-09-20 measurement); `docs/CONCISE_SUMMARY.md`
  (open item 3).

## Context

The design's programme magnitude is `total_eur_m = 40300.0` — the paper's
"40.3 bn at 2019 prices" — and every executed cell of the matrix and the
sectoral family is pinned to it. The raw impulse table behind the
programme vector (`impulses.csv`, 27 annual rows 2024-2050, 71 VGR
sectors) contains no row that equals 40.3 bn: its 2024 row sums to
80,760,315,517 (80.76 bn). The cbase2 session of 2026-09-15 flagged the
factor-2 gap with two candidates: the 2046 row (40,371,886,437, "the
decaying programme's mid-range year") and an exact halving of the 2024 row
(40,380,315,779).

Measured on 2026-09-20 (full horizon): sum 1,590,122,645,111; **mean
58,893,431,300 (58.89 bn = 1.945 % of GDP)**; median 65,873,835,779; peak
80,760,315,517 (2024); tail 18,081,686,729 (2050). Composition is broadly
stable against the 2024 vector (max |Delta-psi| = 0.092 over the horizon,
worst year 2050, when machinery's share rises to 13.5 %). The two
candidate rows near the paper's figure remain the only single rows in the
neighbourhood of 40.3 bn; the horizon average is ~46 % larger than G0.

## Decision (interim, open issue)

**Keep G0 = 40,300 EUR m for the executed generations and the current
manuscript pass. The programme-size question is acknowledged open and is
tracked by this ADR.** No re-anchor, no re-run, no re-mint while the
question is open. The system is homogeneous in the programme level —
shares, signs and rankings are unaffected by G0 and every relative metric
scales with it (recorded in CONCISE_SUMMARY open item 3) — so the interim
course leaves the science untouched and only defers a manuscript-level
number.

Options recorded for settling the question, with their costs:

| Option | Content | Cost |
| --- | --- | --- |
| A. Narrative-only | The manuscript states the average annual impulse (58.9 bn, ~1.9 % of GDP, 2024-2050) as the programme's time-profile and the model's evaluated programme separately (40.3 bn, 1.33 % of GDP); linearity legitimises the scaling statement | one paragraph; no runs |
| B. Re-anchor to the mean | New generation (`matrix_5x3_v10`, 33 cells) at G0 = 58,893; every paper number rescales by 58.893/40.3 = 1.462 | fresh design + preregistration + generation + all paper numbers |
| B'. Level-robustness batch | A handful of headline cells (5 labour rows x F2, GAMMA-F3) at G0 = mean, demonstrating the scaling numerically | ~8 cells |
| C. Re-read the source | The decisive evidence is the price basis and coverage of `impulses.csv`: 2019 constant prices vs nominal per vintage year, and whether the 2024 row double-counts (one of two funding sources — the halving conjecture has never been tested against the source) | document-level, may invalidate the average's meaning |

## Open items

- The price basis of `impulses.csv` (constant 2019 vs year-specific
  nominal) is undetermined; it changes what "the average" means. Unresolved
  since 2026-09-15.
- The halving/funding-split conjecture has not been checked against the
  source document.
- If option B is chosen, which average (horizon mean, 2025-2050 mean,
  mean vector with its 2050 composition drift) must be fixed in the design.
- Manuscript wording for the programme's scale while this ADR is open.

## Consequences

- No `src/`, registry, design or preregistration change. The open item
  stays visible in `docs/CONCISE_SUMMARY.md` (open item 3) and in this
  ADR; no DE record (the question is not abandoned).
- Any future G0 change requires a new design file and preregistration
  (ADR-0006) and, because it changes the programme level of every cell, a
  new generation under ADR-0004.

## Enforcement

- While this ADR is open, no design file's `[programme]` block is edited.
- When the question is settled, this ADR is accepted with the resolution
  (if G0 stays) or superseded by the re-anchoring ADR (if it moves).