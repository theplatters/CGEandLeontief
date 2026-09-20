# ADR-0024 — Programme total: G0 = 40,300 kept while the size question stays open

- **Status:** proposed (records the interim course; the price basis is now
  resolved — see Context; the underlying size question stays open)
- **Date:** 2026-09-20
- **Supersedes:** —
- **Related:** ADR-0003 (registry single source of truth), ADR-0004 (runs),
  ADR-0006 (preregistration — G0 is pinned by the design file's SHA-256);
  `experiments/designs/matrix_5x3_v6.toml` and `matrix_5x3_v9.toml`
  (`[programme] total_eur_m = 40300.0`);
  `cbase2/data_raw/impulses.csv` (the 27 x 71 annual impulse, 2024-2050);
  `cbase2/process_comments.md` (the factor-2 flag, 2026-09-15);
  `docs/log/2026-09.md` (the 2026-09-20 measurement); `docs/CONCISE_SUMMARY.md`
  (open item 3); `paper/main.tex` (the 1.46 deflator footnote);
  `src/BeyondHulten.jl` (`const inflator = 1.46`), `src/core/diagnostics.jl`
  (`impulse_shock`; main-branch `src/util.jl`, revise/revisefinal
  `cbase2/src/core/util.jl`) and `tests/test_model.jl`;
  `experiments/README.md` (§Cell construction).

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
neighbourhood of 40.3 bn; the raw horizon average is ~46 % larger than
G0 — a gap the price-basis resolution below closes.

**Price basis resolved (2026-09-20).** The missing step is in
the repository, in two places. (i) The paper states the deflation explicitly
(`paper/main.tex`, the Hornykewycz.2025 paragraph, main branch and current
tree alike): "about 40.3 bn. EUR in 2019 prices", with the footnote "We
assumed a deflator of 1.46 based on price indices for construction costs as
specified in Hornykewycz.2025. In 2023 prices the shock amounts to about
58 bn. EUR." (ii) The code lineage applies it: the legacy `impulse_shock`
divides the whole impulse table by a scalar `inflator` before use — main
`src/util.jl`, revise/revisefinal `cbase2/src/core/util.jl`, and still in the
current kernel (`src/core/diagnostics.jl`, exported and exercised by
`tests/test_model.jl`) with `const inflator = 1.46` in `src/BeyondHulten.jl`.
Its autonomous-demand component is the horizon mean of the deflated table:
58,893.43 EUR m / 1.46 = 40,337.97 EUR m ≈ 40.3 bn — the paper's figure to
0.1 %.

Three consequences. (i) `impulses.csv` is in **current prices** (the source's
2023 price base): its horizon mean, 58.89 bn, is the paper's "58 bn in 2023
prices", and the design's `total_eur_m = 40300.0` is that mean at 2019
prices. (ii) The factor-2 flag is **peak year vs horizon average** (80.76 vs
58.89) times the **price base** (the 2024 row is 55.3 bn at 2019 prices), not
the halving/funding-split conjecture; the 80.76/2 = 40.38 coincidence is
superseded. (iii) The option-B comparison below mixes price bases: at 2019
prices the mean is 40.34 bn (≈ 1.33 % of GDP), i.e. the current G0 —
re-anchoring to the raw current-price mean would double-count the deflation.

## Decision (interim, open issue)

**Keep G0 = 40,300 EUR m for the executed generations and the current
manuscript pass. The programme-size question is acknowledged open and is
tracked by this ADR.** No re-anchor, no re-run, no re-mint while the
question is open. The system is homogeneous in the programme level —
shares, signs and rankings are unaffected by G0 and every relative metric
scales with it to first order (the sectoral and capacity variants are
nonlinear, so the scaling is approximate; option B' would measure it) — so
the interim course leaves the science untouched and only defers a
manuscript-level number. The **price-basis item is closed** (Context, same
date): the raw impulse is in current (2023) prices, the design's G0 is that
horizon mean at 2019 prices, and the interim course stands with its price
base identified rather than assumed.

Options recorded for settling the question, with their costs:

| Option | Content | Cost |
| --- | --- | --- |
| A. Narrative-only | The manuscript states the average annual impulse **at 2019 prices** (40.34 bn ≈ 1.33 % of GDP, 2024-2050) as the programme's time-profile — the same number as the model's G0, so no reconciliation is needed; the source's current-price mean (58.89 bn) is quoted as the price-base statement | one paragraph; no runs |
| B. Re-anchor to the mean | **A no-op at the correct price base**: the mean at 2019 prices is 40.34 bn ≈ G0. Re-anchoring to the raw current-price mean (58,893) would mix price bases and double-count the 1.46 deflation; withdrawn unless the price base or the composition changes | none as stated; a new generation only if a real level change is chosen |
| B'. Level-robustness batch | Optional: a handful of headline cells at G0 = 40.34 bn (the deflated mean, +0.1 % on G0) to demonstrate level invariance numerically | ~8 cells |
| C. Re-read the source | **Answered in-repo (2026-09-20)**: price basis = current (2023) prices; scalar deflator 1.46; the paper's 40.3 bn is the deflated horizon mean; the halving conjecture is superseded. What remains is the **composition** (ψ is the raw current-price 2024 composition; a constant-2019-price incidence vector needs sectoral deflators) and the unread funding-split question | document-level; the composition is a modelling choice |

## Open items

- **Price basis: resolved 2026-09-20.** Current (2023) prices; the paper's
  40.3 bn is the horizon mean deflated at 1.46; G0 = 40,300 is that figure.
  Evidence: the paper footnote (`paper/main.tex`), the legacy
  `inflator`/`impulse_shock` (`src/BeyondHulten.jl`, `src/core/diagnostics.jl`;
  main `src/util.jl`), and 58,893.43 / 1.46 = 40,337.97 ≈ 40.3 bn.
- **The halving/funding-split conjecture is superseded** by the deflator
  explanation (80.76 / 2 = 40.38 is a numerical coincidence). The
  funding-split question is still unread against the source, but it no longer
  bears on the level.
- **Composition.** ψ is the raw 2024 current-price composition; the legacy
  scalar deflator could not deflate a composition either, so both pipelines
  share the limitation. A constant-2019-price incidence vector needs sectoral
  deflators and would change the incidence results (63.4 % specialised
  construction; the rigid-group rule) — a modelling choice, not housekeeping.
- **Time aggregation.** The legacy code used the first two years (2024-2025)
  for the demand multiplier and the horizon mean for the autonomous component;
  the current pipeline uses the horizon mean as the level, matching the
  paper's "yearly costs". Record in the manuscript.
- If a level change is ever chosen, which average (horizon mean, 2025-2050
  mean, mean vector with its 2050 composition drift) must be fixed in the
  design.
- Manuscript wording for the programme's scale while this ADR is open.

## Consequences

- No `src/`, registry, design or preregistration change. The open item
  stays visible in `docs/CONCISE_SUMMARY.md` (open item 3, v2) and in this
  ADR; no DE record (the question is not abandoned).
- The price base is recorded in this ADR and in `experiments/README.md`
  (§Cell construction); the design files are not edited (their SHA-256 is
  pinned by `registry/preregistration.toml`), so the deflation is
  documented, not re-applied in code.
- Any future G0 change requires a new design file and preregistration
  (ADR-0006) and, because it changes the programme level of every cell, a
  new generation under ADR-0004.

## Enforcement

- While this ADR is open, no design file's `[programme]` block is edited.
- The price base lives in this ADR and in `experiments/README.md`; a future
  generation that wants the deflation explicit in the harness requires a new
  design and preregistration (ADR-0006) and a new generation (ADR-0004).
- When the question is settled, this ADR is accepted with the resolution
  (if G0 stays) or superseded by the re-anchoring ADR (if it moves).