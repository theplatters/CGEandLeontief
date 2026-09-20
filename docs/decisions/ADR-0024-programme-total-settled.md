# ADR-0024 — Programme total: G0 = 40,300 kept (the 2019-price horizon average)

- **Status:** accepted (user decision, 2026-09-20): keep G0 = 40,300 EUR m — the
  deflated horizon mean (40.34 bn to 0.1 %) — as the principal programme; the
  composition item remains tracked (Open items)
- **Date:** 2026-09-20
- **Supersedes:** —
- **Related:** ADR-0003 (registry single source of truth), ADR-0004 (runs),
  ADR-0006 (preregistration — G0 is pinned by the design file's SHA-256);
  `experiments/designs/matrix_5x3_v6.toml` and `matrix_5x3_v9.toml`
  (`[programme] total_eur_m = 40300.0`); `cbase2/data_raw/impulses.csv`
  (the 27 x 71 annual impulse, 2024-2050); `cbase2/process_comments.md`
  (the factor-2 flag, 2026-09-15); `paper/main.tex` (lines 442 and 516: the
  programme statement and its price basis); Hornykewycz, Kapeller, Weber et al.
  (2025), "Carbon neutrality in the residential sector: a general toolbox and
  the case of Germany", npj Climate Action 4(1) (also IFSO WP 41);
  `docs/log/2026-09.md` (the measurements); `docs/CONCISE_SUMMARY.md`
  (open item 3); `src/BeyondHulten.jl` (`const inflator = 1.46`),
  `src/core/diagnostics.jl` (`impulse_shock`; main-branch `src/util.jl`,
  revise/revisefinal `cbase2/src/core/util.jl`) and `tests/test_model.jl`;
  `experiments/README.md` (§Cell construction).

## Context

The design's programme magnitude is `total_eur_m = 40300.0` — the paper's
"40.3 bn at 2019 prices" — and every executed cell of the matrix and the
sectoral family is pinned to it. The raw impulse table behind the programme
vector (`impulses.csv`, 27 annual rows 2024-2050, 71 VGR sectors) contains no
row that equals 40.3 bn: its first row sums to 80,760,315,517 (80.76 bn). The
cbase2 session of 2026-09-15 flagged the factor-2 gap with two candidates: the
2046 row (40,371,886,437) and an exact halving of the 2024 row (40,380,315,779).

Measured on 2026-09-20 (full horizon): sum 1,590,122,645,111; **mean
58,893,431,300 (58.89 bn = 1.945 % of GDP at G0 = mean)**; median
65,873,835,779; peak 80,760,315,517 (2024); tail 18,081,686,729 (2050).
Composition is broadly stable against the 2024 vector (max |Delta-psi| = 0.092
over the horizon, worst year 2050, when machinery's share rises to 13.5 %). The
two candidate rows near the paper's figure remain the only single rows in the
neighbourhood of 40.3 bn; the raw horizon average is ~46 % larger than G0 — a
gap the price-basis resolution below closes.

## The price basis (resolved 2026-09-20)

The factor-2 was a price-basis mismatch, not a data error. The source study
states the programme in **2023 prices**: "an additional yearly investment of
58 bn EUR is needed", "the total investment until 2050 sums up to 3.1 trillion
EUR (in 2023 prices)", and "in the first year of the policy measure, additional
costs of 81 bn EUR are anticipated, about 1.9 % of GDP, which over time
decreases to 0.3 %". The manuscript's own footnote (`paper/main.tex:442`)
supplies the conversion: "We assumed a deflator of 1.46 based on price indices
for construction costs as specified in [Hornykewycz et al. 2025]. In 2023
prices the shock amounts to about 58 bn EUR."

The deflation is also applied in the code lineage: the legacy `impulse_shock`
divides the whole impulse table by a scalar `inflator` before use (main
`src/util.jl`, revise/revisefinal `cbase2/src/core/util.jl`, and still present
in the current kernel — `src/core/diagnostics.jl`, exercised by
`tests/test_model.jl` — with `const inflator = 1.46` in `src/BeyondHulten.jl`).

Three-point verification against `impulses.csv`:

| Figure | Value | Price basis | Match |
| --- | ---: | --- | --- |
| CSV 2024 row | 80.76 bn | 2023 | the study's first-year 81 bn (sectoral disaggregation, rounding) |
| CSV horizon mean | 58.89 bn | 2023 | the study's "58 bn per year" |
| CSV mean / 1.46 | **40.338 bn** | 2019 | `paper/main.tex:516`'s exact "40.337 bn", to 0.002 % |
| CSV 2024 row / 1.46 | 55.32 bn | 2019 | first-year peak, deflated |

Hence the series is denominated in 2023 prices, and 40.3 bn is its 2019-price
level. The horizon average in 2019 prices **is** 40.3 bn, so "keeping the 40"
and "using the average" are the same statement. The two 2026-09-15 candidates
were coincidences, not explanations: the 2046 row (40.372) and the halving
(40.380) both merely sit near the true 2019-price average (40.338).

The deflator is applied **nowhere else in the pipeline**: the harness conveys
only the share vector psi from `impulses.csv` (unit-free, price-basis blind),
and G0 is hard-set in the pinned design; the level basis lived exclusively in
the manuscript footnote, which is why the factor-2 appeared unresolvable. No
code change is required.

Three consequences. (i) `impulses.csv` is in **current prices** (the source's
2023 price base): its horizon mean, 58.89 bn, is the paper's "58 bn in 2023
prices", and the design's `total_eur_m = 40300.0` is that mean at 2019 prices.
(ii) The factor-2 flag is **peak year vs horizon average** (80.76 vs 58.89)
times the **price base** (the 2024 row is 55.3 bn at 2019 prices), not the
halving/funding-split conjecture; the 80.76/2 = 40.38 coincidence is
superseded. (iii) The option-B comparison below mixes price bases: at 2019
prices the mean is 40.34 bn (≈ 1.33 % of GDP), i.e. the current G0 —
re-anchoring to the raw current-price mean would double-count the deflation.

## Decision

**Keep G0 = 40,300 EUR m — the deflated horizon mean (40.34 bn; the 0.09 %
difference is rounding to the paper's 40.3) — as the principal programme for the
executed generations and the manuscript.** Accepted by the operator on
2026-09-20; no re-anchor, no re-run, no re-mint. The manuscript narrative states
the programme's scale from the source study: an average annual impulse of
40.3 bn in 2019 prices (= 58.9 bn in 2023 prices), a first-year peak of about
81 bn (2023 prices, ~1.9 % of GDP), decaying to 0.3 % of GDP by 2050. The level
is unchanged from the executed `matrix_5x3_v9`/`v10` designs, so their
provenance is untouched. The system is homogeneous in the programme level —
shares, signs and rankings are unaffected by G0 and every relative metric
scales with it to first order (the sectoral and capacity variants are
nonlinear, so the scaling is approximate; option B' would measure it). The
**price-basis item is closed** (above): the raw impulse is in current (2023)
prices and the design's G0 is that horizon mean at 2019 prices.

## Options considered, and their resolution

| Option | Content | Resolution |
| --- | --- | --- |
| A. Narrative-only | The manuscript states the average annual impulse **at 2019 prices** (40.34 bn ≈ 1.33 % of GDP, 2024-2050) as the programme's time-profile — the same number as the model's G0, so no reconciliation is needed; the source's current-price mean (58.89 bn) is quoted as the price-base statement | **Chosen** — justified ex post by the price-basis identity: the 2019-price horizon average IS 40.3 bn; one paragraph, no runs |
| B. Re-anchor to the mean | **A no-op at the correct price base**: the mean at 2019 prices is 40.34 bn ≈ G0. Re-anchoring to the raw current-price mean (58,893) would mix price bases and double-count the 1.46 deflation | Not adopted as stated; a new generation only if a real level change is chosen |
| B'. Level-robustness batch | Optional: a handful of headline cells at G0 = 40.34 bn (the deflated mean, +0.1 % on G0) to demonstrate level invariance numerically | Optional; not required given linearity and the identity (~8 cells) |
| C. Re-read the source | Determine the price basis of `impulses.csv` | **Executed (in-repo, 2026-09-20)**: current (2023) prices per Hornykewycz et al. (2025); the halving conjecture resolved as a coincidence (80.76/2 = 40.38 merely sits near 40.338). What remains is the **composition** (Open items) |

Resolution (2026-09-20): **option A's form with the current G0** — the
manuscript states the average annual impulse at 2019 prices (40.34 bn) as the
principal programme, which is the model's G0; option B is not adopted (a no-op
at the correct price base); B' is optional and not required; C is answered in
the price-basis section.

## Open items

- **Principal programme: settled 2026-09-20.** The deflated horizon mean —
  G0 = 40,300 EUR m (40.34 bn to 0.1 %) — is the programme used by every
  generation and by the manuscript. Price basis: current (2023) prices,
  deflator 1.46; evidence above.
- **The halving/funding-split conjecture is superseded** by the deflator
  explanation (80.76 / 2 = 40.38 is a numerical coincidence). The
  funding-split question is still unread against the source, but it no longer
  bears on the level.
- **Composition.** psi is the raw 2024 current-price composition; the legacy
  scalar deflator could not deflate a composition either, so both pipelines
  share the limitation. A constant-2019-price incidence vector needs sectoral
  deflators and would change the incidence results (63.4 % specialised
  construction; the rigid-group rule) — a modelling choice, not housekeeping.
- **Time aggregation.** The legacy code used the first two years (2024-2025)
  for the demand multiplier and the horizon mean for the autonomous component;
  the current pipeline uses the horizon mean as the level, matching the paper's
  "yearly costs". Record in the manuscript.
- If a level change is ever chosen, which average (horizon mean, 2025-2050
  mean, mean vector with its 2050 composition drift) must be fixed in the
  design.
- **Manuscript wording.** State the principal programme as the average annual
  impulse at 2019 prices (40.34 bn ≈ 1.33 % of GDP); the source's current-price
  mean (58.89 bn) may be quoted as the price-base statement. (Part of the
  manuscript pass, open item 5 of `docs/CONCISE_SUMMARY.md`.)

## Consequences

- No `src/`, registry, design or preregistration change: the accepted level
  equals the one already executed, so the `matrix_5x3_v9`/`v10` generations
  keep their provenance. The composition item stays visible in
  `docs/CONCISE_SUMMARY.md` (open item 3, v3) and in this ADR; no DE record
  (the question is not abandoned).
- The price base is recorded in this ADR and in `experiments/README.md`
  (§Cell construction); the design files are not edited (their SHA-256 is
  pinned by `registry/preregistration.toml`), so the deflation is documented,
  not re-applied in code.
- Manuscript and paper-table passes may cite the study's trajectory and both
  price bases (2019: 40.3 bn average, 55.3 bn first-year peak; 2023: 58.9 bn
  average, 80.8 bn first-year peak).
- Any future G0 change requires a new design file and preregistration
  (ADR-0006) and, because it changes the programme level of every cell, a new
  generation under ADR-0004.

## Enforcement

- No design file's `[programme]` block is edited: the accepted level is the
  executed one (G0 = 40,300 EUR m); any future level change is a new design,
  preregistration (ADR-0006) and generation (ADR-0004).
- Wherever "40.3 bn" is stated, the price basis (2019 prices, deflator 1.46 per
  Hornykewycz et al. 2025) is cited.
- The price base lives in this ADR and in `experiments/README.md`; a future
  generation that wants the deflation explicit in the harness requires the
  same route.
- The resolution is recorded here and in `docs/CONCISE_SUMMARY.md`; this ADR
  is superseded only if the level or the price basis moves.
