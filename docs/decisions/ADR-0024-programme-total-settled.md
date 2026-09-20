# ADR-0024 — Programme total: G0 = 40,300 kept (the 2019-price horizon average)

- **Status:** accepted (user ratification, 2026-09-20)
- **Date:** 2026-09-20
- **Supersedes:** —
- **Related:** ADR-0003 (registry single source of truth), ADR-0004 (runs),
  ADR-0006 (preregistration — G0 is pinned by the design file's SHA-256);
  `experiments/designs/matrix_5x3_v6.toml` and `matrix_5x3_v9.toml`
  (`[programme] total_eur_m = 40300.0`); `cbase2/data_raw/impulses.csv`
  (the 27 x 71 annual impulse, 2024-2050); `cbase2/process_comments.md`
  (the factor-2 flag, 2026-09-15); `paper/main.tex` (lines 442 and 516:
  the programme statement and its price basis); Hornykewycz, Kapeller,
  Weber et al. (2025), "Carbon neutrality in the residential sector: a
  general toolbox and the case of Germany", npj Climate Action 4(1) (also
  IFSO WP 41); `docs/log/2026-09.md` (the measurements); `docs/CONCISE_SUMMARY.md`
  (open item 3, left untouched per the operator, 2026-09-20).

## Context

The design's programme magnitude is `total_eur_m = 40300.0` — the paper's
"40.3 bn at 2019 prices" — and every executed cell of the matrix and the
sectoral family is pinned to it. The raw impulse table behind the
programme vector (`impulses.csv`, 27 annual rows 2024-2050, 71 VGR
sectors) contains no row that equals 40.3 bn: its first row sums to
80,760,315,517 (80.76 bn). The cbase2 session of 2026-09-15 flagged the
factor-2 gap with two candidates: the 2046 row (40,371,886,437) and an
exact halving of the 2024 row (40,380,315,779).

Measured on 2026-09-20 (full horizon): sum 1,590,122,645,111; **mean
58,893,431,300 (58.89 bn = 1.945 % of 2019 GDP at G0 = mean)**; median
65,873,835,779; peak 80,760,315,517 (2024); tail 18,081,686,729 (2050).
Composition is broadly stable against the 2024 vector (max |Delta-psi| =
0.092 over the horizon, worst year 2050, when machinery's share rises to
13.5 %).

## The price basis (resolved 2026-09-20)

The factor-2 was a price-basis mismatch, not a data error. The source
study states the programme in **2023 prices**: "an additional yearly
investment of 58 bn EUR is needed", "the total investment until 2050 sums
up to 3.1 trillion EUR (in 2023 prices)", and "in the first year of the
policy measure, additional costs of 81 bn EUR are anticipated, about 1.9 %
of GDP, which over time decreases to 0.3 %". The manuscript's own footnote
(`paper/main.tex:442`) supplies the conversion: "We assumed a deflator of
1.46 based on price indices for construction costs as specified in
[Hornykewycz et al. 2025]. In 2023 prices the shock amounts to about 58 bn
EUR."

Three-point verification against `impulses.csv`:

| Figure | Value | Price basis | Match |
| --- | ---: | --- | --- |
| CSV 2024 row | 80.76 bn | 2023 | the study's first-year 81 bn (sectoral disaggregation, rounding) |
| CSV horizon mean | 58.89 bn | 2023 | the study's "58 bn per year" |
| CSV mean / 1.46 | **40.338 bn** | 2019 | `paper/main.tex:516`'s exact "40.337 bn", to 0.002 % |
| CSV 2024 row / 1.46 | 55.32 bn | 2019 | first-year peak, deflated |

Hence the series is denominated in 2023 prices, and 40.3 bn is its
2019-price level. The horizon average in 2019 prices **is** 40.3 bn, so
"keeping the 40" and "using the average" are the same statement. The two
2026-09-15 candidates were coincidences, not explanations: the 2046 row
(40.372) and the halving (40.380) both merely sit near the true 2019-price
average (40.338).

The deflator is applied **nowhere in the codebase**: the pipeline conveys
only the share vector psi from `impulses.csv` (unit-free, price-basis
blind), and G0 is hard-set in the pinned design. The level basis lived
exclusively in the manuscript footnote, which is why the factor-2 appeared
unresolvable. No code change is required.

## Decision

**Keep G0 = 40,300 EUR m** — the paper's rounded "40.3 bn in 2019 prices"
(the exact 2019-price horizon average is 40.338 bn, 0.09 % above; the
difference is rounding, not substance). The manuscript narrative states the
programme's scale from the source study: an average annual impulse of
40.3 bn in 2019 prices (= 58.9 bn in 2023 prices), a first-year peak of
about 81 bn (2023 prices, ~1.9 % of GDP), decaying to 0.3 % of GDP by 2050;
all reported effects are for the 40.3 bn programme and scale linearly with
it (shares, signs and rankings unaffected by the level). No re-run, no
re-mint, no design change.

## Options considered, and their resolution

| Option | Content | Resolution |
| --- | --- | --- |
| A. Narrative-only | Manuscript states the programme per the source study; model keeps G0 = 40,300 | **Chosen** — justified ex post by the price-basis identity: the 2019-price horizon average IS 40.3 bn |
| B. Re-anchor to the mean | New generation at G0 = 58,893 (2023 prices) | Unnecessary — 58.89 bn is the same programme in 2023 prices; the matrix's G0 is its 2019-price level |
| B'. Level-robustness batch | A few headline cells at G0 = mean | Optional; not needed given linearity and the identity |
| C. Re-read the source | Determine the price basis of `impulses.csv` | **Executed** — 2023 prices per Hornykewycz et al. (2025); the halving conjecture resolved as coincidence (80.76/2 = 40.38 ≈ 40.338 by the shape of the decay profile, not by double counting) |

## Consequences

- G0 stays at 40,300; every executed cell keeps its provenance. The open
  item 3 of `docs/CONCISE_SUMMARY.md` is resolved by this ADR (the summary
  itself is left untouched per the operator, 2026-09-20).
- Manuscript and paper-table passes may cite the study's trajectory and
  the two price bases (2019: 40.3 bn average, 55.3 bn first-year peak;
  2023: 58.9 bn average, 80.8 bn first-year peak).

## Enforcement

- No `[programme]` edit in any design file; a future G0 change requires a
  new ADR, a new design and a new generation (ADR-0006, ADR-0004).
- Wherever "40.3 bn" is stated, the price basis (2019 prices, deflator
  1.46 per Hornykewycz et al. 2025) is cited.