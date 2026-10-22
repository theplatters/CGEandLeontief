# ADR-0037 — Consumption-unit price reporting

- **Status:** accepted
- **Date:** 2026-10-22
- **Supersedes:** ADR-0032's raw-price reporting convention for the pricing
  tables. Stable ids, kernel equations and recorded runs are unchanged.
- **Related:** ADR-0004, ADR-0018, ADR-0028, ADR-0029, ADR-0031, ADR-0032;
  `registry/closures.toml` `[pricing.P1]`,
  `[pricing.P2]`; `experiments/make_pricing_prices.jl`;
  `tests/test_pricing_reporting.jl`

## Context

The recorded mobile/sectoral pricing runs use CPI equal to one; the GAMMA
and DELTA runs use the common wage equal to one and the net external
transfer equal to zero. The demand-only single-wage branch has prices and
CPI equal to one in both normalizations, so the matrix is unaffected. With
the pricing doors open, raw GDP deflators and maximum sectoral deviations
are expressed in different units across these families. Tables 6, 7 and E
previously copied these raw metrics side by side. Their signs and rankings
were therefore not comparisons in consumption units.

## Decision

1. **Report all pricing rows in consumption units.** With `C = CPI(p)`
   and baseline `C0 = 1`, report `p_hat = p/C`, GDP price index
   `D_hat = D/(C/C0)`, `Delta P = 100*(D_hat-1)` and maximum sectoral
   deviation `100*max(abs(p_hat-1))`. `D` is the recorded Törnqvist GDP
   deflator. The CPI uses calibrated consumption shares and the cell's
   consumption elasticity, including DELTA's Leontief limit, matching the
   model's wage-deflation convention. This is a common consumption-unit
   convention with model-specific cost-of-living indices, not an invariant
   physical basket or a monetary inflation prediction.
2. **Keep economics and provenance.** Quantities, employment, consumption
   indices and real GDP are unchanged. Nominal income and the GDP price
   index must be rescaled together; never recompute GDP by combining raw
   income with the normalized deflator. Manifests, index rows and design
   preregistrations remain historical records. There is no kernel edit,
   pin bump or new run generation.
3. **Generate derived prices reproducibly.** `make_pricing_prices.jl`
   reads the saved raw `solution.csv` artifacts and manifests. It verifies
   the reference prices, data hash, sector ordering, recorded cell and
   pricing parameters, and raw maximum deviations. For every fixed-wage
   cell it reconstructs the GDP deflator and verifies its homogeneity and
   real-GDP invariance. The output `revision/tables/pricing_prices.csv`
   records source hashes. Historical whole-design hashes can differ after
   the ADR-0032 BF replacement cell was appended; those hashes are retained,
   and the per-cell parameters are checked against the current design.
4. **Table generation fails on stale inputs.** `make_pricing_door.py`
   consumes the derived price columns and checks the manifest and current
   design hashes. It checks local solution and IO hashes when those files
   exist, allowing table regeneration from the committed sidecar in a
   clean clone. Real-outcome columns still come from their existing sources.
## Consequences

The pricing-table narratives and Appendix D compare GDP prices relative to
consumption prices. The former all-sector sign-reversal explanation in wage
units cannot be carried over unchanged. Maximum sectoral price deviations
must also be normalized. The demand-only matrix numbers remain unchanged.
See the derived sidecar and its cited run ids for the corrected values.

## Enforcement

- `julia --project=. experiments/make_pricing_prices.jl --check`
- `python3 revision/tables/make_pricing_door.py`
- `julia --project=. tests/runtests.jl` (nominal-scale invariance,
  reconstruction of available saved pricing solutions, and stale-source
  rejection by the table generator)
- `julia --project=. scripts/status.jl` and `scripts/check_repo.jl`
- `latexmk -pdf` in `revision/` (a clean temporary output directory may
  be used to avoid stale local auxiliary files)
