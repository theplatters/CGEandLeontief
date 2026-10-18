---
reviewed: "cf414f19419f6fa40bc04d0352326807e7a226cc"
---
# Törnqvist index

The Törnqvist index averages the value shares of the two periods being compared
instead of fixing base-period weights. The project note states it for
consumption quantities $c_{it}$ and prices $p_{it}$, with shares
$s_{it}=p_{it}c_{it}/\sum_j p_{jt}c_{jt}$, as
$\log(Y_1/Y_0)=\sum_i\frac{s_{i0}+s_{i1}}{2}\log(c_{i1}/c_{i0})$;[^1] in the
BeyondHulten paper it is the approximation to the Divisia quantity index of
final demand.[^2]

In the project's kernel the quantity form is the exported
`tornqvist_quantity_index` helper (`src/util.jl`) used by every model variant;
it replaced the fixed-base quantity sum that had previously been reported as
real GDP.[^3] The helper, like the rest of the kernel's `src/`, now lives in the
OpenCGE repository split out of this one in October 2026 and reaches the paper
repo only through the committed root `Manifest.toml` pin.[^6]

The index plays two roles that the consumption-only formula does not separate.
Built on household consumption it is a consumption/welfare aggregator, which is
what the legacy `real_gdp` accessor and `Solution.real_gdp` field still hold;
the kernel's newer API names it `real_consumption(sol)` and reports GDP through
`gdp_income`/`gdp_expenditure` with `gdp_deflator`.[^4] Real GDP is instead
income divided by a Törnqvist price index of the seven signed value components
$V=[C_{gross},\,G+programme,\,I,\,X,\,-M_{final},\,-M_{int},\,-T_{int}]$,
with unit values $|V_j|/Q_j$ and signed value shares $V_j/\sum V$, chain-linked
to the reference solution.[^4]

The BeyondHulten project note records that GDP was moved from a Laspeyres
quantity index to a Törnqvist index as part of the paper revision, and DE-0005
documents the fixed-base sum that the shared helper superseded.[^5]

[^1]: logseq://9a852f6d-9f0f-8f73-8989-ed225acae56b/page/pages/T%C3%B6rnqvist%20index.md — the concept note gives the two-period formula and the value-share definition.
[^2]: ../decisions/ADR-0018-real-gdp-measurement.md — its Related line cites `paper/main.tex` "Törnqvist approximation to the Divisia quantity index of final demand".
[^3]: ../dead-ends/DE-0005-fixed-base-gdp-sum.md — records the exported helper and the dead end it replaced.
[^4]: ../decisions/ADR-0018-real-gdp-measurement.md — Decision 1 defines the Törnqvist price-index deflator, the seven signed components and the chain-linking; Decision 4 names `real_consumption(sol)` and the `gdp_*` measurement API.
[^5]: logseq://9a852f6d-9f0f-8f73-8989-ed225acae56b/page/pages/Beyond%20Hulten.md — project state: "Changed the GDP from a Laspeyres quantity index Index to a Törnqvist index".
[^6]: ../decisions/ADR-0028-split-paper-and-opencge-kernel.md — the kernel (`src/`, `ext/` and kernel tests) became the OpenCGE repository and enters this repo only via the committed root `Manifest.toml` pin.
