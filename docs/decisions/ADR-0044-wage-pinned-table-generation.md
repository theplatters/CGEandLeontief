# ADR-0044 — Regenerate tables on OpenCGE v0.5.0

- **Status:** accepted
- **Date:** 2026-10-23
- **Supersedes:** the GAMMA/DELTA real-wage wording in ADR-0002; historical run records remain immutable.
- **Related:** ADR-0028, ADR-0032, ADR-0037, OpenCGE@ADR-0043; `registry/closures.toml`

## Context

OpenCGE tag `v0.5.0` restores the square nominal-wage system for GAMMA and
DELTA: `w = 1`, `F = 0`, endogenous employment, and moving CPI and real wage.
The extra real-wage restriction in the intervening kernel line is removed.
The release's Project.toml still declares version `0.4.1`; the tag and tree
hash, rather than that metadata alone, identify the kernel generation.

## Decision

Pin the published `v0.5.0` tag through Pkg. Re-execute the table sources with
the same calibration, programme and parameters under new design/run ids:
`matrix_5x3_v11`, `pricing_capacity_rows_v2`, `pricing_capacity_v3`,
`pricing_mobile_v3`, `pricing_markup_v3`, and `pricing_profit_v3`.
The capacity-rows design has one cell per specification: the earlier failed
BF-F1-dm05 attempt and its replacement remain in the historical generation.
Preregister and commit the designs before executing through `experiments/run.jl`.

Regenerate Tables 4–7 and Appendix E only from these runs. Retain the
consumption-unit price convention of ADR-0037 and the welfare definitions;
change GAMMA's table label to fixed nominal wage. Retain the historical
manifests, baseline and failed runs without modification. Do not mix pins
in the new tables or claim a new closure validation status without evidence.

## Consequences

Registry formulations follow the released nominal-wage equations. BETA's
labour supply remains a function of the real wage. New run evidence is
recorded in the scenario registry and generated status board, rather than
restated here. Existing user plotting changes are preserved.

## Enforcement

Run the repository gate before and after each batch, regenerate the board
with zero warnings, verify table inputs against accepted manifests and the
release pin, run `julia --project=. tests/runtests.jl`, and build the revision.
