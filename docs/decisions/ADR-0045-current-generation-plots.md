# ADR-0045 — Regenerate all matrix plots on OpenCGE v0.5.4

- **Status:** accepted
- **Date:** 2026-10-24
- **Supersedes:** ADR-0036's v10 source mapping for the appendix panels.
- **Related:** ADR-0025, ADR-0028, ADR-0036, ADR-0044

## Context

The main-text comparisons already select the ADR-0044 matrix generation,
but the detailed appendix panels retain their historical v10 mapping.
The author requests regeneration with OpenCGE `v0.5.4`. Comparing the
published tags shows only plot changes and accounting documentation:
thin black financing markers, optional comparison annotations, removal
of panel error bars and a two-scatter panel method. Solver, calibration
and accounting code are unchanged.

## Decision

Regenerate the six exploratory matrix figures, fifteen detailed panels and
the range comparison and combined scatter image from `matrix_5x3_v11`.
Pin OpenCGE
`v0.5.4` through Pkg, tree `228312a0668cd3df12b66d6537d9013c1602962f`.
Retain the numerical generation under ADR-0036's presentation-only
exception, conditional on the reproduction baseline passing at 1e-8 and
all fifteen plotted cells reproducing their recorded manifest metrics.
Map `panel_5x3_<L>_<F>.png` to `matrix_5x3-v11-<L>-<F>`, replacing BETA
with `matrix_5x3-v11-BETA-<F>-etas05`. Use the same headline selection for
the exploratory figures. The main-text range and scatter panels retain
their BF/BETA/DELTA selections documented in `experiments/README.md`.
Combine the F1 and F3 scatter selections through
`plot_closure_scatter_panel`, preserving the selected closures and sector
weights. Use one manuscript figure with (a) F1 and (b) F3, with separate
axis limits and unified closure and shock-status legends.
Optional annotations remain at the release default, `:none`.
The financing legend uses the short ids `F1`, `F2` and `F3`.

## Consequences

All regenerated matrix plots use the current run generation. Keep the
existing layout and absence of panel error bars. Use curated sector names
when the local label file is available, otherwise the documented fallback.
Historical run records and the reproduction baseline remain immutable.
No new numerical generation or table-source change is needed. Historical
manifests retain their original kernel pins. The release still declares
package version `0.5.4`; its tag and tree identify the installed code.
The local Julia 1.12.6 Pkg update also records its bundled JLL versions and
runtime metadata in the generated Manifest.

## Enforcement

Use `experiments/plot_matrix.jl` to re-solve and validate every selected
cell against its executed manifest under ADR-0025. Run
`scripts/reproduce.jl --check`, the repository gate and paper suite, and
record any pre-existing failures in the lab log. Never replace the baseline
to make the check pass.
