# ADR-0036 — Print-ready matrix panels (OpenCGE v0.3.2) without a new generation

- **Status:** accepted
- **Date:** 2026-10-22
- **Supersedes:** —
- **Related:** ADR-0004 (runs are immutable manifests; paper figures cite
  `run_id`s), ADR-0025 (matrix visualisation: validated re-solve, figures
  are not runs), ADR-0028 (pin bumps are generation events), ADR-0031
  (precedent: pin bump verified by `reproduce.jl --check`);
  `OpenCGE@ADR-0034`, `OpenCGE@ADR-0035`

## Context

The 15 manuscript panels (`revision/images/panel_5x3_<labour>_<financing>.png`,
revision Section 5.2) were drawn at 1900 × 900 px with 10–12 px text. At
about 14 cm print width that is 2–3 pt, and the figures also carried a
title, a three-line footnote block and a `Source: <run_id>` line. The
drawing code lives in the kernel's GLMakie extension, so changing it means a
kernel release and a pin bump. Under ADR-0028, a pin bump invalidates
the provenance of existing runs unless it is shown not to touch them.

OpenCGE v0.3.1 (`OpenCGE@ADR-0034`) added an opt-in `style = :manuscript`.
v0.3.2 (`OpenCGE@ADR-0035`) reworked it: numbered horizontal programme
bars, one shared legend, integer-percent ticks, and caller-supplied
`sector_labels`. v0.3.0 adds first-order welfare diagnostics. All three
releases add code; none changes the solver, closures, calibration or data
layer.

## Decision

1. **Pin.** The paper repo pins OpenCGE `v0.3.2` (compat `"0.3.2"`),
   moving from `v0.2.1`.
2. **No new generation.** The bump is presentation-only and additive, so
   no design is re-run. The recorded runs and their `v0.2.1` provenance
   stay the citation source. Evidence:
   - `scripts/reproduce.jl --check` passes on the new pin (15/15 within
     1e-8);
   - `experiments/plot_matrix.jl` re-solves all 15 headline cells with zero
     price, quantity and metric deltas against their recorded artifacts.

   A later kernel release that touches the solver, closures, calibration or
   data layer is a generation event again (ADR-0028).
3. **Manuscript panels.** `experiments/plot_matrix.jl --panel-dir` draws
   with `style = :manuscript` and the curated short names in
   `data/sector_short_names.csv` (`label,short`, keyed by
   `data/sector_names.txt`). The exploratory layout stays available in the
   kernel.
4. **Run ids.** The panels and their captions no longer show `run_id`s.
   ADR-0004 is met by this record and by the driver's validation report.
   This ADR records the mapping: `panel_5x3_<L>_<F>.png` is drawn from
   `matrix_5x3-v10-<L>-<F>`, except BETA, which is drawn from
   `matrix_5x3-v10-BETA-<F>-etas05` (ADR-0022 presentation rule).

## Consequences

- The panels are legible at print width, and the §5.2 figures cite no ids.
- Run manifests written after this bump record `v0.3.2` in `[provenance]`.
  Existing manifests keep `v0.2.1`. Tables mixing both pins remain valid
  only because the zero-delta evidence above holds.
- The kernel's first print-layout record was numbered ADR-0032, colliding
  with this repository's ADR-0032. It is renumbered `OpenCGE@ADR-0034`, so
  the next free shared id is ADR-0037.

## Enforcement

`scripts/reproduce.jl --check` and the `plot_matrix.jl` validation report
(15/15, zero deltas) are rerun before any figure regeneration on a new pin.
`scripts/check_repo.jl` and `tests/runtests.jl` must pass on the pinned
kernel.
