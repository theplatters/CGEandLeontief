# ADR-0011 — Conform the frozen cbase2 notebook 01 to the canonical kernel

- **Status:** accepted
- **Date:** 2026-09-17
- **Supersedes:** —
- **Related:** ADR-0001 (amended: one read-only zone file is edited and the
  cbase2 freeze record re-pinned), ADR-0005 (cbase2 backport), ADR-0010
  (retained-sector rebuild), `registry/freeze.toml`, `cbase2/FROZEN.md`

## Context

`cbase2/01_data_wrangling.ipynb` is the frozen submission record of the §4.1
accounting transformation (Steps 0–3): it loads the raw German 2019 use table,
separates domestic from imported uses, forms $\Omega^{D}$ and $\Omega^{raw}$,
decomposes value added into its four components, and writes the seven
`cbase2/data_processed/` artifacts that `02_accounting_consistency.ipynb`
consumes.

Three of its constructs no longer match the codebase:

1. **The schema was declared as hard-coded table positions** (`N = 71`,
   `SEC = 2:72`, `FD = 75:81`). The canonical implementation locates the final-
   demand block and the import / product-tax rows **by label**
   (`src/core/accounting.jl`: `final_demand_columns`, `_FD_COLUMN_NAMES`). The
   label form is not cosmetic: after `retained_io_table`/`retained_dataset`
   (ADR-0010) the sector block is shorter and every position shifts, so
   positional indexing silently reads the wrong columns.
2. **The diagnostic figure was drawn with `Plots.jl`**, which is not a
   dependency of the `BeyondHulten` project (`Project.toml` deps; plotting is
   the GLMakie package extension, `ext/BeyondHultenGLMakieExt.jl`). The cell
   therefore cannot run in the project environment at all.
3. **Nothing connected the notebook to the canonical implementation**, so the
   record could silently drift from the kernel it documents.

The notebook is still the submission's self-contained wrangling record, and its
seven artifacts are **tracked** files (46 tracked files in the zone) consumed by
notebook 02 — so any conformance edit must preserve that contract exactly.

## Decision

- `cbase2/01_data_wrangling.ipynb` is conformed **in place**, as the single
  authorised exception to the `cbase2` freeze, with exactly three amendments:
  - the final-demand block and the import / product-tax rows are resolved by
    label, with the sector range derived from `N`;
  - the figure becomes an unconditional numeric summary plus an **opt-in**
    figure (`CBASE2_PLOT_FIGURE=1`, GLMakie) — the notebook runs headless;
  - an **optional, non-fatal** cross-check cell compares the notebook's own
    objects against `BeyondHulten.generate_data`.
- **The artifact contract is untouched.** Every `data_processed/` file keeps its
  name and content, no other file is written into the zone, and the cross-check
  never blocks the writes. Regeneration on Julia 1.12.7 reproduces all seven
  tracked artifacts **byte-identically**.
- The notebook stays a **self-contained record, not a client of the kernel**.
  The cross-check is guarded (`BeyondHulten` loaded only inside a `try`), so the
  notebook also runs against an environment without the package; a kernel that
  has moved away from the record reports `DIFF` loudly instead of failing.
- The `cbase2` freeze record is **re-pinned** to the conforming commit: new
  `recorded_commit`, `tree_hash`, `tracked_files` and `last_touch_*`, with the
  original freeze commit and tree hash preserved as
  `original_freeze_commit` / `original_freeze_tree_hash`. `cbase2/FROZEN.md`
  carries the amendment.
- The rest of the zone is unchanged: `02_accounting_consistency.ipynb`,
  `03_financing_closures.ipynb`, `cbase2/src/`, `cbase2/plots/`, the raw data
  and the prose records stay read-only.

## Consequences

- The submission record executes in the project environment and demonstrates
  its own equivalence to the canonical kernel — measured at the conformance
  date: $\Delta = 0.0$ on $\Omega^{D}$, $\Omega^{raw}$, `factor_share`, the wage
  component and gross value added.
- The two traps the pipeline history records are now **measured in-notebook**
  rather than described: the off-by-one import read (62'539 EUR m on this table,
  ≈9 % of total imported intermediates) and the row-vector normalisation
  (share rows up to 2.8e15 against a correct maximum of 1.0).
- The re-pinned baseline means the zone's recorded "frozen at" point is the
  conformance commit; the original freeze point remains visible above, in
  `registry/freeze.toml` and in git history.
- Notebooks 02 and 03 still carry the old positional schema and are *not*
  covered by this ADR: conforming them, or any further cbase2 edit, needs a new
  ADR.
- The committed `plots/01_import_intensity.png` is unchanged; regenerating it
  now requires a plotting backend (`CBASE2_PLOT_FIGURE=1`), which is stated in
  the notebook rather than silently dropped.

## Enforcement

- `registry/freeze.toml` `[frozen.cbase2]` with the re-pinned `recorded_commit`
  and `tree_hash`; `scripts/status.jl` warns when the zone's tree hash or
  working tree differs from that record and must render **0 warnings**.
- `cbase2/FROZEN.md` names the amendment, the ADR and the original freeze point.
- `scripts/check_repo.jl` must stay at 0 violations, and the seven tracked
  `data_processed/` artifacts must remain unmodified by a regeneration run.