---
title: "Process Comments: Observations by Pipeline Chapter"
author: "Hermes Agent (Lt. Cmdr Data), for Prof. Dr. J. Kapeller"
date: "2026-09-15"
project: "BFRep (3)BeyondHulten / Metroeconomica revision"
tags: [cbase2, process-comments, pitfalls, data-audit]
version: 1
last-updated: "2026-09-15"
---

This file carries the *commentary*: findings, corrections, deviations, and
open questions, organized in chapters corresponding to the pipeline
notebooks. `documentation.md` stays a neutral, timeless description of the
pipeline as it stands; everything here is dated observation and may be
superseded as the pipeline evolves. Entries are additive; nothing is
silently rewritten -- later entries state what they replace.

# Environment (container)

- **2026-09-15.** The container's Julia setup had drifted: the image was
  built while the `release` channel was 1.13, so the baked depot carried an
  environment named `v1.13` and a kernel named `julia-1.13`, while the
  running binary is **1.12.7** (resolved through the sandbox-home juliaup,
  not the baked install). The 1.13-format Manifest also lost
  `git-tree-sha1` fields on some entries, so packages loaded as "not
  installed". Resolution: depot environment renamed to `v1.12` and
  regenerated natively (`Pkg.update()` + full precompile), stale `v1.11`
  environment and ~3.1 GB of unused sandbox-home depot data removed, kernel
  re-registered as `julia-1.12`, build ARG pinned to 1.12.7 in the image
  Dockerfile. Verified: plain `julia` loads the full stack (CSV 0.10.17,
  DataFrames 1.8.2, Plots 1.41.7, NonlinearSolve 4.32.0).

# Notebook 01 -- Data Wrangling

- **2026-09-15. Finding: silent NaN rows in the domestic share matrix.**
  Three sectors record imported intermediates exceeding their total
  intermediate use; after the domestic remainder is clamped, their
  $\Omega^{D}$ rows are 0/0 = NaN. The parent notebook carried these NaNs
  without notice (it asserted only the $Z^{D}$ column sums). Decision: such
  rows are encoded explicitly as **all-zero** ("no domestic content"),
  consistent with the fully-import-dependent interpretation, and the
  row-sum conventions ($\sum = 1$, or $0$ for these rows) are asserted.
- **2026-09-15. Clamps made audible.** The negative-remainder clamp in the
  domestic matrix and the zero-domestic-fraction clamp (final-demand
  bridge, notebook 02) previously passed silently; both now print a
  warning with the affected count/categories.
- **2026-09-15. Standing assumption (carried, documented).** The
  proportional allocation of imports across supplying sectors remains a
  stated assumption: the source table reports imports by *using* sector
  only. Any future Armington-type extension re-opens this.
- **2026-09-15. Porting note.** `DataFrame(A, :auto)` accepts a matrix
  directly; the intermediate `Matrix(A, :auto)` constructor form is not a
  valid call in the current stack (recorded because it cost one
  validation round).

# Notebook 02 -- Accounting Consistency

- **2026-09-15. Correction: residual documentation.** The parent's
  validation section stated the production-vs-expenditure residual as
  "< 1.5%" while the actual value in this data is **5.387%** and the hard
  gate is 10%. The text now states the real figure. The residual itself is
  a raw-table valuation discrepancy (documented, not fixed).
- **2026-09-15. Validation added.** The parent's bridge computed
  `cons_vec`/`cons_share` (the `consumption_share` calibration) without
  any checks; these are now asserted finite, non-negative, and summing to
  one. $\Omega^{raw}$ row sums are re-checked at the bridge.
- **2026-09-15. Open question: shock magnitude.** The defensive load of
  `impulses.csv` (27 rows x 75 columns) sums row 1's numeric columns to
  EUR 80.76 bn, while the paper states ~EUR 40.3 bn at 2019 prices. The
  factor ~2 discrepancy is unresolved: candidates are the row-1 year not
  being 2019, category double-counting in the column sum, or a units
  difference. To be checked before the financing notebooks rely on the
  impulse.
- **2026-09-15. Convention carried.** The expenditure side sums domestic
  final demand **including exports** (the parent pipeline's documented
  convention); it is part of what the 5.387% residual absorbs and is kept
  unchanged for comparability.

# Notebooks 03--08 (financing, closures, preregistration, validation, matrix, sensitivity)

- No observations yet; chapters fill as the notebooks are built. Expected
  first entry for 03: the F1/F2/F3 equation design decisions and their
  rationale, recorded before implementation.
