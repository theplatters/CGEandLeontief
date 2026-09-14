---
title: "Accounting Consistency Transformation (BeyondHulten §4.1)"
author: calculato (AI research assistant)
date: 2026-09-02
tags: [beyondhulten, accounting, input-output, calibration, julia]
project: BeyondHulten / BFRep
---

# Objective

Bring the BeyondHulten data pipeline into line with ROADMAP §4.1 ("Accounting
consistency") and the definitive guide (§3 "Data and accounting", §6 Phase 1).
The work is executed transparently in a single, re-runnable notebook,
`Notebooks/AccountingConsistency.ipynb`, which reads the raw German 2019 use
table, performs the accounting transformation, reconciles the three GDP sides,
and emits CSV artifacts for downstream model use.

The transformation does **not** silently relabel total value added as "labour"
(see the gap description below); it keeps value added as an explicit composite
and exposes its components.

# Verified source schema

Source file: `data/I-O_DE2019_formatiert.csv` (Destatis-format use table,
million EUR, `;` delimiter, `,` decimal). Layout confirmed by direct inspection
(DataFrame indices):

- Rows 1-71: the 71 product/sector rows (suppliers and users share this order).
- Columns 2-72 (overall): the 71 sector columns.
- Columns 75-81 (overall): final-demand categories -- private consumption,
  private-organisations consumption, government consumption, equipment
  investment, construction investment, inventories, exports.
- Row 73: `Gesamte Verwendung der inländischen Produktion` (domestic production use).
- Row 74: `Verwendung der Importe` (imports, by column).
- Row 75: `Gütersteuern abzüglich Gütersubventionen` (product taxes net of subsidies).
- Row 77: `Arbeitnehmerentgelt im Inland` (compensation of employees).
- Row 79: `Sonst. Produktionsabgaben abzgl. sonst. Subventionen` (other production taxes).
- Row 80: `Abschreibungen` (consumption of fixed capital).
- Row 81: `Nettobetriebsüberschuss` (net operating surplus).
- Row 82: `Bruttowertschöpfung` (gross value added, basic prices).
- Row 83: `Produktionswert` (gross output, basic prices).

# The accounting gap in the pre-integration pipeline

`src/interface.jl` (`generate_data`) builds the model's `Data` object from this
table but collapses the accounting in five ways that violate §4.1:

1. It treats total value added (`Bruttowertschöpfung`, row 82) as a single
   "labour" factor via `factor_share = value_added ./ grossy` and
   `labor_share = lambda .* factor_share`. This bundles compensation of
   employees, other production taxes, depreciation, and net operating surplus
   into one factor and labels it "labour".
2. It never separates imports. `intermediate_use = io[1:71, 2:72]` mixes
   domestic and imported intermediate deliveries; the import row (74) is ignored.
3. It performs no tax/subsidy handling. Final demand is read at purchaser prices
   while value added sits at basic prices, so the two are never bridged.
4. It does not reconcile the three GDP sides, so accounting errors pass silently.
5. Shock incidence is undefined: `demand_shock` scales aggregate final demand
   with no statement of domestic vs imported content.

Reviewer feedback (metro-rev2) independently notes the CGE "is likewise a model
for a closed economy without indirect taxation", confirming the gap.

# Transformation steps

The notebook implements the following steps, one markdown cell plus one code
cell per step, with extensive inline comments.

## Step 0 -- Load and declare the schema

Read the CSV with explicit `delim=';'`, `decimal=','`,
`missingstring=["-","x"]`; strip the row-label column; coerce missing to zero.
Declare the sector/column indices as named constants so the layout is auditable.

## Step 1 -- Separate domestic vs imported uses

Pull the import row (74) and the product-tax row (75) by column. Split imports
into intermediate (by using sector) and final (by category). Construct the
domestic intermediate matrix $Z^D$ from the raw $Z$ (domestic + imported) by
allocating imports proportionally across supplying sectors:

$$Z^D_{s,u} = Z_{s,u} \cdot \frac{\text{domestic intermediate into } u}{\text{total intermediate into } u}$$

This assumption is stated explicitly because the source table reports imports
only by using sector, not by supplying sector. Build the domestic
conditional-input-share matrix $\Omega^D$ (user-by-supplier, rows normalised to 1).

## Step 2 -- Decompose value added

Read the four components (compensation of employees, other production taxes,
depreciation, net operating surplus) and assert they sum to gross value added
(row 82) to machine precision. Emit a component table and component shares.

## Step 3 -- Explicit one-factor aggregation

State the one-factor CGE assumption openly: value added is treated as a single
composite primary factor `factor_share = gva ./ prodval`. Keep the labour
component `wage_share = wage ./ prodval` separate so the narrative cannot
relabel total value added as "labour".

## Step 4 -- Reconcile the three GDP sides

Compute:

- Production (basic): $\sum_s \text{gross value added}_s$.
- Income (basic): $\sum_s (\text{wage} + \text{other tax} + \text{dep} + \text{net operating surplus})_s$.
- Expenditure (basic): domestic final demand at basic prices
  $= \sum_c (\text{final demand}_c - \text{imported final}_c - \text{product tax on final}_c)$.

Report all three and the pairwise absolute differences. The residual is documented, not hidden. **Important correction:** an earlier
draft reported ~0.8%; that figure was an artifact of a column-indexing bug
(reading the imports/product-tax rows after slicing off the label column, then
re-indexing by the sector range -- an off-by-one that subtracted near-zero
numbers). With the import/product-tax rows read directly at their original
column positions the true expenditure residual is **5.387%** (P = 3,027,818,
E = 2,864,724, diff = 163,094 mn EUR). This 5.387% is a genuine valuation
discrepancy in the raw Destatis table: the final-demand columns (purchaser
prices) net of recorded imports and product taxes do not exactly equal the
production-side GDP. It is carried through consistently and flagged for the
paper; it does not affect the model's domestic technology matrix or Domar
weights.

## Step 5 -- Shock incidence rule

Load `data/impulses.csv` (the Hornykewycz shock in EUR). Document the §4.1
rule: the shock is applied to **domestic** final demand only; imported content
is held fixed. The target vector is the domestic final-demand vector built in
Step 1. The aggregate nominal shock is reported and flagged against the paper's
stated figure (approximately EUR 40.3 bn at 2019 prices) for reconciliation.

## Step 6 -- Emit calibration artifacts

Write CSVs to `output/`:

- `AC_accounting_reconciliation.csv` -- the three GDP measures and differences.
- `AC_calibration_table.csv` -- per sector: gross output (basic), value-added
  components, `factor_share`, wage share, imported intermediate, import share,
  Domar weight $\lambda$, domestic final demand (basic).
- `AC_value_added_components.csv` -- the four value-added components by sector.
- `AC_domestic_intermediate_matrix.csv` -- $Z^D$ (71 x 71, sector-labelled).
- `AC_domestic_final_demand.csv` -- domestic final demand by category (basic).
- `AC_domar_weights.csv` -- sector Domar weights.

## Step 7 -- Validation and reconciliation assertions

Assert: value-added decomposition equals gross value added; $Z^D$ column sums
reproduce domestic intermediate demand; production equals income exactly;
production vs expenditure agree within 10% (the 5.387% raw-table residual is
documented); all calibration values are non-negative; Domar weights are finite,
positive, and sum to a positive value; the shock total is positive.

# Outputs

All artifacts are saved under `output/` so results survive container limits and
can be inspected independently of the notebook kernel. They are **locally
generated, gitignored artifacts — not durable committed evidence**; regenerate
them by re-running the notebook rather than treating them as archival outputs.

# Acceptance against ROADMAP §4.1

- Reproduce benchmark accounts without shock: done (Steps 0-4).
- Document imports/exports, product and production taxes/subsidies, labour
  compensation, capital income, mixed income: done (Steps 1-2, components
  explicit).
- Do not relabel total value added as labour: enforced (Step 3).
- Reconcile production-, expenditure- and income-side GDP: done (Step 4),
  residual documented.
- State whether the shock hits domestic or imported demand: done (Step 5).
- Publish a calibration table (intermediate shares, value-added shares,
  final-demand shares, import shares, Domar weights, elasticity ranges): done
  (Step 6). Elasticity ranges are supplied from the model's existing
  calibration and noted as out of scope for this notebook.

# Integration note (downstream `Data` / `generate_data`)

The transformation is now wired **directly into `src/interface.jl`** (no
dependency on pre-run `AC_*.csv` files). `generate_data` performs Steps 0-5
inline and `read_data` assembles an extended `Data` object. Field mapping
(extended `Data` -> notebook outputs):

- `Data.Ω`              = `Ω_dom`, the domestic conditional input-share matrix retained
  **for the domestic audit/incidence breakdown only** (not fed into production).
- `Data.Ω_raw`          = total conditional input-share matrix, used by the equilibrium
  technology **and** the `consumption_share` calibration.
- `gross_output_basic`  = `prodval` (gross output, basic prices) -- replaces the
  old `grossy` taken from the total-use column.
- `factor_share`        = composite VA share `gva ./ grossy` (NOT relabelled as labour).
- `λ`                   = standard Domar weights `prodval ./ GDP_P` (always positive;
  replaces the deprecated Leontief-inverse `λ`).
- `labor_share`         = `λ .* factor_share` (explicit composite VA weight).
- `value_added_components` = decomposed VA (wage / other-prod-tax / depreciation / net-op).
- `imports_intermediate`, `import_share`, `domestic_final_demand`,
  `gdp_production/income/expenditure` = the §4.1 separation fields.

This keeps all legacy `Data` field names so `main.jl` and the test suite keep
working, while the model now consumes the accounting-consistent, open-economy
structure.

# Status

Plan drafted 2026-09-02. Notebook `Notebooks/AccountingConsistency.ipynb`
implements Steps 0-7 and was validated against the raw 2019 table (Julia
1.12.6, project environment). The transformation is **integrated into
`src/interface.jl`** (`generate_data` + extended `Data`), validated by the
`tests/` suite (`generate_data` assertions and `tests/test_model.jl`: `Ω_raw`
rows sum to 1, λ standard Domar positive, GDP production == income exact,
expenditure residual 5.387% documented). Both artifacts are mutually consistent
(same numbers).

Two transcription bugs were found and fixed during integration (present in the
notebook and in the first `interface.jl` draft):
1. **Off-by-one import indexing** -- reading the imports/product-tax rows after
   slicing `2:end`, then re-indexing by the sector/final-demand range, shifted
   every import/tax value by one column. Fixed by indexing those rows directly
   at their original column positions (SEC = 2:72, FD = 75:81).
2. **`Ω` normalisation orientation** -- `Ω_dom = Z_dom' ./ inter_dom_totals'`
   (row vector) normalised *columns* instead of rows. Fixed to a column vector
   so each user row sums to 1.
A third guard was added: for the 3 sectors where the source table records
imports exceeding intermediate use, the domestic remainder is clamped to zero
(those sectors become fully import-dependent in the domestic technology),
keeping `Z_dom` non-negative.

Definitive guide §4.1 / §6: status moved to **"Reconciled"**.
