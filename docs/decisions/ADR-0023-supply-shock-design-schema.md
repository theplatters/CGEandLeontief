# ADR-0023 — Supply-shock design schema for the identification arm

- **Status:** accepted (user ratification, 2026-09-20)
- **Date:** 2026-09-20
- **Supersedes:** —
- **Related:** ADR-0002 (naming), ADR-0006 (preregistration), ADR-0007
  (WIP and repository gate), ADR-0018 (real GDP measurement; the "schema
  item" recorded in its Consequences), ADR-0022 (Sectoral labour markets,
  whose Consequences record that "the supply-shock schema ADR must take a
  later number"); `docs/ETAs.md` (v3, the arm's workplan); DE-0011 (the
  numeraire is not the lever); `paper/equivalence.tex` (Lemma 1: a demand
  channel to real wages requires technology, i.e. `A != 1`);
  `experiments/probes/supply_identification.jl` and
  `experiments/probes/probe15_sectoral_supply.jl` (the read-only pilots);
  `experiments/run.jl` (the four hard-coded `Shocks(...)` call sites);
  `registry/closures.toml [labor.BETA]` open gates.

## Context

`experiments/run.jl` hard-codes the no-shock

    shocks = Shocks(ones(N), ones(N), zeros(N))

in `build_reference`, `build_cell_model` and `solve_cell` (lines 294, 311,
~508 and ~543 at HEAD), so no design can carry a supply-side scenario. The
arm that identifies `eta_s` cannot run without a per-cell supply-shock
source and magnitude from the design file.

Why the supply arm exists at all: under demand-only shocks the price block
is demand-free, the real wage sits at its anchor, and `eta_s` is
unidentified — the executed `matrix_5x3_v9` ladder is a sensitivity band,
not an estimate (`paper/equivalence.tex` Lemma 1; measured in
`matrix_5x3_v9`, and in the demand-only `BETA ≡ ALPHA` degeneration).
Technology is the only channel that moves `w/P` (DE-0011: the numeraire is
not the lever). ADR-0022 kept the supply arm as its identification
complement and deferred the schema decision to this ADR.

The arm changes the experiment, not the closure: nothing in `src/` moves,
so no executed generation's provenance is invalidated and no re-minting is
needed (contrast ADR-0022).

## Options

| Option | Mechanism | Assessment |
| --- | --- | --- |
| A. Per-design `[shock]` block | A `[shock]` table in the design file read by `run.jl`, with `kind = "none" | "sectoral" | "programme" | "uniform"`, a magnitude, and a target rule; default `kind = "none"` reproduces today's behaviour bit-for-bit | **Recommended.** One schema change, existing designs untouched, manifests self-describing |
| B. Per-cell shock fields | Shock source/magnitude repeated on every cell | Noisier manifests; higher risk of drift between cells that should be identical; no benefit over A because the arm's shock is a design-level property |
| C. One scenario TOML per shock design, convention only | 45 small tomls with no schema change | Multiplies design files and bypasses the manifest contract; the schema change is needed anyway for `tests/test_run_manifest.jl` |

## Decision

Adopt option A: the design schema gains a `[shock]` block, defaulting to
`kind = "none"`, so that every existing design file is unchanged and
behaves bit-identically.

1. **Schema.** `[shock]` in the design file:

   ```toml
   [shock]
   kind = "none"            # "none" | "sectoral" | "programme" | "uniform"
   magnitude = 1.2          # the shock level (for "sectoral": the A_i of the target)
   alpha = 0.10             # for "programme": A_i = 1 + alpha * psi_i on the programme's sectors
   targets = []             # for "sectoral": 1-based sector indices; empty + kind="programme"/"uniform" = rule-based
   note = ""
   ```

   `experiments/run.jl` builds `Shocks(A, ones(N), zeros(N))` with `A`
   resolved from the block; `kind = "none"` keeps `A = ones(N)`.
2. **Reference.** The reference continuation stays no-shock; real GDP and
   the welfare index remain "against the baseline" (the ETAs convention).
   A supply arm design therefore measures a *different experiment* (a
   technology shock) against the same baseline — never mixed into the
   demand-only matrix.
3. **Id convention.** Follows ADR-0002, with the shock id in the variant
   slot: `<design>-<labor>-<financing>[-<shockid>][-<etas>]`, e.g.
   `supply_etas_s1-BETA-F2-etas2`. One design file per shock id
   (`supply_etas_s1`, `supply_etas_prog`, `supply_etas_unif`) keeps the
   manifests readable (ETAs v3).
4. **Identification vehicle.** The arm's primary identification claim runs
   on the **scalar single-wage BETA** cells (2N+2, `eta_s` a scalar — the
   v9 matrix form), because the preregistered signature is scalar:
   `ln L / ln (w/P)` recovers the input `eta_s` to `1e-6`. Each design
   additionally carries the **sectoral uniform-vector** cells (3N+1,
   `eta_s_vec = fill(eta_s, N)`) as the robustness arm: there the
   identification is a *vector* (`eta_s,i = ln(L_i/Lbar_i)/ln(w_i/Pi(p))` with `Pi(p)` the CPI, per
   sector) and no scalar signature applies — reported as per-sector
   recovery, not as one number. The two rows are different systems and are
   never merged into one table.
5. **Manifest.** The manifest records the resolved shock vector (the
   `[shock]` fields plus the derived `A` vector) so that every cell is
   reproducible without re-reading the design.
6. **Canary and determinacy at `A != 1`** (the ETAs open gates, checked
   not assumed): the ADR-0010 external-account identity and the ADR-0017
   round-gain criterion have only been asserted at `p = 1`. They are
   asserted on the first cells of the first design as diagnostics and
   recorded; if either fails, a follow-up ADR is written rather than
   threshold-fitted. Note probe15 row set 3 has already measured the
   rigid-corner nesting bit-exact under `A[1] = 1.2`, and financing
   neutrality F2 = F3 extends to supply shocks (recorded in
   `registry/closures.toml`).

## Consequences

- Three preregistered designs, 45 cells each (5 labour rows x 3 financing
  columns x 3 shock magnitudes), a few seconds per cell, one batch per
  design under the WIP limit (ADR-0007). `supply_etas_s1` first (the
  pilot's design), then `_prog` (the economically relevant case: the
  programme's own sectors), then `_unif`.
- The arm delivers two publishable quantities: the identification theorem
  (the map `eta_s -> (w/P, L, real GDP)` is monotone and invertible under a
  supply shock, degenerate under demand-only) and the sensitivity band of
  the headline to `eta_s` per financing column (ETAs §What "estimating
  eta_s" means here).
- No `src/` change, no new generation: existing runs keep their provenance.
- `tests/test_run_manifest.jl` gains the shock keys for new manifests; a
  regression cell with `kind = "none"` must reproduce a committed manifest
  bit-for-bit.

## Enforcement

- A preregistered design per shock id (ADR-0006), registry rows
  `planned -> executed`, manifests + `runs/index.csv` rows, flow table
  `paper/tables/supply_etas_flows.md` citing run ids (ADR-0004).
- The board (`scripts/status.jl`, 0 warnings) and the repository gate
  (`scripts/check_repo.jl`, 0 violations) after every batch.
- `docs/ETAs.md` v3 is the executing annexe and is updated, not re-derived.